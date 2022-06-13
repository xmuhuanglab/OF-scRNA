# .libPaths('/cluster/huanglab/hhuang/app/envs/R3.6.2/lib/R/library/')


library(optparse)
library(Seurat)
library(tidyverse)
library(simspec)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input dataset file directory", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="tnf_css.h5ad", 
              help="output file name [default= %default]", metavar="character"), 
  make_option(c("-w", "--workdir"), type="character", default=NULL, 
              help="work directory, save results", metavar="character"), 
  make_option(c("-b", "--batch"), type="character", default=NULL, 
              help="batch tag", metavar="character"),
  make_option(c("-c", "--css_cluster_resolution"), type="numeric", default=0.6, 
              help="cluster resolution to css", metavar="numeric"),
  make_option(c("-s", "--seurat_resolution"), type="numeric", default=0.6, 
              help="cluster resolution to seurat clustering", metavar="numeric"), 
  make_option(c("-k", "--k_neighbor"), type="numeric", default=20, 
              help="k neighbor for clustering", metavar="numeric")
)




opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}



library(sceasy)
library(reticulate)
library(Seurat)


# reticulate::use_condaenv('/cluster/huanglab/hhuang/app/envs/scvelo/', required = T)

loompy <- reticulate::import('loompy')

anndata <- reticulate::import('anndata')

# dir = '/cluster/huanglab/hhuang/project/Cancer_dev/Liver/OF-scRNA/results/1.Identification/scanpy2seurat'
# file.path(dir, 'X.csv')

ReadScanpy = function(dir, hvg = T, umap = F, pca =T, css = F, diffmap = F) {
  
  # matrix
  mat = readr::read_csv(file.path(dir, 'X.csv'), col_names = F)
  mat = t(mat)
  
  meta = read.csv(file.path(dir, 'obs.csv'), row.names = 1)
  
  # meta = meta %>% as.data.frame() %>% column_to_rownames(colnames(meta)[1])
  
  var = read.csv(file.path(dir, 'var.csv'), row.names = 1)
  
  colnames(mat) <- rownames(meta)
  rownames(mat) <- rownames(var)
  
  obj = CreateSeuratObject(counts = mat, meta.data = meta)
  
  # 
  if (hvg == T) {
    VariableFeatures(obj) <- rownames(var)[var$highly_variable == 'True']
  }
  
  dimen <- readr::read_csv(file.path(dir, 'obsm.csv'))
  
  if (pca == T) {
    pca = dimen %>% dplyr::select(starts_with('X_pca'))
    pca = as.data.frame(pca)
    rownames(pca) <- colnames(mat)
    obj[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(pca), key = "PC_", assay = DefaultAssay(obj))
  }
  
  if (umap == T) {
    umap = dimen %>% dplyr::select(X_umap1, X_umap2)
    umap = as.data.frame(umap)
    rownames(umap) <- colnames(mat)
    obj[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(obj))
    
  }
  
  if (diffmap == T) {
    diffmap = dimen %>% dplyr::select(X_diffmap1, X_diffmap2)
    diffmap = as.data.frame(diffmap)
    rownames(diffmap) <- colnames(mat)
    obj[["diffmap"]] <- CreateDimReducObject(embeddings = as.matrix(diffmap), key = "diffmap_", assay = DefaultAssay(obj))
    
  }
  
  if (css == T) {
    css = dimen %>% dplyr::select(starts_with('X_css'))
    css = as.data.frame(css)
    rownames(css) <- colnames(mat)
    obj[["css"]] <- CreateDimReducObject(embeddings = as.matrix(css), key = "CSS_", assay = DefaultAssay(obj))
  }
  
  obj <- ScaleData(obj)
  return(obj)
}


setwd(opt$workdir)
getwd()

# fibro = ReadScanpy(dir = opt$input)

# load("/cluster/huanglab/hhuang/project/Cancer_dev/Liver/OF-scRNA/data/fibro.Rdata")

runCSS <- function(dir = opt$input, label_tag = NULL, css_cluster_resolution = 0.6, seurat_resolution = 0.55, k = 20
                   , plot_by = label_tag, outFile = opt$output){
  message("Read Scanpy data")
  obj <- ReadScanpy(dir = dir)
  
  message('Run Cluster similarity spectrum algorithm to correct batch effect')
  
  if (is.null(label_tag) == T) {
    stop("Please provide label_tag")
  }

  obj <- cluster_sim_spectrum(obj, label_tag = label_tag, cluster_resolution = css_cluster_resolution, k = k)
  obj <- RunUMAP(obj, reduction = "css", dims = 1:ncol(Embeddings(obj, "css")), seed.use = 43)
  obj <- FindNeighbors(obj, reduction = "css", dims = 1:ncol(Embeddings(obj, "css")))
  obj <- FindClusters(obj, resolution = seurat_resolution)
  
  
  pdf('UMAP-figure.pdf')
  
  if (is.null(plot_by) == T) {
    
  }else {
    p <- lapply(plot_by, function(x) { DimPlot(obj, group.by = x) })
    print(p)
  }
  dev.off()
  
  sceasy::convertFormat(obj, from="seurat", to="anndata",
                      outFile=opt$output)
  save(obj, file = 'Fibro_css.Rdata')
  return(obj)
}


fibro = runCSS(opt$input, label_tag = opt$batch, css_cluster_resolution = opt$css_cluster_resolution, seurat_resolution = opt$seurat_resolution, plot_by = c('batch', 'seurat_clusters'), k = opt$k_neighbor)



















