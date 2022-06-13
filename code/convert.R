# .libPaths('/cluster/huanglab/hhuang/app/envs/R3.6.2/lib/R/library/')

library(sceasy)
library(reticulate)
library(Seurat)
library(optparse)
library(tidyverse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="input dataset file directory", metavar="character"),
  make_option(c("-c", "--celltype"), type="character", default=NULL, 
              help="celltype", metavar="character"), 
  make_option(c("-w", "--work_dir"), type="character", default=NULL, 
              help="work dir", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


setwd(opt$work_dir)
getwd()


# reticulate::use_condaenv('/cluster/huanglab/hhuang/app/envs/scvelo/', required = T)

loompy <- reticulate::import('loompy')

anndata <- reticulate::import('anndata')

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


obj <- ReadScanpy(dir = opt$file, hvg = T, umap = T, pca = T, css = T, diffmap = F)

saveRDS(obj, file = paste0(opt$celltype, '_final.rds'))

