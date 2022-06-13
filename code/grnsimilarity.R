# .libPaths('/cluster/huanglab/hhuang/app/envs/R4/lib/R/library/')
# .libPaths('/cluster/huanglab/hhuang/app/envs/R3.6.2/lib/R/library/')

library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(SCENIC)
library(AUCell)
library(RcisTarget)
# library(SCopeLoomR)
library(KernSmooth)
library(BiocParallel)
library(ggplot2)
library(data.table)
library(grid)
library(ComplexHeatmap)

options(stringsAsFactors = F)

library(optparse)
library(Seurat)
library(tidyverse)


option_list = list(
  make_option(c("-w", "--work_dir"), type="character", default=NULL, 
              help="work directory", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input dataset", metavar="character"),
  make_option(c("-l", "--run_SCENIC"), type="character", default=FALSE, 
              help="run SCENIC", metavar="character"), 
  make_option(c("-d", "--database"), type="character", default=NULL, 
              help="Rcistarget database", metavar="character"), 
  make_option(c("-n", "--n_cores"), type="numeric", default=10, 
              help="Cores to run", metavar="numeric"),
  make_option(c("-c", "--minCountsPerGene_rate"), type="numeric", default=0.001, 
              help="Filtering gene by min counts", metavar="numeric"), 
  make_option(c("-s", "--minSamples_rate"), type="numeric", default=0.01, 
              help="Filtering gene by min samples", metavar="numeric"),
  make_option(c("-q", "--tumor_query"), type="character", default=NULL, 
              help="tumor query", metavar="character"),
  make_option(c("-r", "--fetal_ref"), type="character", default=NULL, 
              help="fetal ref", metavar="character"),
  make_option(c("-y", "--width"), type="numeric", default=8, 
              help="width", metavar="numeric"),
  make_option(c("-h", "--height"), type="numeric", default=10, 
              help="height", metavar="numeric")
)


opt_parser = OptionParser(option_list=option_list, add_help_option = F);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

setwd(opt$work_dir)


obj = readRDS(opt$input)

if (opt$run_SCENIC == 'yes') {
  
  
  # rm(mtx)
  exprMat = GetAssayData(obj, slot = 'data', assay = 'RNA') %>% as.matrix()
  cellInfo = obj@meta.data
  cellInfo$subtype
  dim(exprMat)
  
  dir.create("int")
  saveRDS(cellInfo, file="int/cellInfo.Rds")
  
  # Color to assign to the variables (same format as for NMF::aheatmap)
  colVars <- list(subtype = c(colorRampPalette(brewer.pal(9, 'Set1'))(length(names(table(cellInfo$subtype))))))
  names(colVars$subtype) = names(table(cellInfo$subtype))
  
  saveRDS(colVars, file="int/colVars.Rds")
  # plot.new(); legend(0,1, fill=colVars$subtype, legend=names(colVars$subtype))
  
  library(SCENIC)
  org <- "hgnc" # or hgnc, or dmel
  dbDir <- opt$databse # RcisTarget databases location
  # myDatasetTitle <- "Fibro" # choose a name for your analysis
  # data(defaultDbNames)
  # dbs <- defaultDbNames[[org]]
  scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=list.files(dbDir, pattern = 'hg38'), datasetTitle=myDatasetTitle, nCores=10) 
  
  # Modify if needed
  scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
  scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
  # Databases:
  # scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr"="mm9-tss-centered-5kb-10species.mc8nr.feather")
  # scenicOptions@settings$db_mcVersion <- "v8"
  
  # Save to use at a later time...
  saveRDS(scenicOptions, file="int/scenicOptions.Rds")
  
  # (Adjust minimum values according to your dataset)
  genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                             minCountsPerGene=3*opt$minCountsPerGene_rate*ncol(exprMat),
                             minSamples=ncol(exprMat)*opt$minSamples_rate)
  
  # load("/cluster/huanglab/hhuang/project/Cancer_dev/Liver/Results/scanpy/Compare/Trajectory/trajectory_gene/final-gene.Rdata")
  # interestingGenes <- final.genes
  # # any missing?
  # interestingGenes[which(!interestingGenes %in% genesKept)]
  
  exprMat_filtered <- exprMat[genesKept, ]
  save(exprMat_filtered, file = 'exprMat_filtered.Rdata')
  dim(exprMat_filtered)
  rm(exprMat)
  
  runCorrelation(exprMat_filtered, scenicOptions)
  
  
  # Run GENIE3
  runGenie3(exprMat_filtered, scenicOptions)
  ### Build and score the GRN
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
  
  # Optional: Binarize activity
  # aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
  # savedSelections <- shiny::runApp(aucellApp)
  # newThresholds <- savedSelections$thresholds
  # scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
  # saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
  scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
  # tsneAUC(scenicOptions, aucType="AUC") # choose settings
  
  # Export:
  # saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo")) # Temporary, to add to loom
  # export2loom(scenicOptions, exprMat_filtered)
  
  saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
}

# Correlation
cellInfo <- readRDS("./int/cellInfo.Rds")
# cellInfo = obj@meta.data
scenicOptions <- readRDS("./int/scenicOptions.Rds")

minPerc <- .5
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$subtype),
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]

binary_cor = cor(binaryActPerc_subset)

my.color = colorRampPalette(brewer.pal(9, 'RdBu'))(3)

library(ggcorrplot)
pdf('./correlation-heatmap.pdf')
ggcorrplot(binary_cor,  hc.order = F, type = "lower", tl.cex = 18, 
           outline.col = "white", method = 'square', show.diag = F
           , colors = c("#6D9EC1", "white", "#E46726"), lab = F) + 
  theme( axis.text = element_text(face = 'bold', size = 15))
dev.off()
# binaryRegulonActivity <- readRDS("./int/4.1_binaryRegulonActivity.Rds")
# binaryRegulonActivity_nonDupl <- readRDS("./int/4.2_binaryRegulonActivity_nonDupl.Rds")
# regulonAUC <- readRDS("./int/3.4_regulonAUC.Rds")



# RSS
binaryRegulonActivity_nonDupl <- readRDS("./int/4.2_binaryRegulonActivity_nonDupl.Rds")
rss <- calcRSS(AUC=binaryRegulonActivity_nonDupl, cellAnnotation=cellInfo[colnames(binaryRegulonActivity_nonDupl), "subtype"])
PlotRSS <- function(rss, zThreshold=1, varName="cellType"){
  rssNorm <- scale(rss)    
  rssNorm[rssNorm < 0] <- 0    
  rss.df <- reshape2::melt(rss)    
  colnames(rss.df) <- c("Topic", varName, "RSS")    
  rssNorm.df <- reshape2::melt(rssNorm)    
  colnames(rssNorm.df) <- c("Topic", varName, "Z") 
  rss.df <- base::merge(rss.df, rssNorm.df) 
  return(rss.df)
}

# write.csv(rss, file = 'rss.csv')
rss_df <- PlotRSS(rss)
rss_df = rss_df[rss_df$Z > 1, ]
write.csv(rss_df, file = 'rss-table.csv')

table(rss_df$cellType)

before_name = names(table(rss_df$cellType))

rss_df$RSS[rss_df$RSS < 0.1] <- 0

a = rss_df %>% dplyr::filter(RSS > 0.1 & (cellType == opt$tumor_query | cellType == opt$fetal_ref) ) %>% pull(Topic)

rss_df = rss_df[rss_df$Topic %in% a, ]
table(rss_df$cellType)
write.csv(rss_df, file = 'tmp.csv')

filter_name = names(table(rss_df$cellType))

ll = setdiff(before_name, filter_name)

aa = lapply(ll, function(x) {
  x = data.frame(Topic = rep(a, 1), 
                 cellType = c( rep(x, length(a) ))
                 , RSS = c(rep(0, length(a))),
                 Z = c(rep(0, length(a)) ) )
  
  
})

rss_df = bind_rows(aa, rss_df)


# rss_df = rbind(rss_df, data.frame(Topic = rep(a, 2), 
#                                   cellType = c( rep('apFibro', length(a) )
#                                                 , rep('vCAF', length(a) ) )
#                                   , RSS = c(rep(0, length(a))
#                                             , rep(0, length(a))), 
#                                   Z = c(rep(0, length(a))
#                                         , rep(0, length(a))) ) )


rss_df$cellType = factor(rss_df$cellType, levels = c(opt$tumor_query, opt$fetal_ref, setdiff(before_name, c(opt$tumor_query, opt$fetal_ref))))

colnames(rss_df)[4] <- 'Z-score'
pdf('RSS-filtering.pdf', height = opt$height, width = opt$width)
ggplot(rss_df, mapping = aes(x = cellType, y = Topic)) + geom_point(mapping = aes(color = `Z-score`, size = RSS)) + 
  scale_color_gradient2() + theme_bw(base_size = 15) + scale_size(breaks = c(0.1, 0.2, 0.3, 0.4)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text = element_text(face = 'bold', size = 18)) + ylab(label = '') + xlab('')
dev.off()











