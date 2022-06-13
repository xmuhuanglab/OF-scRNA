# .libPaths('/cluster/huanglab/hhuang/app/envs/R3.6.2/lib/R/library/')


library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(optparse)
library(Seurat)
library(tidyverse)
library(patchwork)

option_list = list(
  make_option(c("-w", "--work_dir"), type="character", default=NULL, 
              help="work directory", metavar="character"),
  
  make_option(c("-t", "--tumor_data"), type="character", default=NULL, 
              help="tumor dataset", metavar="character"),
  
  make_option(c("-n", "--normal_data"), type="character", default=FALSE, 
              help="normal dataset", metavar="character"), 
  
  make_option(c("-f", "--fetal_data"), type="character", default=NULL, 
              help="fetal database", metavar="character"), 
  
  make_option(c("-g", "--group.by"), type="character", default=NULL, 
              help="group by", metavar="character"), 
  
  make_option(c("-m", "--type"), type="character", default=NULL, 
              help="method to compute communication", metavar="character"),
  
  make_option(c("-r", "--trim"), type="numeric", default=0.1, 
              help="trim", metavar="numeric"),
  
  make_option(c("-q", "--tumor_type"), type="character", default=NULL, 
              help="tumor type", metavar="character"),
  
  make_option(c("-p", "--normal_type"), type="character", default=NULL, 
              help="normal type", metavar="character"),
  
  make_option(c("-e", "--fetal_type"), type="character", default=NULL, 
              help="fetal type", metavar="character"),
  
  make_option(c("-l", "--run_CellChat"), type="character", default=NULL, 
              help="run Cellchat?", metavar="character")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$tumor_data)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

setwd(opt$work_dir)
getwd()

run_cellchat <- function(obj, group.by = NULL, tag = NULL, type = 'truncatedMean', trim = 0.1){
  
  data.input <- GetAssayData(obj, assay = "RNA", slot = "data") # normalized data matrix
  # labels <- Idents(obj)
  # meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
  meta = obj@meta.data
  
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = group.by)
  saveRDS(cellchat, file = paste0(tag, '.cellchat.orig.rds'))
  
  groupSize <- as.numeric(table(cellchat@idents))
  
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  # showDatabaseCategory(CellChatDB)
  
  # use a subset of CellChatDB for cell-cell communication analysis
  # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
  # use all CellChatDB for cell-cell communication analysis
  CellChatDB.use <- CellChatDB.human # simply use the default CellChatDB
  
  # set the used database in the object
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  
  future::plan("multisession", workers = 15) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # project gene expression data onto PPI network (optional)
  
  cellchat <- projectData(cellchat, PPI.human)
  
  cellchat <- computeCommunProb(cellchat, type = type, trim = trim)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  
  cellchat <- aggregateNet(cellchat)
  
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  saveRDS(cellchat, file = paste0(tag, '.cellchat.final.rds'))
  return(cellchat)
}

# opt$normal_data$group.by = gsub(normal_type, tumor_type, opt$normal_data$group.by)
# opt$fetal_data$group.by = gsub(fetal_type, tumor_type, opt$normal_data$group.by)

tumor_data <- readRDS(opt$tumor_data)
tumor_data
fetal_data <- readRDS(opt$fetal_data)
fetal_data
normal_data <- readRDS(opt$normal_data)
normal_data


if (opt$run_CellChat == 'yes') {
  message('Run CellChat')
  tumor.cellchat.final <- run_cellchat(tumor_data, tag = 'tumor', group.by = opt$group.by, type = opt$type, trim = opt$trim)
  normal.cellchat.final <- run_cellchat(normal_data, tag = 'normal', group.by = opt$group.by, type = opt$type, trim = opt$trim)
  fetal.cellchat.final <- run_cellchat(fetal_data, tag = 'fetal', group.by = opt$group.by, type = opt$type, trim = opt$trim)
  
}

tumor.cellchat.final <- readRDS("./tumor.cellchat.final.rds")
normal.cellchat.final <- readRDS("./normal.cellchat.final.rds")
fetal.cellchat.final <- readRDS("./fetal.cellchat.final.rds")

object.list <- list(Tumor = tumor.cellchat.final,  'Adj Normal' = normal.cellchat.final, Fetal = fetal.cellchat.final)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf('outgoing-compare.pdf', height = 5, width = 15)
par(mfrow = c(1,3))
n1 <- netVisual_circle(tumor.cellchat.final@net$count, weight.scale = T
                 , label.edge= T, edge.weight.max = weight.max[2]
                 , edge.width.max = 12
                 , title.name = paste0("Number of interactions in ", 'Tumor')
                 , sources.use = opt$tumor_type, vertex.label.cex = 1, edge.label.cex = 0.8)

n2 <- netVisual_circle(normal.cellchat.final@net$count, weight.scale = T
                 , label.edge= T, edge.weight.max = weight.max[2]
                 , edge.width.max = 12
                 , title.name = paste0("Number of interactions in ", 'Normal')
                 , sources.use = opt$normal_type, vertex.label.cex = 1, edge.label.cex = 0.8)

n3 <- netVisual_circle(fetal.cellchat.final@net$count, weight.scale = T
                 , label.edge= T, edge.weight.max = weight.max[2]
                 , edge.width.max = 12
                 , title.name = paste0("Number of interactions in ", 'Fetal')
                 , sources.use = opt$fetal_type, vertex.label.cex = 1, edge.label.cex = 0.8)

dev.off()


idents.use = c('CD8 T', 'Mac', 'DCs', 'CD4 T', 'NK', 'Endo', 'Hepa')


sub_tumor_cellchat = subsetCellChat(tumor.cellchat.final, idents.use = c(idents.use, opt$tumor_type))

sub_normal_cellchat = subsetCellChat(normal.cellchat.final, idents.use = c(idents.use, opt$normal_type))

sub_fetal_cellchat = subsetCellChat(fetal.cellchat.final, idents.use = c(idents.use, opt$fetal_type))

# ----------- merge objects
df.tumor = subsetCommunication(sub_tumor_cellchat, sources.use = opt$tumor_type, targets.use = idents.use )
df.normal = subsetCommunication(sub_normal_cellchat, sources.use = opt$normal_type, targets.use = idents.use )
df.fetal = subsetCommunication(sub_fetal_cellchat, sources.use = opt$fetal_type, targets.use = idents.use )



LR = setdiff(intersect(paste(df.tumor$interaction_name, df.tumor$target, sep = '|'), paste(df.fetal$interaction_name, df.fetal$target, sep = '|')), paste(df.normal$interaction_name, df.normal$target, sep = '|'))

LR = sapply(strsplit(LR, split = '|', fixed = T), '[[', 1)

pdf('LR-bubble.pdf')
netVisual_bubble(tumor.cellchat.final, sources.use = opt$tumor_type, targets.use =  idents.use
                 , pairLR.use = data.frame(interaction_name = LR))
netVisual_bubble(fetal.cellchat.final, sources.use = opt$fetal_type, targets.use =  idents.use
                 , pairLR.use = data.frame(interaction_name = LR))
dev.off()












