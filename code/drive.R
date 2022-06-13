.libPaths('/cluster/huanglab/hhuang/app/envs/R3.6.2/lib/R/library/')

library(nichenetr)
library(tidyverse)
library(Seurat)
library(optparse)

option_list = list(
  make_option(c("-w", "--work_dir"), type="character", default=NULL, 
              help="work directory", metavar="character"),
  
  make_option(c("-d", "--data"), type="character", default=NULL, 
              help="data", metavar="character"),
  
  make_option(c("-s", "--sender"), type="character", default=NULL, 
              help="sender", metavar="character"), 
  
  make_option(c("-r", "--receiver"), type="character", default=NULL, 
              help="receiver", metavar="character"), 
  
  make_option(c("-l", "--ligand_target_matrix"), type="character", default=NULL, 
              help="ligand_target_matrix", metavar="character"), 
  
  make_option(c("-n", "--lr_network"), type="character", default=NULL, 
              help="lr_network", metavar="character"),
  
  make_option(c("-c", "--lfc_cutoff"), type="numeric", default=0.5, 
              help="lfc_cutoff", metavar="numeric"),
  
  make_option(c("-e", "--weighted_networks"), type="character", default=NULL, 
              help="weighted_networks", metavar="character"),
  
  make_option(c("-o", "--condition_colname"), type="character", default=NULL, 
              help="condition_colname", metavar="character"),
  
  make_option(c("-i", "--condition_oi"), type="character", default=NULL, 
              help="condition_oi", metavar="character"),
  
  make_option(c("-f", "--condition_reference"), type="character", default=NULL, 
              help="condition_reference", metavar="character"),
  
  make_option(c("-m", "--organism"), type="character", default=NULL, 
              help="organism", metavar="character"),
  
  make_option(c("-u", "--width"), type="numeric", default=NULL, 
              help="width", metavar="numeric"),
  
  make_option(c("-h", "--height"), type="numeric", default=NULL, 
              help="height", metavar="numeric"),
  make_option(c("-b", "--label"), type="character", default=NULL, 
              help="label", metavar="character")
  
)


opt_parser = OptionParser(option_list=option_list, add_help_option = F);
opt = parse_args(opt_parser);

if (is.null(opt$data)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

setwd(opt$work_dir)
getwd()



# ----------------------------------------
load(opt$lr_network)
load(opt$ligand_target_matrix)
load(opt$weighted_networks)


data = readRDS(opt$data)

Idents(data) <- data@meta.data[, opt$label]

nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = data, lfc_cutoff = opt$lfc_cutoff, 
  receiver = opt$receiver, 
  condition_colname = opt$condition_colname, condition_oi = opt$condition_oi, condition_reference = opt$condition_reference, 
  sender = opt$sender, 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = opt$organism)


pdf(paste0(opt$sender, '_', opt$receiver, '_ligand_activity_target.pdf'), width = opt$width, height = opt$height)
nichenet_output$ligand_target_heatmap
dev.off()

saveRDS(nichenet_output, file = 'nichenet_output.rds')





