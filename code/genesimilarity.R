.libPaths('/cluster/huanglab/hhuang/app/envs/R3.6.2/lib/R/library/')


library(optparse)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)


option_list = list(
  make_option(c("-w", "--work_dir"), type="character", default=NULL, 
              help="work directory", metavar="character"),
  make_option(c("-t", "--tumor_normal_marker"), type="character", default=NULL, 
              help="tumor & normal markers path", metavar="character"),
  make_option(c("-f", "--fetal_marker"), type="character", default=NULL, 
              help="fetal markers path", metavar="character"), 
  make_option(c("-m", "--full_tumor_normal_marker"), type="character", default=NULL, 
              help="full tumor & normal markers path", metavar="character"),
  make_option(c("-n", "--full_fetal_marker"), type="character", default=NULL, 
              help="full fetal markers path", metavar="character"), 
  make_option(c("-q", "--tumor_query"), type="character", default=NULL, 
              help="tumor query", metavar="character"),
  make_option(c("-r", "--fetal_ref"), type="character", default=NULL, 
              help="fetal ref", metavar="character")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



if (is.null(opt$tumor_normal_marker) | is.null(opt$fetal_marker)){
  print_help(opt_parser)
  stop("Please supplied marker files", call.=FALSE)
}

setwd(opt$work_dir)


tn.markers = read.csv(opt$tumor_normal_marker, check.names = F, row.names = 1, stringsAsFactors = F)
fetal.markers = read.csv(opt$fetal_marker, check.names = F, row.names = 1, stringsAsFactors = F)

mtx = matrix(nrow = ncol(tn.markers), ncol = ncol(fetal.markers), dimnames = list(
  rownames = colnames(tn.markers), colnames = colnames(fetal.markers)
))

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


for (i in rownames(mtx)) {
  for (k in colnames(mtx)) {
    
    mtx[i, k] = jaccard( tn.markers[, i], fetal.markers[, k] )
    
  }
}

my.color = colorRampPalette(brewer.pal(9,'Reds'))(256)

pdf(paste0('top', nrow(tn.markers), '-markers-jaccard.pdf'), width = 6, height = 3)
Heatmap(t(mtx), col = my.color, cluster_rows = F, cluster_columns = F, border = T, heatmap_legend_param = list(title = 'Similarity')
        , row_names_gp = gpar(fontsize = 16), column_names_gp = gpar(fontsize = 16) )
dev.off()


full_tn_markers = read.csv(opt$full_tumor_normal_marker, stringsAsFactors = F, check.names = F, row.names = 1)
full_fetal_markers = read.csv(opt$full_fetal_marker, stringsAsFactors = F, check.names = F, row.names = 1)

tn_index = grep(pattern = '_p', colnames(full_tn_markers), ignore.case = T)
fetal_index = grep(pattern = '_p', colnames(full_fetal_markers), ignore.case = T)

for (i in tn_index) {
  full_tn_markers[, i-1][full_tn_markers[, i] >= 0.05 | full_tn_markers[, i+1] <= 1] <- NA
}

for (i in fetal_index) {
  full_fetal_markers[, i-1][full_fetal_markers[, i] >= 0.05 | full_fetal_markers[, i+1] <= 1] <- NA
}

full_tn_markers = full_tn_markers %>% dplyr::select(ends_with('_n'))
full_fetal_markers = full_fetal_markers %>% dplyr::select(ends_with('_n'))

colnames(full_tn_markers) = gsub('_n', '', colnames(full_tn_markers))
colnames(full_fetal_markers) = gsub('_n', '', colnames(full_fetal_markers))


oncofetal = setdiff(setdiff(intersect(full_tn_markers[, opt$tumor_query], full_fetal_markers[, opt$fetal_ref]), unique(unlist(as.list(full_tn_markers[, -which(colnames(full_tn_markers) == opt$tumor_query)]))) ), unique(unlist(as.list(full_fetal_markers[, -which(colnames(full_fetal_markers) == opt$fetal_ref)]))) )

non_oncofetal = setdiff(full_tn_markers$mCAF, c(unique(unlist(as.list(full_fetal_markers))), unique(unlist(as.list(full_tn_markers[, -which(colnames(full_tn_markers) == opt$tumor_query)])))) )[1:length(oncofetal)]


write.csv(oncofetal, file = 'oncofetal.csv')
write.csv(non_oncofetal, file = 'non_oncofetal.csv')

# iCAF = setdiff(tn.markers$iCAF, c(unique(unlist(as.list(fetal.markers))), unique(unlist(as.list(tn.markers[, -2]))) )  )[1:14]
# vCAF = setdiff(tn.markers$vCAF, c(unique(unlist(as.list(fetal.markers))), unique(unlist(as.list(tn.markers[, -1]))) )  )[1:14]
# apCAF = setdiff(tn.markers$apCAF, c(unique(unlist(as.list(fetal.markers))), unique(unlist(as.list(tn.markers[, -5]))) )  )[1:14]









