# .libPaths('/cluster/huanglab/hhuang/app/envs/R3.6.2/lib/R/library/')

library(survival)
library(survminer)
library(GSVA)
library(optparse)

option_list = list(
  make_option(c("-w", "--work_dir"), type="character", default=NULL, 
              help="work directory", metavar="character"),
  
  make_option(c("-g", "--genelist"), type="character", default=NULL, 
              help="oncofetal gene list", metavar="character"),
  
  make_option(c("-r", "--data_rna"), type="character", default=NULL, 
              help="rna data", metavar="character"),
  
  make_option(c("-c", "--data_clinical"), type="character", default=NULL, 
              help="clinical data", metavar="character"),
  
  make_option(c("-t", "--data_tag"), type="character", default=NULL, 
              help="dataset name", metavar="character")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$genelist)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

setwd(opt$work_dir)
getwd()

# 
# tn.markers = read.csv(opt$tn_markers, check.names = F, row.names = 1, stringsAsFactors = F)
# fetal.markers = read.csv(opt$fetal_markers, check.names = F, row.names = 1, stringsAsFactors = F)

oncofetal = read.csv(opt$genelist, check.names = F, row.names = 1, stringsAsFactors = F)[, 1]

data <- readRDS(opt$data_rna)
data_surv <- readRDS(opt$data_clinical)



# ---------------------------- gsva
res = gsva(as.matrix(data), list(oncofetal = oncofetal) )
range(res)

all(colnames(res) == data_surv$barcode)
data_surv_rna = cbind(data_surv, t(res) )

dataset = opt$data_tag
data_surv.cut <- surv_cutpoint(
  data_surv_rna,
  time = "time",
  event = "status",
  variables = 'oncofetal'
)


data_surv.cat <- surv_categorize(data_surv.cut)

fit <- survfit(as.formula(paste('Surv(time, status)~', 'oncofetal')), data = data_surv.cat)

cox.res = coxph(as.formula(paste('Surv(time, status)~', 'oncofetal')), data = data_surv.cat)

sum.surv<- summary(cox.res)
c_index <- round(sum.surv$concordance[1], digits = 3)

pval = strsplit(as.character(surv_pvalue(fit)[4]), split = ' = ')[[1]][2]
# pval = round(surv_pvalue(fit)[2], 5)

g <- ggsurvplot(fit, risk.table = F, pval = F, conf.int = F, surv.median.line = 'hv', palette = 'Set1', pval.size = 7
                , ggtheme = theme_survminer(font.main = c(25, "plain", "black"),
                                            font.submain = c(15, "plain", "black"),
                                            font.x = c(18, "plain", "black"),
                                            font.y = c(18, "plain", "black"),
                                            font.caption = c(15, "plain", "black"),
                                            font.tickslab = c(15, "plain", "black"),
                                            legend = c("right"),
                                            font.legend = c(19, "plain", "black"))
                # ,pval.coord = c(0, 0.2)
)

g$plot <- g$plot + annotate('text', label = paste0('n(High) = ', fit$n[1]),  size = 7, x = 15, y = 0.4) +
  annotate('text', label = paste0('n(Low) = ', fit$n[2]), x = 15, y = 0.3, size = 7) + xlab(label = 'Time(months)')+
  ylab(label = 'Overall survival') +
  annotate('text', label = paste0('C-Index = ', c_index), x = 15, y = 0.1, size = 7) + ggtitle(label = dataset) + 
  annotate('text', label = paste0('P = ', pval), x = 15 , y = 0.2, size = 7)

pdf(paste0('oncofetal', '-gsva-', dataset, '-KM-plot.pdf'), onefile = F, width = 9, height = 5)
print(g)
dev.off()

