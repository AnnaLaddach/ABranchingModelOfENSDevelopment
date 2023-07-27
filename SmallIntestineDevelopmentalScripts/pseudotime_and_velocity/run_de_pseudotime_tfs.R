#load libraries
library(slingshot)
library(RColorBrewer)
library(SingleCellExperiment)
library(ggplot2)
library(monocle)
library(gplots)
source("de_pseudotime.R")

##read trajectory
adata = readRDS('../../output/slingshot/adata_start_louvain_semisup.rds')
tfs = read.table("/home/anna/software/Antler-master/inst/extdata/Annotations/TF_mmusculus_150722.dat",stringsAsFactors = F)
tfs = tfs$gene_name

get_name = function(name){
  try(return(strsplit(name,'\\.')[[1]][1]))
}
gene_names = as.character(lapply(as.character(rownames(adata)),get_name))

tfs_selected = rownames(adata)[gene_names %in% tfs]
adata = adata[tfs_selected,]
##subset for different trajectories
#trajectory from early progenitor to adult neuron
adata_neural = adata[,!is.na(adata$slingPseudotime_1)]
adata_neural = adata_neural[,order(adata_neural$slingPseudotime_1)]


#trajectory from early progenitor to early neuron
adata_neural_early = adata[,!is.na(adata$slingPseudotime_3)]
adata_neural_early = adata_neural_early[,order(adata_neural_early$slingPseudotime_3)]


##run DE for trajectories over embryonic cells
de_pseudotime(adata_neural[,adata_neural$slingPseudotime_1 > 85 & !(adata_neural$time_point %in% c("P25","P60"))],
             adata_neural[,adata_neural$slingPseudotime_1 > 85  & !(adata_neural$time_point %in% c("P25","P60"))]$slingPseudotime_1,
            "../../output/DE_pseudotime/TF/neural_bp_85_emb_10_", all_genes = T) #neural trajectory from branch point
# de_pseudotime(adata_neural[,!(adata_neural$time_point %in% c("P25","P60"))], adata_neural[,!(adata_neural$time_point %in% c("P25","P60"))]$slingPseudotime_1, "../../output/DE_pseudotime/TF/neural_emb_10_", all_genes = T)
# de_pseudotime(adata_neural_early[,!(adata_neural_early$time_point %in% c("P25","P60"))],
# adata_neural_early[,!(adata_neural_early$time_point %in% c("P25","P60"))]$slingPseudotime_3,
# "../../output/DE_pseudotime/TF/neural_early_emb_10_", all_genes = T)
de_pseudotime(adata_neural_early[,adata_neural_early$slingPseudotime_3 > 60 & !(adata_neural_early$time_point %in% c("P25","P60"))],
              adata_neural_early[,adata_neural_early$slingPseudotime_3 > 60 & !(adata_neural_early$time_point %in% c("P25","P60"))]$slingPseudotime_3,
              "../../output/DE_pseudotime/TF/neural_early_bp_emb_10", all_genes = T) #early neural trajectory from branch point


