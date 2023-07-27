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
adata = adata[rowData(adata)["highly_variable"]$highly_variable,]

##subset for different trajectories
#trajectory from early progenitor to adult neuron
adata_neural = adata[,!is.na(adata$slingPseudotime_1) & slingCurveWeights(SlingshotDataSet(adata))[,1] == 1]
adata_neural = adata_neural[,order(adata_neural$slingPseudotime_1)]

#trajectory from early progenitor to adult glia
adata_glial = adata[,!is.na(adata$slingPseudotime_2) & slingCurveWeights(SlingshotDataSet(adata))[,2] == 1]
adata_glial = adata_glial[,order(adata_glial$slingPseudotime_2)]


#trajectory from early progenitor to early neuron
adata_neural_early = adata[,!is.na(adata$slingPseudotime_3) & slingCurveWeights(SlingshotDataSet(adata))[,3] == 1]
adata_neural_early = adata_neural_early[,order(adata_neural_early$slingPseudotime_3)]

##run DE for trajectories over all cells
# de_pseudotime(adata_neural, adata_neural$slingPseudotime_1, "../../output/DE_pseudotime/neural_all", all_genes = T)
# de_pseudotime(adata_neural[,adata_neural$slingPseudotime_1 > 85], adata_neural[,adata_neural$slingPseudotime_1 > 85]$slingPseudotime_1,
#              "../../output/DE_pseudotime/neural_bp_85_all", all_genes = T) #neural trajectory from branch point
# de_pseudotime(adata_glial, adata_glial$slingPseudotime_2, "../../output/DE_pseudotime/glial_all", all_genes = T)
# de_pseudotime(adata_neural_early, adata_neural_early$slingPseudotime_3, "../../output/DE_pseudotime/neural_early_all", all_genes = T)
# de_pseudotime(adata_neural_early[,adata_neural_early$slingPseudotime_3 > 38], adata_neural_early[,adata_neural_early$slingPseudotime_3 > 38]$slingPseudotime_3,
#             "../../output/DE_pseudotime/neural_early_bp_all", all_genes = T) #early neural trajectory from branch point

# ##run DE for trajectories over embryonic cells
# de_pseudotime(adata_neural[,adata_neural$slingPseudotime_1 > 85 & !(adata_neural$time_point %in% c("P25","P60"))],
#             adata_neural[,adata_neural$slingPseudotime_1 > 85  & !(adata_neural$time_point %in% c("P25","P60"))]$slingPseudotime_1,
#             "../../output/DE_pseudotime/GM/neural_bp_85_emb_10_GM_stringent", all_genes = T) #neural trajectory from branch point
# de_pseudotime(adata_neural[,!(adata_neural$time_point %in% c("P25","P60"))], adata_neural[,!(adata_neural$time_point %in% c("P25","P60"))]$slingPseudotime_1, "../../output/DE_pseudotime/GM/neural_emb_10_GM_stringent", all_genes = T)
de_pseudotime(adata_glial[,!(adata_glial$time_point %in% c("P25","P60"))], adata_glial[,!(adata_glial$time_point %in% c("P25","P60"))]$slingPseudotime_2, "../../output/DE_pseudotime/GM/glial_emb_10_GM_stringent", all_genes = T)
de_pseudotime(adata_neural_early[,!(adata_neural_early$time_point %in% c("P25","P60"))],
adata_neural_early[,!(adata_neural_early$time_point %in% c("P25","P60"))]$slingPseudotime_3,
"../../output/DE_pseudotime/GM/neural_early_emb_10_GM_stringent", all_genes = T)
# de_pseudotime(adata_neural_early[,adata_neural_early$slingPseudotime_3 > 38 & !(adata_neural_early$time_point %in% c("P25","P60"))],
#               adata_neural_early[,adata_neural_early$slingPseudotime_3 > 38 & !(adata_neural_early$time_point %in% c("P25","P60"))]$slingPseudotime_3,
#               "../../output/DE_pseudotime/neural_early_bp_emb", all_genes = T) #early neural trajectory from branch point


