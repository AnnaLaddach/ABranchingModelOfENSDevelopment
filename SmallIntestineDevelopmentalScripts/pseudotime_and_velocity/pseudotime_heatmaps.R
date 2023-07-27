## load libraries
library(slingshot)
library(RColorBrewer)
library(SingleCellExperiment)
library(ggplot2)
library(monocle)
library(gplots)
source("de_pseudotime.R")


get_name = function(name){
  try(return(strsplit(name,'\\.')[[1]][1]))
}

## function to create heatmap for differentially expressed genes over pseudotime
plot_genes_over_pt = function(cell_type, stage, tf = T){
  sc = paste0("../../output/DE_pseudotime/GM/",cell_type, stage,"_10_GM_smooth_curves_pseudotime.csv")
  genes = read.csv(paste0("../../output/DE_pseudotime/GM/",cell_type, stage,"_10_GM_de_pseudotime.csv"),stringsAsFactors = F)
  genes$name = as.character(lapply(as.character(genes$gene_short_name),get_name))
  genes = genes[genes$qval < 0.05,]
  tf_lab = ""
  if (tf){
    genes = genes[genes$name %in% tfs$gene_name,]
    tf_lab = "_tf"
  }
  genes = genes$gene_short_name
  n_genes = 30
  if (length(genes) < 30){
    n_genes = length(genes)
  }
  genes = genes[1:n_genes]
  gene_order = plot_smooth_curves(sc, genes, paste0("../../output/DE_pseudotime/GM/top_30/",cell_type, stage, tf_lab, "_10_GM_top_30"), file = T)
}

tfs = read.table("/home/anna/software/Antler-master/inst/extdata/Annotations/TF_mmusculus_150722.dat",stringsAsFactors = F)

plot_genes_over_pt("glial", "_emb",F)
plot_genes_over_pt("neural", "_emb",F)
plot_genes_over_pt("neural_bp_85", "_emb",F)
plot_genes_over_pt("neural_early", "_emb",F)
plot_genes_over_pt("neural_early_bp", "_emb",F)

plot_genes_over_pt("glial", "_all", T)
plot_genes_over_pt("neural", "_all",T)
plot_genes_over_pt("neural_bp", "_all",T)
plot_genes_over_pt("neural_early", "_all",T)
plot_genes_over_pt("neural_early_bp", "_all",T)



plot_genes_over_pt("glial", "_emb",T)
plot_genes_over_pt("neural", "_emb",T)
plot_genes_over_pt("neural_bp_85", "_emb",T)
plot_genes_over_pt("neural_early", "_emb",T)
plot_genes_over_pt("neural_early_bp", "_emb",T)

plot_genes_over_pt("glial", "_all")
plot_genes_over_pt("neural", "_all")
plot_genes_over_pt("neural_bp", "_all")
plot_genes_over_pt("neural_early", "_all")
plot_genes_over_pt("neural_early_bp", "_all")


