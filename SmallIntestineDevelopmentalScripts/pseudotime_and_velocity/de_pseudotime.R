library(monocle)
library(gplots)
library(magrittr)


de_pseudotime = function(adata, pseudotime, name, all_genes = T, selected_genes = c()){
  
  ## normalise pseudotime to range 100
  pseudotime_norm = pseudotime %>% {100*((. - min(.))/(max(.) - min(.)))}
  
  ## format data for monocle
  mat = counts(adata)
  mat = mat[rowSums(mat > 0) > 9,]
  print(dim(mat))
  anno = data.frame(colData(adata))
  feature = data.frame(rownames(mat))
  rownames(feature) = rownames(mat)
  names(feature) = 'gene_short_name'
  anno$Pseudotime = pseudotime_norm  
  
  ## create monocle dataset
  currHSMM <- monocle::newCellDataSet(mat,
                                      phenoData = AnnotatedDataFrame(anno), 
                                      featureData = AnnotatedDataFrame(feature),
                                      lowerDetectionLimit = 0.5,
                                      expressionFamily = negbinomial.size())
  currHSMM <- estimateSizeFactors(currHSMM)
  currHSMM <- estimateDispersions(currHSMM,remove_outliers = F) 
  cds_subset = currHSMM
  
  ## choose whether to run only on highly variable genes
  if (all_genes == F){
    cds_subset = cds_subset[rowData(adata)$highly_variable,]
  }
  if (length(selected_genes) > 0){
    cds_subset = cds_subset[selected_genes,]
  }
  
  new_data = data.frame(Pseudotime = seq(min(pseudotime_norm), max(pseudotime_norm), length.out = 100)) 
  sc = monocle::genSmoothCurves(cds_subset, cores=1, trend_formula = "~sm.ns(Pseudotime, df = 3)", relative_expr = T, new_data = new_data)
  sc = na.omit(sc)
  write.csv(sc, paste0(name,"_smooth_curves_pseudotime.csv"), quote = F)
  
  ## normalise curves to range from 0 to 1
  sc_norm = sc/rowMeans(sc)
  sc_norm = (sc- rowMin(sc))/(rowMax(sc)-rowMin(sc)) 
  
  ## differential expression testing over pseudotime
  diff_test_res <- differentialGeneTest(cds_subset,
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
  diff_test_res = diff_test_res[order(diff_test_res$qval),]
  write.csv(diff_test_res,paste0(name,"_de_pseudotime.csv"), quote = F)
}

plot_smooth_curves = function(sc, genes, name, rowv = T, file = F){
  if (file){
    sc = read.table(sc, sep = ",")
    rownames(sc) = sc$V1
    sc$V1 = NULL
  }
  sc = as.matrix(sc)
  genes = genes[genes %in% rownames(sc)]
  sc = sc[genes,]
  sc_norm = (sc- rowMin(sc))/(rowMax(sc)-rowMin(sc)) 
  png(paste0(name,"_heatmap_pseudotime.png"), height = 10, width = 4, units = "in", res= 300)
  breaks = seq(0,1,length = 101)
  color.palette = colorRampPalette(c("blue", "white", "red"))
  h = heatmap.2(sc_norm, Colv = F, Rowv = rowv, trace = 'none',col = color.palette, breaks = breaks,
            cexRow = 1,labCol = F)

  dev.off()
  return(rownames(sc)[h$rowInd])
}
