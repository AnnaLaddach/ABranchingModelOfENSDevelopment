
opt.clust = function(dat, min.clusters, max.clusters, output){
  #based on a function written by Jamie MacPherson
  # read data
  print(dat)
  # determine optimal number of clusters
  nc = NbClust(dat,
               min.nc=min.clusters,
               max.nc=max.clusters,
               method="kmeans")
  
  # bar plot of optimal clusters
  nc.plt = barplot(table(nc$Best.n[1,]),
                   xlab="Numer of Clusters",
                   ylab="Number of Criteria",
                   colour = 'lightblue')
  
  # write barplot to pdf file
  pdf(output)
  nc.plt = barplot(table(nc$Best.n[1,]),
                   xlab="Number of Clusters",
                   ylab="Number of Criteria",
                   colour = 'lightblue')
  print(nc.plt)
  dev.off()
}