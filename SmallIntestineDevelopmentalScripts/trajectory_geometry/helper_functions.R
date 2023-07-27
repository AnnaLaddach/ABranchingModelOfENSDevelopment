devtools::load_all("/home/anna/Documents/TrajectoryGeometry")
library(TrajectoryGeometry)


runMultiPcs = function(nPCs,pseudotime, attributes){
  results = list()
  for (i in 2:nPCs){
    print(i)
    results[[i]] = analyseSingleCellTrajectory(attributes[,1:i], 
                                              pseudotime, 
                                              nSamples = 1000, 
                                              randomizationParams = c('byPermutation','permuteWithinColumns'), 
                                              statistic = "mean", 
                                              N = 1)
  }
  return(results)
}