devtools::load_all("/home/anna/Documents/TrajectoryGeometry")
library(TrajectoryGeometry)
load("neuralPseudoTime.rda")
load("glialPseudoTime.rda")
neuralPseudoTime = neuralPseudoTime[!is.na(neuralPseudoTime)]
glialPseudoTime = glialPseudoTime[!is.na(glialPseudoTime)]
neuralAttributes = readRDS("ensNeuralAttributes.rds")
glialAttributes = readRDS("ensGlialAttributes.rds")
set.seed(42)


for (i in 2:15){
  print(i)
  branchPointResults = analyseBranchPoint(neuralAttributes[,1:i], 
                                          neuralPseudoTime,
                                          randomizationParams = c('byPermutation',
                                                                  'permuteWithinColumns'), 
                                          statistic = "mean",
                                          start = 60,
                                          stop = 100,
                                          step = 5,
                                          nSamples = 1000, 
                                          N = 1)
  
  saveRDS(branchPointResults,paste0("./ensBranchPoint/neuralBranchPointResults_",i,".rds"))
}

for (i in 2:15){
  print(i)
  branchPointResults = analyseBranchPoint(glialAttributes[,1:i], 
                                          glialPseudoTime,
                                          randomizationParams = c('byPermutation',
                                                                  'permuteWithinColumns'), 
                                          statistic = "mean",
                                          start = 60,
                                          stop = 100,
                                          step = 5,
                                          nSamples = 1000, 
                                          N = 1)
  
  saveRDS(branchPointResults,paste0("./ensBranchPoint/glialBranchPointResults_",i,".rds"))
}





