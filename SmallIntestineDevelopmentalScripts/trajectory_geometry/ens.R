source("helper_functions.R")
set.seed(42)
ensNeuralPseudotime = readRDS("ensNeuralPseudotimeNormalised.rds")
ensGlialPseudotime = readRDS("ensGlialPseudotimeNormalised.rds")
ensNeuralAttributes = readRDS("ensNeuralAttributes.rds")
ensGlialAttributes = readRDS("ensGlialAttributes.rds")
ensNeuralResults = runMultiPcs(15, ensNeuralPseudotime, ensNeuralAttributes)
saveRDS(ensNeuralResults,"ensNeuralResults.rds")
ensGlialResults = runMultiPcs(15, ensGlialPseudotime, ensGlialAttributes)
saveRDS(ensGlialResults,"ensGlialResults.rds")
