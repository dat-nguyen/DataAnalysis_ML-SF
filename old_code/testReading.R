###############################################################
# preprocessing CASF data for MLSF Project, August 2015
###############################################################

# for splitting training set and test set of CASF data
createTrainingTest <- function(CASFset) {
  #training   = read.table("/home/dat/2007.csv")
  #training   = training[,1]
  #refined   = read.csv("/home/dat/WORK/DB/DESCRIPTORS/CASF/CASF07_refined_SIFt.csv", na.strings=c(".", "NA", "", "?"))
  #write.table(newData  , file = output, sep = ",", row.names = FALSE)
  #filterDataRow(dataPath="/home/dat/WORK/DB/DESCRIPTORS/CASF/CASF07_refined_SIFt.csv", listPath="/home/dat/2007.csv", output="/home/dat/test.csv", nameIndex=14)  
  for (desc in c("_elementsv2-SIFt", "_elementsv2-SIFt-xscore")) {
    refinedPath   = paste(path2Name, CASFset, desc, ".csv", sep = "")
    corePath      = paste(path2Name, substr(CASFset, 0, 6), "_core_PDB.csv", sep = "")
    testPath      = paste(path2Name, substr(CASFset, 0, 6), "_test", desc, ".csv", sep = "")
    trainingPath     = paste(path2Name, substr(CASFset, 0, 6), "_training", desc, ".csv", sep = "")
    #filterDataRow(refinedPath, corePath, output = testPath)
    nameIndex = 1
    selectedList = read.table(corePath)
    selectedList = selectedList[,1]
    data = read.csv(refinedPath, na.strings=c(".", "NA", "", "?"))
    testData      = data[  data[, nameIndex] %in% selectedList ,]    
    trainingData  = data[!(data[, nameIndex] %in% selectedList),]    
    # write to file    
    write.table(testData    , file=testPath, sep=",", quote=FALSE, row.names=FALSE)
    write.table(trainingData, file=trainingPath, sep=",", quote=FALSE, row.names=FALSE)
  }
  
}

for (CASFset in c("CASF07_refined", "CASF12_refined", "CASF13_refined", "CASF14_refined") ) {
  #mergeData(CASFset)
  createTrainingTest(CASFset)
}