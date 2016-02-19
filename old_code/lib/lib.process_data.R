#####################################################################
# library for manipulating/processing data
#####################################################################

# make all values of the col column from data unique
# data as dataframe
makeUniqueList <- function(data, col = 1) {
  # rounding the value
  data[,col] = round(data[,col],6)
  for (i in seq(length(data[,col]))) {
    if (sum(data[i, col] == data[, col]) > 1) { # if exist more than one entry with the same value
      data[i, col] = data[i, col] + 0.000001 # increase the original valueby 0.000001 (could be dealt with 1mil entries)
    }
  }
  return (data)
}

# check two data frames with a given index whether two data frames have common rows
checkCommonRows <- function(dataA, dataB, index) {
  if (length(intersect(dataA[,index], dataB[,index])) > 0) {
    return (TRUE)
  }
  else return (FALSE)
}

# create complements from 2 data frames with a given index row (= dataA - dataB). Pay attention to the order of two sets.
createRowsComplementsData <- function(dataA, dataB, index) {
  newData = dataA[!(dataA[,index] %in% dataB[,index]),]
  return (newData)
}

# create a new dataset from a given dataset and a list of selected rows
filterDataRow <- function(dataPath, listPath, output, nameIndex=1) {
  selectedList = read.table(listPath)
  selectedList = selectedList[,1]
  data = read.csv(dataPath)
  newData  = data[data[, nameIndex] %in% selectedList,]
  
  write.table(newData, file=output, sep=",", quote=FALSE, row.names=FALSE)
}
  
# read validation data from cvs file at outputDir, if another file exists then row binding two datas 
readAndCombiValiData <- function(outputDir = "/home/dat/WORK/output/", dataFile1, dataFile2="") {
  evalData 	= read.csv(paste(outputDir,dataFile1,".csv", sep=""))
  if (dataFile2 != "") { # if another dataFile2 exists, then combine rows of the 2 datasets
    evalData2 	= read.csv(paste(outputDir,dataFile2,".csv", sep=""))
    evalData	= rbind(evalData, evalData2)
    # \TODO: (fix) quick and dirty to remove FXa from CSAR set
    FXa = read.table("/home/dat/WORK/output/FXa.txt")
    FXa = FXa[,1]
    evalData	= evalData[!evalData[, 1] %in% FXa,]
  }
  return (evalData)
}

# merge validation data after column
mergeScoresData <- function(outputDir = "/home/dat/WORK/output/", dataSet1, dataFile2) {	
  for (method in METHODS) {
    dataSet2 	= read.csv(paste(outputDir,"ML-Scores/",dataFile2,method,".csv", sep=""))
    mergeData	= merge(dataSet1, dataSet2[,1:2], by.x=1 , by.y=1)
    colnames(mergeData)[length(mergeData[1,])] = method
    dataSet1 	= mergeData
  }
  names(mergeData)[1] = "PDB"
  write.csv(mergeData, file=paste(outputDir,dataFile2,"allScores.csv",sep=""),row.names = FALSE)
  return (mergeData)
}

# create full desc depends on cutoff and binsize value (only matters for desc=credo AND binsize=12)
concatCutoffAndBinsize <- function(desc, cutoff=12, binsize=0) {
  if (desc == "credo" && binsize == 12) {
      fullDesc = paste(desc,"_c",cutoff, sep="")
  }
  else {
      fullDesc = paste(desc,"_c",cutoff,"b",binsize, sep="")
  }
  return (fullDesc)
}

# concatenate path for reading descriptors
#concatPathWithDataName <- function(path, prefix, dataName, desc, cutoff, binsize) {
#  fullPath = paste(path,prefix,dataName,"_",concatCutoffAndBinsize(desc, cutoff, binsize),".csv", sep="")
#  return (fullPath)
#}

# concatenate path for reading descriptors
concatPath <- function(path, prefix, cutoff, binsize, desc, postfix="") {
  fullPath = paste(path,prefix,concatCutoffAndBinsize(desc, cutoff, binsize),postfix,".csv", sep="")
  return (fullPath)
}

# create training and test data sets from checkCommonRows of 2 given sets
# the order of the given sets is IMPORTANT
# 1. set = smaller set, which is typically test set
# 2. set = bigger set, which is typically base for training set
# USAGE: 
#source("lib/lib.process_data.R")
#createTrainingTestData(path = "/home/dat/WORK/dev/data/CASF12_raw/", prefix = "CASF12_", desc = "elements", binsize = 1)
#createTrainingTestData(path = "/home/dat/WORK/dev/data/CASF12_raw/", prefix = "CASF12_", desc = "sybyl", binsize = 1)
#createTrainingTestData(path = "/home/dat/WORK/dev/data/CASF12_raw/", prefix = "CASF12_", desc = "credo", binsize = 1)
#
createTrainingTestData <- function(path, prefix, dataName = c("core", "refined"), desc = "elements", cutoff = 12, binsize = 0, nameIndex) {
  testData 	  = read.csv(paste(path,prefix,dataName[1],"_",concatCutoffAndBinsize(desc, cutoff, binsize),".csv", sep=""), na.strings=c(".", "NA", "", "?"))
  trainData 	= read.csv(paste(path,prefix,dataName[2],"_",concatCutoffAndBinsize(desc, cutoff, binsize),".csv", sep=""), na.strings=c(".", "NA", "", "?"))
  # if nameIndex is not given then use the last col as PDB_ID
  if (is.na(nameIndex)) { 
    nameIndex   = length(testData[1,])
  }
  newTrainData = createRowsComplementsData(trainData, testData, nameIndex)
  # write to file
  write.table(testData		  , file = paste(path,prefix,"test_"    ,concatCutoffAndBinsize(desc, cutoff, binsize),".csv", sep=""), sep = ",", row.names = FALSE)
  write.table(newTrainData	, file = paste(path,prefix,"training_",concatCutoffAndBinsize(desc, cutoff, binsize),".csv", sep=""), sep = ",", row.names = FALSE)
}

# create training and test data sets from checkCommonRows of 2 given sets
# USAGE: 
#source("/Users/knight/dev/classification_R/R/lib.data.r")
#createTrainingTestDataCSAR(path = "/Users/knight/dev/classification_R/data/CSAR_raw/", prefix = "CSAR", trainPath = "/Users/knight/dev/classification_R/data/CASF12_raw/CASF12_refined", desc = "elements", binsize = 12)
#createTrainingTestDataCSAR(path = "/Users/knight/dev/classification_R/data/CSAR_raw/", prefix = "CSAR", trainPath = "/Users/knight/dev/classification_R/data/CASF12_raw/CASF12_refined", desc = "sybyl", binsize = 12)
#createTrainingTestDataCSAR(path = "/Users/knight/dev/classification_R/data/CSAR_raw/", prefix = "CSAR", trainPath = "/Users/knight/dev/classification_R/data/CASF12_raw/CASF12_refined", desc = "credo", binsize = 12)
createTrainingTestDataCSAR <- function(path, prefix, trainPath, desc = "elements", cutoff = 12, binsize = 12) {
  core1   = read.csv(paste(trainPath,prefix,"set1_",concatCutoffAndBinsize(desc, cutoff, binsize),".csv", sep=""), na.strings=c(".", "NA", "", "?"))
  core2   = read.csv(paste(trainPath,prefix,"set2_",concatCutoffAndBinsize(desc, cutoff, binsize),".csv", sep=""), na.strings=c(".", "NA", "", "?")) 
  refine  = read.csv(paste(trainPath,prefix="",dataName="_",concatCutoffAndBinsize(desc, cutoff, binsize),".csv", sep=""), na.strings=c(".", "NA", "", "?"))
  core    = rbind(core1, core2)
  nameIndex = length(core[1,])
  # \TODO: quick and dirty to remove FXa from CSAR set
  FXa     = read.table("/Users/knight/dev/classification_R/data/FXa.txt")
  FXa 		= FXa[,1]
  core		= core[!core[,nameIndex] %in% FXa, ]

  testData  = core
  trainData = refine
  nameIndex = length(testData[1,])
  source("/Users/knight/dev/classification_R/R/lib.data.r")
  newTrainData = createRowsComplementsData(trainData, testData, nameIndex)
  
  # write to file
  write.table(testData		  ,file = paste(path,prefix,"13_test_"    ,concatCutoffAndBinsize(desc, cutoff, binsize),".csv", sep=""), sep = ",", row.names = FALSE)
  write.table(newTrainData	,file = paste(path,prefix,"13_training_",concatCutoffAndBinsize(desc, cutoff, binsize),".csv", sep=""), sep = ",", row.names = FALSE)
}

# combinate 2 training data sets to 1 set, 
# \RETURN the combi set
# USAGE:
#createMixTrainingTestData(path = "/home/dat/WORK/DB/DESCRIPTORS/", testPrefix = "CASF12", desc = "credo", binsize = 0)
#createMixTrainingTestData(path = "/home/dat/WORK/DB/DESCRIPTORS/", testPrefix = "CASF12", desc = "elements", binsize = 0)
#createMixTrainingTestData(path = "/home/dat/WORK/DB/DESCRIPTORS/", testPrefix = "CASF12", desc = "sybyl", binsize = 0)

createMixTrainingData <- function(path, prefix = c("CASF12/CASF12_training_","CASF13/CASF13_training_"), desc = "elements", cutoff = 12, binsize = 12) {	
	data1 = read.csv(concatPath(path, prefix[1], cutoff, binsize, desc), na.strings=c(".", "NA", "", "?"))
	data2 = read.csv(concatPath(path, prefix[2], cutoff, binsize, desc), na.strings=c(".", "NA", "", "?"))
	
	nameIndex = length(data1[1,])
	
	complementsData1 = createRowsComplementsData(data1, data2, nameIndex)
	
	mixData = rbind(complementsData1, data2)
	
	return (mixData)
}

# create training and test data sets from mixing training sets and checkCommonRows of 2 given sets
# USAGE: 
#source("/home/dat/WORK/dev/R/lib.data.r")
#createMixTrainingTestData(path = "/home/dat/WORK/dev/data/", testPrefix = "CASF12", desc = "elements", binsize = 12)
#createMixTrainingTestData(path = "/home/dat/WORK/dev/data/", testPrefix = "CASF12", desc = "credo", binsize = 12)
#createMixTrainingTestData(path = "/home/dat/WORK/dev/data/", testPrefix = "CASF12", desc = "sybyl", binsize = 12)
createMixTrainingTestData <- function(path, prefix = c("CASF12/CASF12_training_","CASF13/CASF13_training_"), testPrefix = "CASF12", desc = "elements", cutoff = 12, binsize = 12) {
	totalTestPrefix = paste(testPrefix,"/",testPrefix,"_test_",sep="")
	testData = read.csv(concatPath(path, totalTestPrefix, cutoff, binsize, desc), na.strings=c(".", "NA", "", "?"))
	trainData = createMixTrainingData(path, prefix, desc, cutoff, binsize)
	
	nameIndex = length(testData[1,])
	newTrainData = createRowsComplementsData(trainData, testData, nameIndex)
	
	newTestPrefix = paste(testPrefix,"/mix",testPrefix,"_",sep="")
	write.table(newTrainData	, file = paste(path,prefix=newTestPrefix,dataName="training_",concatCutoffAndBinsize(desc, cutoff, binsize),".csv", sep=""), sep = ",", row.names = FALSE)
  #	write.table(testData		, file = concatPathWithDataName(path, prefix = newTestPrefix, "test", desc, cutoff, binsize), sep = ",", row.names = FALSE)
}


# create combinations data from merging 2 data sets
#
# USAGE: 
#createCombiData(path = "/home/dat/WORK/dev/rfscore/bin/", prefix = "CASF12_core_", descName = c("elementsv2", "SIFt"))
#createCombiData(path = "/home/dat/WORK/dev/rfscore/bin/", prefix = "CASF12_refined_", descName = c("elementsv2", "SIFt"))
#
createCombiData <- function(path, prefix, cutoff = 12, binsize = 0, descName = c("elementsv2", "SIFt"), postfix="") {
  data1 = read.csv(concatPath(path, prefix, cutoff, binsize, descName[1], postfix), na.strings=c(".", "NA", "", "?"))
  data2 = read.csv(concatPath(path, prefix, cutoff, binsize, descName[2], postfix), na.strings=c(".", "NA", "", "?"))
  # \IMPORTANT remote the col which contains pKd information when combine the standard CASF data
  if (postfix == "") { # remove the pKd col
    data2 = data2[,-1]
    postfix = paste(postfix, "_process", sep="")
  }
  nameIndex1 = length(data1[1,])
  nameIndex2 = length(data2[1,])
  mergeData12   = merge(data1, data2, by.x = nameIndex1, by.y = nameIndex2)
  colnames(mergeData12)[2] = "pKd.pKi"  
  # write to file
  write.table(mergeData12, file = paste(path,prefix,descName[1],"-",descName[2],"_c",cutoff,"b",binsize, postfix,".csv", sep=""), sep = ",", row.names = FALSE)  
}

# process the test set by adding docking scores/pIC50 to 2. col (like standard CASF data structure and make this scores unique
#           
# USAGE: 
#createCombiData(path = "/home/dat/WORK/dev/rfscore/bin/", prefix = "CASF12_core_", descName = c("elementsv2", "SIFt"))
#
processTestData <- function(path, prefix, cutoff = 12, binsize = 0, descName, postfix, scoreFilename) {  
  data = read.csv(concatPath(path, prefix, cutoff, binsize, descName, postfix), na.strings=c(".", "NA", "", "?"))
  #scores = read.csv(paste(path, "JMJ", postfix, ".csv", sep=""))
  #data = read.csv(concatPath(path, prefix, cutoff, binsize, descName, paste(postfix,"",sep="")), na.strings=c(".", "NA", "", "?"))
  scores = read.csv(paste(path, scoreFilename, postfix, ".csv", sep=""))
  
  if (length(scores[1,]) > 2) { # if scores has more than 2 cols, then remove the last one (not-used)
    scores = scores[, -length(scores[1,])]
  }
  
  colnames(scores)[2] = "pKd.pKi"  
  scores = makeUniqueList(scores, col = 2)
  
  mergeData   = merge(scores, data, by.x = 1, by.y = 1)
  colnames(mergeData)[1] = 'pdb'
  # write to file
  write.table(mergeData, file = gsub(".csv", "_process.csv", concatPath(path, prefix, cutoff, binsize, descName, postfix)), sep = ",", row.names = FALSE)  
  
}
  

# create combinations data from merging 3 data sets
#
# USAGE: 
#source("/Users/knight/dev/classification_R/lib.data.r")
#createMergeData(path = "/Users/knight/dev/classification_R/data/CASF13/", prefix = "CASF13_test_", binsize = 1)
#createMergeData(path = "/Users/knight/dev/classification_R/data/CASF13/", prefix = "CASF13_training_", binsize = 1)
#createMergeData(path = "/Users/knight/dev/classification_R/data/CSAR12/", prefix = "CSAR12_test_", binsize = 12)
#createMergeData(path = "/Users/knight/dev/classification_R/data/CSAR12/", prefix = "CSAR12_training_", binsize = 12)
#
createMergeData <- function(path, prefix, cutoff = 12, binsize = 12, descName = c("elements", "sybyl", "credo")) {
  data1 = read.csv(concatPath(path, prefix, cutoff, binsize, descName[1]), na.strings=c(".", "NA", "", "?"))
  data2 = read.csv(concatPath(path, prefix, cutoff, binsize, descName[2]), na.strings=c(".", "NA", "", "?"))
  data3 = read.csv(concatPath(path, prefix, cutoff, binsize, descName[3]), na.strings=c(".", "NA", "", "?"))
  # \IMPORTANT remote the col which contains pKd information 
  data2_removePkD = data2[,-1]
  data3_removePkD = data3[,-1]
  nameIndex1 = length(data1[1,])
  nameIndex2 = length(data2_removePkD[1,])
  nameIndex3 = length(data3_removePkD[1,])
  mergeData12   = merge(data1, data2_removePkD, by.x = nameIndex1, by.y = nameIndex2)
  mergeData13   = merge(data1, data3_removePkD, by.x = nameIndex1, by.y = nameIndex3)
  mergeData23   = merge(data2, data3_removePkD, by.x = (nameIndex2+1), by.y = nameIndex3)
  mergeData123  = merge(mergeData12, data3_removePkD, by.x = 1, by.y = nameIndex3)
  # write to file
  write.table(mergeData12, file = paste(path,prefix,descName[1],"-",descName[2],"_c",cutoff,"b",binsize,".csv", sep=""), sep = ",", row.names = FALSE)
  write.table(mergeData13, file = paste(path,prefix,descName[1],"-",descName[3],"_c",cutoff,"b",binsize,".csv", sep=""), sep = ",", row.names = FALSE)
  write.table(mergeData23, file = paste(path,prefix,descName[2],"-",descName[3],"_c",cutoff,"b",binsize,".csv", sep=""), sep = ",", row.names = FALSE)
  write.table(mergeData123, file = paste(path,prefix,descName[1],"-",descName[2],"-",descName[3],"_c",cutoff,"b",binsize,".csv", sep=""), sep = ",", row.names = FALSE)
}
