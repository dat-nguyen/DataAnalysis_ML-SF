#source("lib/lib.validation_extern.R")
source("lib/lib.process_data.R")
source("lib/configuration.R")
path = "/home/dat/WORK/DB/DESCRIPTORS/CASF/"

processTestDataWithoutScore <- function(path, target, descName, postfix) {  
  data    = read.csv(paste(path, target, descName, postfix, ".csv", sep=""), na.strings=c(".", "NA", "", "?"))
  ligandNames = data[,1]
  print(length(ligandNames))
  scores  = data.frame("pdb" = ligandNames, "pKd.pKi" = 1:length(ligandNames))
  
  mergeData   = merge(scores, data, by.x = 1, by.y = 1)
  colnames(mergeData)[1] = 'pdb'
  # write to file
  write.table(mergeData, file = paste(path, target, descName, postfix, "_process.csv", sep=""), sep = ",", row.names = FALSE)  
  
}
# create combination datas from combining 2 types of descriptors (elementsv2 + SIFt/SIFtv2)
runCombiData <- function() {
  
  # only for CASF data
  #for (pre in c("CASF12_core_", "CASF13_core_", "CASF14_core_", "CASF12_refined_", "CASF13_refined_", "CASF14_refined_") ) {
  for (CASFset in c("CASF12_refined_", "CASF13_refined_", "CASF14_refined_") ) {
    #print(CASFset)
    #createCombiData(path = path2Name, prefix = CASFset, descName = c("elementsv2", "SIFt"))
    #write.table(merge3DescCASF(path, CASFset), file = paste(path, CASFset, concatCutoffAndBinsize("elementsv2-SIFt"), "-xscore_process.csv", sep=""), sep = ",", row.names = FALSE)  

  }

  # for all test data
  #for (target in c("3PDQ_new_", "T36_JMJ_Xray_") ) {
  #target = "DIG10.2_"
  #path2Name = "/home/dat/WORK/DB/DESCRIPTORS/DIG10.2/confgen/"
          
    for (testset in dockingMethods) {      
    #for (post in c("_actives", "_inactives", "_test_gold")) {
      #print(testset)
      testset = paste('_', testset, sep='')
#      createCombiData(path = path2Name, prefix = target, descName = c("elementsv2", "SIFt"), postfix = testset)      
#      processTestData(path = path2Name, prefix = target, descName = "elementsv2-SIFt", postfix = testset, scoreFilename="cs-confgen_all")      
#      write.table(merge3DescCASF(path2Name, target, testset), file = paste(path2Name, target, concatCutoffAndBinsize("elementsv2-SIFt-xscore"), testset,".csv", sep=""), sep = ",", row.names = FALSE)        
#      processTestData(path = path2Name, prefix = target, descName = "elementsv2-SIFt-xscore", postfix = testset, scoreFilename="cs-confgen_all")
      
#      createCombiData(path = path2Name, prefix = target, descName = c("elementsv2", "SIFt"), postfix = post)
#      processTestData(path = path2Name, prefix = target, descName = "  elementsv2-SIFt", postfix = post)
    }
  target = "DUD-E_"
  path2Name = "/home/dat/WORK/DB/DESCRIPTORS/DUD-E/"
  for (post in c("_RENI", "_FGFR1")) { #, "_ADA"
    createCombiData(path = path2Name, prefix = target, descName = c("elementsv2", "SIFt"), postfix = post) 
    processTestDataWithoutScore(path2Name, target, "elementsv2-SIFt_c12b0", postfix = post)
    write.table(merge3DescCASF(path2Name, target, post), file = paste(path2Name, target, concatCutoffAndBinsize("elementsv2-SIFt-xscore"), post, ".csv", sep=""), sep = ",", row.names = FALSE)   
    processTestDataWithoutScore(path2Name, target, "elementsv2-SIFt-xscore_c12b0", postfix = post)
  }
    
}

# create only-metallo data for training purpose
createMetallo <- function(path, prefix, desc, cutoff = 12, binsize = 0) {
  dataPath = paste(path,prefix,"refined_",concatCutoffAndBinsize(desc, cutoff, binsize),".csv", sep="")
  listPath = paste(path,"ID_",prefix,"metallo_",concatCutoffAndBinsize(desc, cutoff, binsize),".csv", sep="")
  output   = paste(path,prefix,"metallo_",concatCutoffAndBinsize(desc, cutoff, binsize),".csv", sep="")
  
  filterDataRow(dataPath, listPath, output, nameIndex=1)  
}

# create training and test set by checking the common rows of two sets and remove it from the training set
for (pre in dataSets ) {
#  print(pre)
  #path = "/home/dat/WORK/DB/DESCRIPTORS/New/"
  #  createTrainingTestData(path, prefix= pre, desc = "elementsv2-SIFt", nameIndex = 1)
  #  createTrainingTestData(path, prefix= pre, desc = "elementsv2-SIFtv2", nameIndex = 1)
}

for (prefix in dataSets ) {
#  print(prefix)
  for (desc in c("elementsv2-SIFt", "elementsv2-SIFtv2") ) {
    #createMetallo(path, prefix, desc)
  }
}

# 
mergeDesc <- function(path1, path2, nameIndex1, nameIndex2, removeIndex = 0) {
  data1 = read.csv(path1, na.strings=c(".", "NA", "", "?"))
  data2 = read.csv(path2, na.strings=c(".", "NA", "", "?"))
  if (missing(nameIndex1)) # if the index is missing, then it's the last col of the set
    nameIndex1 = length(data1[1,])
  if (missing(nameIndex2)) # if the index is missing, then it's the last col of the set
    nameIndex2 = length(data2[1,])
  if (removeIndex > 0) { # sometimes we need to remove a col @ removeIndex @ the second dataset
    data2 = data2[,-removeIndex]
    if (nameIndex2 > removeIndex) # if the col to be merged from second set is right from the removed col, need to adjust its position as well
      nameIndex2 = nameIndex2-1
  }
  mergeData = merge(data1, data2, by.x = nameIndex1, by.y = nameIndex2)
  
  return (mergeData)
}

# 
merge3DescCASF <- function(path, CASFset, testset = "") {
  desc1 = paste(path, CASFset, concatCutoffAndBinsize("elementsv2-SIFt"), testset, '.csv', sep='')
  desc2 = paste(path, CASFset, 'xscore', testset, '.csv', sep='')
  mergeData = mergeDesc(desc1, desc2, nameIndex1=1,nameIndex2=1)

  return (mergeData)
}
createRMSDdata <- function() {
  # only for CASF data
  #path = "/home/dat/WORK/DB/DESCRIPTORS/CASF/RMSD/"
  path = "/Users/knight/Dropbox/data/RMSD/"
  #for (pre in c("CASF12_core_", "CASF13_core_", "CASF14_core_", "CASF12_refined_", "CASF13_refined_", "CASF14_refined_") ) {
  for (CASFset in c("CASF14_refined_RMSD_") ) {
    #print(CASFset)
    #createCombiData(path, prefix = CASFset, descName = c("elementsv2", "SIFt"))
    #write.table(merge3DescCASF(path, CASFset), file = paste(path, CASFset, concatCutoffAndBinsize("elementsv2-SIFt"), "-xscore_process.csv", sep=""), sep = ",", row.names = FALSE)

  }
}
createRMSDdata()
#runCombiData()
