###############################################################
# preprocessing CASF and test data for MLSF Project, Feb-2016
###############################################################
source("libs/libProcessData.R")
###############################################################
createTrainingData <- function(CASFset) {
  for (trainingSet in TRAINING_SETS) {
    mergeDesc(RMSD_TRAINING_PATH, CASFset, trainingSet)
  }
}
###############################################################
createTestData <- function(dirPath, testNamePrefix) {
  path1 = paste(dirPath, testNamePrefix, "_", DESCRIPTORS[1], ".csv", sep="")
  path2 = paste(dirPath, testNamePrefix, "_", DESCRIPTORS[2], ".csv", sep="")
  data1 = read.csv(path1, na.strings=c(".", "NA", "", "?"))
  data2 = read.csv(path2, na.strings=c(".", "NA", "", "?"))      
  mergePath12 = paste(PROCESSED_DATA_PATH, testNamePrefix, "_", DESCRIPTORS[1], "-", DESCRIPTORS[2], ".csv", sep="")
  mergeData12 = mergeRowData(path1, path2, createFakeValue = TRUE)
  write.table(mergeData12, file=mergePath12, sep=",", quote=FALSE, row.names=FALSE)
}
###############################################################
createTestDataAll <- function() {
  ############### anti-cancer targets from Fidele ###############
  #target = TARGET_LIST_FIDELE[1]
  #target_db = TARGET_DB_FIDELE[1]
  #pose  = POSES_GEN_LIST[2]
  for (target in TARGET_LIST_FIDELE) {
    for (target_db in TARGET_DB_FIDELE)  {
      for (pose in POSES_GEN_LIST_FIDELE) {
        testNamePrefix  = paste(target, "-", target_db, "-", pose, sep="")
        RMSD_testPath   = paste(RMSD_TEST_PATH, "Fidele/", sep="")
        createTestData(RMSD_testPath, testNamePrefix)
      }
    }
  }
  ############### Sirt2 from Berin ###############
  testNamePrefix  = "Aminothiazoles_71_PB_md1_rst_Dat"
  RMSD_testPath   = paste(RMSD_TEST_PATH, "Sirt2/", sep="")
  createTestData(RMSD_testPath, testNamePrefix)
  ############### SmHDAC8 from Jelena ###############
  testNamePrefix  = "SmHDAC8_inhibitors"
  RMSD_testPath   = paste(RMSD_TEST_PATH, "SmHDAC8/", sep="")
  createTestData(RMSD_testPath, testNamePrefix)
  ############### KDMs from Prof. Sippl ###############
  RMSD_testPath   = paste(RMSD_TEST_PATH, "KDMs/", sep="")
  for (testSet in TEST_SETS_KDM) {
    testNamePrefix  = paste("T36_JMJ_Xray_",testSet, sep="")
    createTestData(RMSD_testPath, testNamePrefix)
  }
  ############### DIG from CSAR2013 ###############
  RMSD_testPath   = paste(RMSD_TEST_PATH, "DIG10.2/", sep="")
  for (testSet in c("goldscore", "XP", "SP")) {
    testNamePrefix  = paste("DIG_",testSet, sep="")
    createTestData(RMSD_testPath, testNamePrefix)
  }
  
}      
###############################################################
combiAllTestData <- function(testPath, combiName) {
  pattern = paste(DESCRIPTORS[1], "-", DESCRIPTORS[2], ".csv", sep="")
  allFiles = list.files(testPath, pattern)
  # create the base data for the first file
  file = allFiles[1]
  path = paste(testPath, file, sep="")
  combiData = read.csv(path, na.strings=c(".", "NA", "", "?"))
  # find the basename
  namePrefix = gsub(paste("_", pattern, sep=""), "", file)
  # add the basename to the ID for identifying
  combiData[,1] = paste(namePrefix, "/", combiData[,1], sep="")
  # merge all other data to the base data (by column)
  for (file in allFiles[-1]) {
    newPath = paste(testPath, file, sep="")
    newData = read.csv(newPath, na.strings=c(".", "NA", "", "?"))
    # find the basename
    namePrefix = gsub(paste("_", pattern, sep=""), "", file)
    # add the basename to the ID for identifying
    newData[,1] = paste(namePrefix, "/", newData[,1], sep="")
    print(newPath)
    combiData = rbind(combiData, newData)
  }
  combiPath = paste(testPath, paste(combiName, "_", pattern, sep=""), sep="")
  write.table(combiData, file = combiPath, sep = ",", quote=FALSE, row.names = FALSE)
  
}
###############################################################

#createTrainingData(CASF_SETS[4])
#createTestDataAll()
#combiAllTestData(PROCESSED_DATA_PATH, "targets_2016-02-23")
#combiAllTestData("/home/dat/WORK/DB/DESCRIPTORS/Processed/Fidele/", "targets_Fidele_2016-02-23")

