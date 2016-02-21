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
  path1 = paste0(dirPath, testNamePrefix, "_", DESCRIPTORS[1], ".csv")
  path2 = paste0(dirPath, testNamePrefix, "_", DESCRIPTORS[2], ".csv")
  data1 = read.csv(path1, na.strings=c(".", "NA", "", "?"))
  data2 = read.csv(path2, na.strings=c(".", "NA", "", "?"))      
  mergePath12 = paste0(PROCESSED_DATA_PATH, testNamePrefix, "_", DESCRIPTORS[1], "-", DESCRIPTORS[2], ".csv")
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
        testNamePrefix  = paste0(target, "-", target_db, "-", pose)
        RMSD_testPath   = paste0(RMSD_TEST_PATH, "Fidele/")
        createTestData(RMSD_testPath, testNamePrefix)
      }
    }
  }
  ############### Sirt2 from Berin ###############
  RMSD_testPath   = paste0(RMSD_TEST_PATH, "Sirt2/")
  testNamePrefix  = "Aminothiazoles_71_PB_md1_rst_Dat"
  createTestData(RMSD_testPath, testNamePrefix)
  ############### SmHDAC8 from Jelena ###############
  RMSD_testPath   = paste0(RMSD_TEST_PATH, "SmHDAC8/")
  testNamePrefix  = "SmHDAC8_inhibitors"
  createTestData(RMSD_testPath, testNamePrefix)
  ############### KDMs from Prof. Sippl ###############
  RMSD_testPath   = paste0(RMSD_TEST_PATH, "KDMs/")
  for (testNamePrefix in TEST_SETS_KDM) {
    createTestData(RMSD_testPath, testNamePrefix)
  }
  ############### DIG from CSAR2013 ###############
  RMSD_testPath   = paste0(RMSD_TEST_PATH, "DIG10.2/")
  for (testNamePrefix in TEST_SETS_DIG) {
    createTestData(RMSD_testPath, testNamePrefix)
  }
  
}      
###############################################################
combiAllTestData <- function(testPath, combiName) {
  pattern = paste0(DESCRIPTORS[1], "-", DESCRIPTORS[2], ".csv")
  allFiles = list.files(testPath, pattern)
  # create the base data for the first file
  file = allFiles[1]
  path = paste0(testPath, file)
  combiData = read.csv(path, na.strings=c(".", "NA", "", "?"))
  # find the basename
  namePrefix = gsub(paste0("_", pattern), "", file)
  # add the basename to the ID for identifying
  combiData[,1] = paste0(namePrefix, "/", combiData[,1])
  # merge all other data to the base data (by column)
  for (file in allFiles[-1]) {
    newPath = paste0(testPath, file)
    newData = read.csv(newPath, na.strings=c(".", "NA", "", "?"))
    # find the basename
    namePrefix = gsub(paste0("_", pattern), "", file)
    # add the basename to the ID for identifying
    newData[,1] = paste0(namePrefix, "/", newData[,1])
    print(newPath)
    combiData = rbind(combiData, newData)
  }
  combiPath = paste0(testPath, paste0(combiName, "_", pattern))
  write.table(combiData, file = combiPath, sep = ",", quote=FALSE, row.names = FALSE)
  
}
###############################################################

#createTrainingData(CASF_SETS[4])
#createTestDataAll()
#combiAllTestData(PROCESSED_DATA_PATH, "targets_2016-02-23")
#combiAllTestData("/home/dat/WORK/DB/DESCRIPTORS/Processed/Fidele/", "targets_Fidele_2016-02-23")

