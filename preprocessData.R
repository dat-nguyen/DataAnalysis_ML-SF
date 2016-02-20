###############################################################
# preprocessing CASF data for MLSF Project, Feb-2016
###############################################################
RMSD_TRAINING_PATH  = "/home/dat/WORK/DB/DESCRIPTORS/RMSD/"
RMSD_TEST_PATH      = "/home/dat/WORK/DB/DESCRIPTORS/"

PROCESSED_DATA_PATH = "/home/dat/WORK/DB/DESCRIPTORS/Processed/"

DESCRIPTORS         = c("elementsv2", "SIFt")
TRAINING_SETS       = c("sampling_clusters10", "sampling_100")
CASF_SETS           = c("CASFv2007", "CASFv2012-refined", "CASFv2013-refined", "CASFv2014-refined")
TARGET_LIST_FIDELE  = c("5P21", "4BBG", "3KKP", "1GS4", "2X9E", "3E37")
TARGET_DB_FIDELE    = c("africa", "npact")
POSES_GEN_LIST_FIDELE  = c("GS", "XP")

TEST_SETS_KDM = c("test_gold", "actives", "inactives")

source("libs/libProcessData.R")


createTrainingData <- function(CASFset) {
  for (trainingSet in TRAINING_SETS) {
    mergeDesc(RMSD_TRAINING_PATH, CASFset, trainingSet)
  }
}

#createTrainingData(CASF_SETS[3])

createTestData <- function(dirPath, testNamePrefix) {
  path1 = paste(dirPath, testNamePrefix, "_", DESCRIPTORS[1], ".csv", sep="")
  path2 = paste(dirPath, testNamePrefix, "_", DESCRIPTORS[2], ".csv", sep="")
  data1 = read.csv(path1, na.strings=c(".", "NA", "", "?"))
  data2 = read.csv(path2, na.strings=c(".", "NA", "", "?"))      
  mergePath12 = paste(PROCESSED_DATA_PATH, testNamePrefix, "_", DESCRIPTORS[1], "-", DESCRIPTORS[2], ".csv", sep="")
  mergeData12 = mergeRowData(path1, path2, createFakeValue = TRUE)
  write.table(mergeData12, file=mergePath12, sep=",", quote=FALSE, row.names=FALSE)
}

createTestDataAll <- function() {
  ############### anti-cancer targets from Fidele ###############
  #target = TARGET_LIST_FIDELE[1]
  #target_db = TARGET_DB_FIDELE[1]
  #pose  = POSES_GEN_LIST[2]
  for (target in TARGET_LIST_FIDELE) {
    for (target_db in TARGET_DB_FIDELE)  {
      for (pose in POSES_GEN_LIST_FIDELE) {
        testNamePrefix  = paste(target, "-", target_db, "-", pose, sep="")
        RMSD_testPath   = paste(RMSD_TEST_PATH, "Fidele", sep="")
        #createTestData(RMSD_testPath, testNamePrefix)
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
    print(testNamePrefix)
    createTestData(RMSD_testPath, testNamePrefix)
  }
  
}      

      #path1 = paste(RMSD_testPath, target, "-", target_db, "-", POSES_GEN_LIST[1], "_", DESCRIPTORS[2], ".csv", sep="")
      #path2 = paste(RMSD_testPath, target, "-", target_db, "-", POSES_GEN_LIST[2], "_", DESCRIPTORS[2], ".csv", sep="")
 #     data1 = read.csv(path1, na.strings=c(".", "NA", "", "?"))
#      data2 = read.csv(path2, na.strings=c(".", "NA", "", "?"))      
#      mergePath12 = paste(RMSD_testPath, target, "-", target_db, "-", pose, "_", DESCRIPTORS[1], "-", DESCRIPTORS[2], ".csv", sep="")
      #mergeColData12 = mergeColData(path1, path2)

createTestDataAll()

