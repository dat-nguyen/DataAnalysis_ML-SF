###############################################################
# preprocessing CASF data for MLSF Project, Feb-2016
###############################################################
RMSD_TRAINING_PATH  = "/home/dat/WORK/DB/DESCRIPTORS/RMSD/"

DESCRIPTORS         = c("elementsv2", "SIFt")
TRAINING_SETS       = c("sampling_clusters10", "sampling_100")
CASF_SETS           = c("CASFv2007", "CASFv2012-refined", "CASFv2013-refined", "CASFv2014-refined")
TARGET_LIST_FIDELE  = c("5P21", "4BBG", "3KKP", "1GS4", "2X9E", "3E37")
TARGET_DB_FIDELE    = c("africa", "npact")
POSES_GEN_LIST      = c("GS", "XP")

source("libs/libProcessData.R")

RMSD_testPath       = "/home/dat/WORK/DB/DESCRIPTORS/Fidele/"

createTrainingData <- function(CASFset) {
  for (trainingSet in TRAINING_SETS) {
    mergeDesc(RMSD_TRAINING_PATH, CASFset, trainingSet)
  }
}

#createTrainingData(CASF_SETS[2])

createTestData <- function() {
  #target = TARGET_LIST_FIDELE[1]
  #target_db = TARGET_DB_FIDELE[1]
  #pose  = POSES_GEN_LIST[2]
  for (target in TARGET_LIST_FIDELE) {
    for (target_db in TARGET_DB_FIDELE)  {
      for (pose in POSES_GEN_LIST) {
        path1 = paste(RMSD_testPath, target, "-", target_db, "-", pose, "_", DESCRIPTORS[1], ".csv", sep="")
        path2 = paste(RMSD_testPath, target, "-", target_db, "-", pose, "_", DESCRIPTORS[2], ".csv", sep="")
        data1 = read.csv(path1, na.strings=c(".", "NA", "", "?"))
        data2 = read.csv(path2, na.strings=c(".", "NA", "", "?"))      
        mergePath12 = paste(RMSD_testPath, "processed/", target, "-", target_db, "-", pose, "_", DESCRIPTORS[1], "-", DESCRIPTORS[2], ".csv", sep="")
        mergeData12 = mergeRowData(path1, path2, createFakeValue = TRUE)
        write.table(mergeData12, file=mergePath12, sep=",", quote=FALSE, row.names=FALSE)
      }
    }
  }
}      
      #path1 = paste(RMSD_testPath, target, "-", target_db, "-", POSES_GEN_LIST[1], "_", DESCRIPTORS[2], ".csv", sep="")
      #path2 = paste(RMSD_testPath, target, "-", target_db, "-", POSES_GEN_LIST[2], "_", DESCRIPTORS[2], ".csv", sep="")
 #     data1 = read.csv(path1, na.strings=c(".", "NA", "", "?"))
#      data2 = read.csv(path2, na.strings=c(".", "NA", "", "?"))      
#      mergePath12 = paste(RMSD_testPath, target, "-", target_db, "-", pose, "_", DESCRIPTORS[1], "-", DESCRIPTORS[2], ".csv", sep="")
      #mergeColData12 = mergeColData(path1, path2)

#createTestData()

