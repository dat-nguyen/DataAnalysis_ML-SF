# constants file
###############################################################
DESC_PATH           = "/home/dat/WORK/DB/DESCRIPTORS/"
RMSD_TRAINING_PATH  = paste0(DESC_PATH, "RMSD/")
RMSD_TEST_PATH      = DESC_PATH

PROCESSED_DATA_PATH = paste0(DESC_PATH, "/Processed/")

DESCRIPTORS         = c("elementsv2", "SIFt")
DESC_SETS           = c("elementsv2-SIFt")
TRAINING_SETS       = c("sampling_clusters10", "sampling_100")
ML_METHODS          = c("RoF_REPT", "RoF_RoT", "RF_")
CASF_SETS           = c("CASFv2007", "CASFv2013-refined", "CASFv2014-refined") #"CASFv2012-refined", 
###############################################################
TARGET_LIST_FIDELE  = c("5P21", "4BBG", "3KKP", "1GS4", "2X9E", "3E37")
TARGET_DB_FIDELE    = c("africa", "npact")
POSES_GEN_LIST_FIDELE  = c("GS", "XP")

TEST_SETS_KDM = c("T36_JMJ_Xray_test_gold", "T36_JMJ_Xray_actives", "T36_JMJ_Xray_inactives")
TEST_SETS_DIG = c("DIG_XP", "DIG_SP", "DIG_goldscore")
###############################################################
TRUE_RANK_DIG10.2 = c(7, 8, 6, 4, 3, 5, 2, 1, 10, 9)

RESULT_PATH       = "/home/dat/WORK/RESULTS/2016-02-23/"
SPLIT_RESULT_PATH = paste0(RESULT_PATH, "splitted_results/")
TRUE_SCORES_PATH  = paste0(RESULT_PATH, "true_scores/")

###############################################################
