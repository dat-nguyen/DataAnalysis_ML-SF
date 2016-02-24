# constants file
###############################################################
DESC_PATH           = "/home/dat/WORK/DB/DESCRIPTORS/"
RMSD_TRAINING_PATH  = paste0(DESC_PATH, "RMSD/")
RMSD_TEST_PATH      = DESC_PATH

OUTPUT_RMSD         = "/home/dat/WORK/output/RMSD/"

PROCESSED_DATA_PATH = paste0(DESC_PATH, "Processed/")

DESCRIPTORS         = c("elementsv2", "SIFt")
DESC_SETS           = c("elementsv2-SIFt")
TRAINING_SETS       = c("sampling_clusters10", "sampling_100")
TRAINING_SETS_NAME  = c("setA", "setB")
ML_METHODS          = c("RoF_REPT", "RoF_RoT", "RF_")

CASF_NAMES          = c("2007", "2012", "2013", "2014")
CASF_YEARS          = c("v2007", "v2012-refined", "v2013-refined", "v2014-refined")
#CASF_YEARS          = c("v2012-refined")
CASF_SETS           = paste0("CASF", CASF_YEARS)
###############################################################
TARGET_LIST_FIDELE  = c("5P21", "4BBG", "3KKP", "1GS4", "2X9E", "3E37")
TARGET_DB_FIDELE    = c("africa", "npact")
POSES_GEN_LIST_FIDELE  = c("GS", "XP")

FideleList = sapply(TARGET_LIST_FIDELE, paste, sapply(TARGET_DB_FIDELE, paste, POSES_GEN_LIST_FIDELE, sep="-"), sep="-") 

TEST_SETS_KDM = c("T36_JMJ_Xray_test_gold", "T36_JMJ_Xray_actives", "T36_JMJ_Xray_inactives")
TEST_SETS_DIG = c("DIG_XP", "DIG_SP", "DIG_goldscore")

TARGET_LIST = c("Aminothiazoles_71_PB_md1_rst_Dat", "SmHDAC8_inhibitors", 
                TEST_SETS_KDM, 
                TEST_SETS_DIG)
###############################################################
TRUE_RANK_DIG10.2 = c(7, 8, 6, 4, 3, 5, 2, 1, 10, 9)

RESULT_PATH             = "/home/dat/WORK/RESULTS/2016-02-23/"
SPLIT_RESULT_PATH       = paste0(RESULT_PATH, "splitted_results/")
PREDICTED_SCORES_PATH   = paste0(RESULT_PATH, "predicted_scores/")
TRUE_SCORES_PATH        = paste0(RESULT_PATH, "true_scores/")

###############################################################
