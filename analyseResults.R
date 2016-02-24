###############################################################
# analyse the results for MLSF Project, Feb-2016
###############################################################
source("libs/libProcessData.R")
###############################################################
splitAllResults <- function() {
  for (CASFset in CASF_SETS) {
    for (trainSet in TRAINING_SETS) {
      for (method in ML_METHODS) {
        classifyMethod = paste(CASFset, "_", trainSet, "_", method, sep="")
        splitWekaResults(IDPath = PROCESSED_DATA_PATH, classifyMethod, combiName="targets_2016-02-23", TARGET_LIST)
      }
    }
  } 
}
###############################################################
splitAllResults_Fidele <- function() {
#  FideleList = c()
#  for (target in TARGET_LIST_FIDELE) 
#    for (target_db in TARGET_DB_FIDELE)  
#      for (pose in POSES_GEN_LIST_FIDELE) 
#          FideleList = c(FideleList, paste(target, target_db, pose, sep="-"))
  for (CASFset in CASF_SETS) 
    for (trainSet in TRAINING_SETS)
      for (method in ML_METHODS) {
        classifyMethod = paste(CASFset, "_", trainSet, "_", method, sep="")
        splitWekaResults(IDPath = paste0(PROCESSED_DATA_PATH, "Fidele/"), classifyMethod, combiName="targets_Fidele_2016-02-23", FideleList)
      }
   
}
###############################################################
createRMSDhistogram <- function(RMSDpath, pattern="sampling_cluster", mainText) {
  allFiles = list.files(RMSDpath, pattern=pattern)
  # create empty dataframe
  RMSDall = read.csv(text="x,RMSD")
  for (CSVfile in allFiles) {
    RMSD = read.csv(paste0(RMSDpath, CSVfile), header=FALSE)
    RMSDall = rbind(RMSDall, RMSD)
  }
#  write.table(RMSDall, file=paste0(OUTPUT_RMSD,"2014all.csv"))
  hist(RMSDall[,2], xlab="RMSD", col="lightblue", main=paste0(mainText, "\nNumber of samples: ", length(rownames(RMSDall))) )
}
###############################################################
plotRMSDhistogram <- function(CASFyear=CASF_YEARS[1]) {
  pdf(paste0(RESULT_PATH, "Hist_CASF",CASFyear,".pdf"), width = 12, height = 9)
  par(mfrow = c(1,3))
  RMSDpath = paste0(OUTPUT_RMSD, CASFyear, "/_sampling/")
  #createRMSDhistogram(RMSDpath, pattern="1a", mainText=paste0("CASF",CASFyear," - Set A") )
  createRMSDhistogram(RMSDpath, pattern="sampling_cluster", mainText=paste0("CASF",CASFyear," - Set A") )
  createRMSDhistogram(RMSDpath, pattern="sampling_100"    , mainText=paste0("CASF",CASFyear," - Set B") )
  RMSDpath = paste0(OUTPUT_RMSD, CASFyear, "/_pool/")
  createRMSDhistogram(RMSDpath, pattern="RMSD", mainText=paste0("CASF",CASFyear," - All poses") )
  dev.off()
  
}
###############################################################
# testing 
###############################################################
#classifyMethod = "CASFv2007_sampling_clusters10_RoF_RoT"
#splitWekaResults(RESULT_PATH, classifyMethod, combiName="targets_2016-02-23", TARGET_LIST)
###############################################################
# main part
###############################################################
#splitAllResults()
#splitAllResults_Fidele()
#for (CASFyear in CASF_YEARS) 
#plotRMSDhistogram(CASF_YEARS[4]) 
