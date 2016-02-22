###############################################################
# analyse the results for MLSF Project, Feb-2016
# for target HDAC
###############################################################
source("libs/libProcessData.R")
###########################################################################
#test <- function() {
  classifyMethod = "CASFv2014-refined_sampling_clusters10_RF_"
  print(FideleList)
  target = FideleList[22]
  # for Glide score
  splitChar = '_'

  ligandList = read.csv(paste0(PREDICTED_SCORES_PATH, target, ".txt"), header=FALSE)
  
  scoreFile   = paste0(SPLIT_RESULT_PATH, "Fidele/", classifyMethod, "_", target,".csv")
  scoreList   = read.csv(scoreFile, na.strings=c(".", "NA", "", "?"))
  # modify the ligand name because of splitting problem
  ligandList[,1] = scoreList[,1] 
  
  ligandList    = merge(ligandList, scoreList, by.x = 1, by.y = 1)
  colnames(ligandList)[1] = "ligand"

  print(sum(ligandList[,1] %in% scoreList[,1]))
  #print(sum(scoreList[,1] %in% ligandList[,1]))
  #print( scoreList[!(scoreList[,1] %in% ligandList[,1]),] )
  print(length(ligandList[,1]))
  #print(cor(ligandList[,2], ligandList[,3], method="spearman"))
  print(cor(abs(ligandList[,2]), ligandList[,3], method="pearson"))
  #hist(ligandList[,3])
test <- function() {  
  # for Goldscore
  target = FideleList[21]
  ligandList = read.csv(paste0(PREDICTED_SCORES_PATH, target, ".txt"), header=FALSE)
  # modify the ligand name because of splitting problem
  ligandList[,1] = gsub(' $', '', ligandList[,1])
  ligandList[,1] = gsub(' ', splitChar, ligandList[,1])
  
  scoreFile   = paste0(SPLIT_RESULT_PATH, "Fidele/", classifyMethod, "_", target,".csv")
  scoreList   = read.csv(scoreFile, na.strings=c(".", "NA", "", "?"))
  # modify the ligand name because of splitting problem
  scoreList[,1] = gsub(pattern=paste0(splitChar,"[0-9]+$"), replacement="", scoreList[,1])
  
  ligandList    = merge(ligandList, scoreList, by.x = 1, by.y = 1)
  colnames(totalList)[1] = "ligand"
  
  print(sum(ligandList[,1] %in% scoreList[,1]))
  #print(sum(scoreList[,1] %in% ligandList[,1]))
  #print( scoreList[!(scoreList[,1] %in% ligandList[,1]),] )
  print(length(ligandList[,1]))
  #print(cor(ligandList[,2], ligandList[,3], method="spearman"))
  print(cor(ligandList[,2], ligandList[,3], method="pearson"))
}
test()
###########################################################################
#analysis_test_set_all()
#analysis_test_set_all(typeOfScore="NormA")
#analysis_test_set_all(typeOfScore="NormB")