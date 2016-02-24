###############################################################
# analyse the results for MLSF Project, Feb-2016
# for target DIG
###############################################################
source("libs/libProcessData.R")
source("libs/libAnalysis.R")
###############################################################
sortAndRemoveSuffix <- function(scorelist, scoreIndex = 3, splitChar = '_', decreasing = TRUE) {
  scorelist[,scoreIndex] = abs(scorelist[,scoreIndex])
  # sorting the list
  if (decreasing) {
    scorelist = scorelist[ order(-scorelist[,scoreIndex]), ]
  }
  else {
    scorelist = scorelist[ order(scorelist[,scoreIndex]), ]
  }
  charPos = regexpr(splitChar, scorelist[, 1]) - 1
  scorelist[,1] = substr(scorelist[, 1], 1, charPos)
  
  return (scorelist)
}
###########################################################################
# ranking based on best pose
###########################################################################
rank_bestpose <- function(scorelist, scoreIndex = 3) {
  poselist = sortAndRemoveSuffix(scorelist, scoreIndex)
  library(hash)
  keep_pose = hash(unique(poselist[,1]), 1) 
  index = 0
  for (pose in poselist[,1]) {
    index = index + 1
    if (keep_pose[[pose]] > 0) {
      keep_pose[[pose]] = keep_pose[[pose]] - 1
    }
    else {
      poselist = poselist[-c(index),]     
      index = index - 1
    }
  }
  poselist = poselist[ order(poselist[,1]), c(1,scoreIndex)]
  return (poselist)
}
###############################################################
calcSpearmanCorRMSD <- function(target, classifyMethod, numberOfPoses = 1) {
  resultFile = paste(SPLIT_RESULT_PATH, classifyMethod, "_", target, ".csv", sep="")
  #print(resultFile)
  scoreFile = paste(PREDICTED_SCORES_PATH, target, ".csv", sep="")
  #print(scoreFile)
  result = mergeRowData(resultFile, scoreFile, nameIndex1 = 1, nameIndex2 = 1)
  #print(result)
  #RMSDvalues = result[,2]
  #scores     = abs(result[,3])
  #normalizedRMSD    = (RMSDvalues - min(RMSDvalues)) / (max(RMSDvalues) - min(RMSDvalues))
  #normalizedScores  = (scores - min(scores)) / (max(scores) - min(scores))
  # score type A
  #newScores = scores - normalizedRMSD
  # score type B
  #newScores = normalizedScores - normalizedRMSD
  #result[,3] = newScores
  #topPoses = rank_bestpose(result, scoreIndex = 3)
  
  topPoses = calc_average_bestpose(list_top_bestpose(result, numberOfPoses, IDindex = 1, decreasing = FALSE), scoreIndex = 3, IDindex = 1) 
  topPoses[,2] = abs(topPoses[,2])    
  rankTopPoses = rank(-topPoses[,2])
  return ( round(cor(TRUE_RANK_DIG10.2, rankTopPoses, method="spearman"), digits = 3) )
}
###############################################################
analyse_result_MLscoring <- function(CASFset, trainSet, numberOfPoses = 3) {
  result = data.frame(row.names=ML_METHODS)
  for (method in ML_METHODS) {
    suffix = "_RMSD_DIG10.2"  
    classifyMethod = paste(CASFset, trainSet, method, sep="_")
    #print(classifyMethod)
    #IDFile = paste("DIG10.2_", desc, sep="")
    test = sapply(TARGET_LIST[6:8], calcSpearmanCorRMSD, classifyMethod, numberOfPoses = numberOfPoses)
    result = rbind(result, test)
  }
  colnames(result) = TARGET_LIST[6:8]
  rownames(result) = paste(CASFset, trainSet, ML_METHODS, sep="_")
  return (result)
}
###############################################################
runCalcSpearman <- function() {
#  for (desc in DESC_SETS) {
  
    for (CASFset in CASF_SETS) {
      for (trainSet in TRAINING_SETS) {
        for (poses in c(1)) {
          endResult = analyse_result_MLscoring(CASFset, trainSet, numberOfPoses=poses) 
          write.table(endResult, file = paste0("Spearman_", CASFset, trainSet, "_DIG_Pose_", poses, ".csv"), sep=",")
        }
      }
    }
#  }
}
###############################################################
analyse_result_MLscoring_new <- function(CASFset, numberOfPoses = 1) {
  result = data.frame(row.names=ML_METHODS)
  for (trainSet in TRAINING_SETS) 
    for (method in ML_METHODS) {
      suffix = "_RMSD_DIG10.2"  
      classifyMethod = paste(CASFset, trainSet, method, sep="_")
      #print(classifyMethod)
      #IDFile = paste("DIG10.2_", desc, sep="")
      test = sapply(TARGET_LIST[6:8], calcSpearmanCorRMSD, classifyMethod, numberOfPoses = numberOfPoses)
      result = rbind(result, test)
    }
  colnames(result) = TARGET_LIST[6:8]
  # for shorter name
  rownames(result) = sapply(TRAINING_SETS_NAME, paste, ML_METHODS, sep="_")
  
  return (result)}
###############################################################
runCalcSpearman_new <- function() {
  #  for (desc in DESC_SETS)
  for (CASFset in CASF_SETS)
      for (poses in c(1)) {
        endResult = analyse_result_MLscoring_new(CASFset, numberOfPoses=poses) 
        write.table(endResult, file = paste0(RESULT_PATH, "DIG_Spearman_", CASFset, ".csv"), sep=",")
      }
}
###############################################################
plotHeatmap <- function(CASFset) {
  dataSetsShortName = sapply(TRAINING_SETS_NAME, paste, ML_METHODS, sep="_")
  dataSetsShortName = c("original scores", dataSetsShortName)

  heatmap = read.csv(paste0(RESULT_PATH, "DIG_Spearman_", CASFset,".csv"))
  heatmap = rbind("original scores" = c(0.672, 0.66, .551), heatmap)
  #rownames(heatmap)[length(dataSetsShortName)] = 
  image(1:ncol(heatmap), 1:length(dataSetsShortName), z = t(heatmap), axes = FALSE, xlab = "", ylab = "", main= CASFset, cex.lab=2, cex.main = 2)
  axis(1, 1:ncol(heatmap), colnames(heatmap))
  axis(2, 1:length(dataSetsShortName), rownames(heatmap), labels=FALSE, las=2)
  at = 1:length(dataSetsShortName)
  mtext(side = 2, text = dataSetsShortName[at], at = 1:length(dataSetsShortName), line=0.5, las=2, cex=0.7)
  
  for (x in 1:ncol(heatmap))
    for (y in 1:nrow(heatmap))
      text(x, y, round(heatmap[y,x], digits = 4)*100, cex=1)
  
}
###############################################################
#analyseDIG_test <- function() {
  test = read.csv(file = paste0(SPLIT_RESULT_PATH, "CASFv2014-refined_sampling_100_RoF_RoT_DIG_SP.csv") )
  score = read.csv(file = paste0(PREDICTED_SCORES_PATH, "DIG_SP.csv") )
  test = merge(test, score, by.x=1, by.y=1)
  #print(test)
  RMSDvalues = test[,2]
  scores     = abs(test[,3])
  normalizedRMSD    = (RMSDvalues - min(RMSDvalues)) / (max(RMSDvalues) - min(RMSDvalues))
  normalizedScores  = (scores - min(scores)) / (max(scores) - min(scores))
  # score type A
  #newScores = scores - normalizedRMSD
  # score type B
  #newScores = normalizedScores - normalizedRMSD
  #test[,3] = newScores
  rankpose = rank_bestpose(test)
  # sorting after index
  #rankpose = rankpose[ order(rankpose[,1]), c(1,3)]
  print(rankpose)
  print(rank(-rankpose[,2]))
  print(cor(TRUE_RANK_DIG10.2, rank(-rankpose[,2]), method="spearman"))
  
  #rankpose = list_top_bestpose(test, numberOfPoses = 1, IDindex = 1, decreasing = FALSE)
  #print(rankpose)
  rankpose = calc_average_bestpose(list_top_bestpose(test, numberOfPoses = 1, IDindex = 1, decreasing = FALSE), scoreIndex = 3, IDindex = 1)
  rankpose[,2] = abs(rankpose[,2])
  print(rankpose)
  print(rank(-rankpose[,2]))
  print(cor(TRUE_RANK_DIG10.2, rank(-rankpose[,2]), method="spearman"))
#}
#analyseDIG_test()
plotAll <- function() {
  pdf(paste0(RESULT_PATH, "DIG_Spearman_1.pdf"), width = 12, height = 9)
  par(mfrow = c(1,2))
  plotHeatmap(CASF_SETS[1])
  plotHeatmap(CASF_SETS[2])
  dev.off()
  pdf(paste0(RESULT_PATH, "DIG_Spearman_2.pdf"), width = 12, height = 9)
  par(mfrow = c(1,2))
  plotHeatmap(CASF_SETS[3])
  plotHeatmap(CASF_SETS[4 ])
  dev.off()
}
#plotAll()
###############################################################

###############################################################
#print(calcSpearmanCorRMSD("DIG_SP", "CASFv2014-refined_sampling_100_RoF_RoT"))
#print(analyse_result_MLscoring(CASF_SETS[2], TRAINING_SETS[1]))
#runCalcSpearman_new()

