###############################################################
# analyse the results for MLSF Project, Feb-2016
###############################################################
source("libs/libProcessData.R")
matchPattern = c("Aminothiazoles_71_PB_md1_rst_Dat", "SmHDAC8_inhibitors", 
                 TEST_SETS_KDM, 
                 TEST_SETS_DIG)

###############################################################
splitAllResults <- function() {
  for (CASFset in CASF_SETS) {
    for (trainSet in TRAINING_SETS) {
      for (method in ML_METHODS) {
        classifyMethod = paste(CASFset, "_", trainSet, "_", method, sep="")
        splitWekaResults(RESULT_PATH, classifyMethod, combiName="targets_2016-02-23", matchPattern)
      }
    }
  } 
}
#splitAllResults()
###############################################################
analyseDIG_test <- function() {
  test = read.csv(file = "/home/dat/WORK/RESULTS/2016-02-23/splitted_results/CASFv2013-refined_sampling_100_RF__DIG_SP.csv")
  score = read.csv(file = "/home/dat/WORK/RESULTS/2016-02-23/true_scores/DIG_SP.csv")
  test = merge(test, score, by.x=1, by.y=1)
  #print(test)
  rankpose = list_top_bestpose(test, numberOfPoses = 3, IDindex = 1, decreasing = FALSE)
  print(rankpose)
  
  rankpose = calc_average_bestpose(list_top_bestpose(test, numberOfPoses = 3, IDindex = 1, decreasing = FALSE), scoreIndex = 3, IDindex = 1)
  rankpose[,2] = abs(rankpose[,2])
  print(rankpose)
  print(rank(-rankpose[,2]))
  print(cor(TRUE_RANK_CSAR13, rank(-rankpose[,2]), method="spearman"))
}
###############################################################
calcSpearmanCorRMSD <- function(target, classifyMethod, numberOfPoses = 1) {
  resultFile = paste(SPLIT_RESULT_PATH, classifyMethod, "_", target, ".csv", sep="")
  print(resultFile)
  scoreFile = paste(TRUE_SCORES_PATH, target, ".csv", sep="")
  print(scoreFile)
  
  result = mergeRowData(resultFile, scoreFile, nameIndex1 = 1, nameIndex2 = 1)
  topPoses = calc_average_bestpose(list_top_bestpose(result, numberOfPoses, IDindex = 1, decreasing = FALSE), scoreIndex = 3, IDindex = 1) 
  topPoses[,2] = abs(topPoses[,2])    
  rankTopPoses = rank(-topPoses[,2])
  return (cor(TRUE_RANK_DIG10.2, rankTopPoses, method="spearman"))
}

#print(calcSpearmanCorRMSD("DIG_SP", "CASFv2013-refined_sampling_clusters10_RoF_RoT"))
###############################################################
analyse_result_MLscoring <- function(CASFset, trainSet, numberOfPoses = 3) {
  result = data.frame(row.names=dockingMethods)
  for (method in ML_METHODS) {
    suffix = "_RMSD_DIG10.2"  
    classifyMethod = paste(CASFset, "_", trainSet, "_", method, sep="")
    print(classifyMethod)
    #IDFile = paste("DIG10.2_", desc, sep="")
    test = sapply(matchPattern[6:8], calcSpearmanCorRMSD, classifyMethod, numberOfPoses = numberOfPoses)
    result = rbind(result, test)
  }
  colnames(result) = matchPattern[6:8]
  #rownames(result) = paste(trainingData, methods, desc, sep="_")
  return (result)
}
print(analyse_result_MLscoring(CASF_SETS[2], TRAINING_SETS[1]))
###############################################################
runCalcSpearman <- function() {
#  for (desc in DESC_SETS) {
    for (CASFset in CASF_SETS) {
      for (trainSet in TRAINING_SETS) {
        for (poses in c(1,3)) {
          endResult = analyse_result_MLscoring(CASFset, trainSet, numberOfPoses=poses) 
          write.table(endResult, file = paste("Spearman_refined_RMSD_DIG10.2_", desc, "_Pose_", poses, ".csv",sep=""), sep=",")
        }
      }
    }
#  }
}
###############################################################
#classifyMethod = "CASFv2007_sampling_clusters10_RoF_RoT"
#classifyMethod = "CASFv2013-refined_sampling_clusters10_RoF_RoT"
#classifyMethod = "CASFv2013-refined_sampling_100_RoF_RoT"
#classifyMethod = "CASFv2013-refined_sampling_100_RF_"
#splitWekaResults(RESULT_PATH, classifyMethod, combiName="targets_2016-02-23", matchPattern)