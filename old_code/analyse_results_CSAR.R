source("old_code/lib/configuration.R")
source("old_code/lib/lib.analysis.R")

resultPath = "/home/dat/WORK/RESULTS/2015-06-23/"

IDPath = "/home/dat/WORK/DB/DESCRIPTORS/June-2015/process/"

###############################################################
calcSpearmanCor <- function(target, resultFile, IDFile, numberOfPoses = 3) {
  resultFile = paste(resultFile, "_", target, ".csv", sep="")
  #print(resultFile)
  result = readWekaResult(resultPath, resultFile)
  IDFile = paste(IDFile, target, "process.csv", sep="_")
  #print(IDFile)
  result = addID(result, IDPath, IDFile)

  topPoses = calc_average_bestpose(list_top_bestpose(result, numberOfPoses, IDindex = 3), IDindex = 3) 
  topPoses[,2] = abs(topPoses[,2])
  rankTopPoses = rank(-topPoses[,2])
  return (cor(TRUE_RANK_CSAR13, rankTopPoses, method="spearman"))
}
###############################################################
calcSpearmanCorRMSD <- function(target, resultFile, IDFile, numberOfPoses = 3) {
  resultFile = paste(resultFile, "_", target, ".csv", sep="")
  #print(resultFile)
  result = readWekaResult(resultPath, resultFile)
  IDFile = paste(IDFile, target, "process.csv", sep="_")
  #print(IDFile)
  result = addID(result, IDPath, IDFile)
  
  topPoses = calc_average_bestpose(list_top_bestpose(result, numberOfPoses, IDindex = 3, decreasing = FALSE), scoreIndex = 1, IDindex = 3) 
  topPoses[,2] = abs(topPoses[,2])    
  rankTopPoses = rank(-topPoses[,2])
  return (cor(TRUE_RANK_CSAR13, rankTopPoses, method="spearman"))
}
###############################################################
# with replate \NOT WORKING
calcSpearmanCorRMSDWithReplace <- function(target, resultFile, IDFile, numberOfPoses = 3) {
  resultFile = paste(resultFile, "_", target, ".csv", sep="")
  predictFile = "CASF12_refined_DIG10.2_RoF-REP_elementsv2-SIFt_c12b0_XP.csv"
  
  predict = readWekaResult(resultPath, predictFile)
  test = readWekaResult(resultPath, resultFile)
  # replace the actual with the predicted by same score
  #test$actual = predict$predicted
  #print(resultFile)
  result = readWekaResult(resultPath, resultFile)
  IDFile = paste(IDFile, target, "process.csv", sep="_")
  #print(IDFile)
  result = addID(result, IDPath, IDFile)
  
  topPoses = calc_average_bestpose(list_top_bestpose(result, numberOfPoses, IDindex = 3, decreasing = FALSE), scoreIndex = 1, IDindex = 3) 
  topPoses[,2] = abs(topPoses[,2])    
  rankTopPoses = rank(-topPoses[,2])
  return (cor(TRUE_RANK_CSAR13, rankTopPoses, method="spearman"))
}
###############################################################
analyse_result_MLscoring <- function(trainingData, desc, numberOfPoses = 1, RMSD=FALSE) {
  result = data.frame(row.names=dockingMethods)
  for (trainingMethod in methods) {
    if (RMSD) {
      suffix = "refined_RMSD_DIG10.2"  
    } else {
      suffix = "refined_DIG10.2"  
    }
    
    resultFile = paste(trainingData, suffix, trainingMethod, desc, sep="_")
    print(resultFile)
    IDFile = paste("DIG10.2_", desc, sep="")
    if (!RMSD) {
      test = sapply(dockingMethods, calcSpearmanCor, resultFile, IDFile, numberOfPoses = numberOfPoses)
    }
    else {
      test = sapply(dockingMethods, calcSpearmanCorRMSD, resultFile, IDFile, numberOfPoses = numberOfPoses)
    }
    result = rbind(result, test)
  }
  colnames(result) = dockingMethods
  rownames(result) = paste(trainingData, methods, desc, sep="_")
  return (result)
}
###############################################################
runCalcSpearman <- function(RMSD=FALSE) {
  for (desc in descSets) {
    for (poses in c(1,3)) {
      if (!RMSD) {
        endResult = lapply(dataSets, analyse_result_MLscoring, desc=desc, numberOfPoses=poses, RMSD)  
        endResult = do.call(rbind.data.frame, endResult)        
      } else {
        endResult = analyse_result_MLscoring("CASF14", desc=desc, numberOfPoses=poses, RMSD) 
      }
      
      if (RMSD) {
        write.table(endResult, file = paste("Spearman_refined_RMSD_DIG10.2_", desc, "_Pose_", poses, ".csv",sep=""), sep=",")
      }
      else {
        write.table(endResult, file = paste("Spearman_DIG10.2_", desc, "_Pose_", poses, ".csv",sep=""), sep=",")  
      }
      
    }
  }
}
###############################################################
#heatmap = heatmap[nrow(heatmap):1,]

#mainLabel = paste(metricName, "on", evalSet, "set")
#image(1:6, 1:6, z = t(heatmap), axes = FALSE, xlab = paste("Descriptor:"), ylab = "", main = "main", cex.lab=2, cex.main = 2)
#axis(1, 1:length(methods), colnames(heatmap))
#axis(2, 1:length(dataSets), rownames(heatmap))
#for (x in 1:ncol(heatmap))
#  for (y in 1:nrow(heatmap))
#    text(x, y+0.4, round(heatmap[y,x], digits = 3), cex = 1)
###############################################################
# ranking based on average of numberOfPoses best poses
analyse_result_original_poses <- function(numberOfPoses) {
  for (method in dockingMethods) {
    print(method)
    poseList = read.csv(paste(IDPath, "cs-confgen_all_", method, ".csv", sep=""))
    topPoses = calc_average_bestpose(list_top_bestpose(poseList, numberOfPoses, IDindex = 1), IDindex = 1) 
    #print(topPoses)
    rankTopPoses = rank(-topPoses[,2])
    print(rankTopPoses)
    print(cor(TRUE_RANK_CSAR13, rankTopPoses, method="spearman"))
  }
}
#analyse_result_original_poses(1)
#analyse_result_original_poses(3)
#runCalcSpearman(RMSD=TRUE)

testing <- function(target) {
  resultFile = paste("CASF12_refined_DIG10.2_RoF-REP_elementsv2-SIFt_c12b0_", target, ".csv", sep="")
  test = readWekaResult(resultPath, resultFile)
  #print(resultFile)
  IDFile = paste("DIG10.2_elementsv2-SIFt_c12b0", target, "process.csv", sep="_")
  #print(IDFile)
  test = addID(test, IDPath, IDFile)
  
  #print(rank_bestpose(test, IDindex = 3))
  #x = list_top_bestpose(test, numberOfPoses = 1, IDindex = 3)
  
  rankpose = calc_average_bestpose(list_top_bestpose(test, numberOfPoses = 3, IDindex = 3), scoreIndex = 2, IDindex = 3)
  rankpose[,2] = abs(rankpose[,2])
  print(rankpose)
  print(rank(-rankpose[,2]))
  print(cor(TRUE_RANK_CSAR13, rank(-rankpose[,2]), method="spearman"))
  
}

for (target in dockingMethods) {
  #testing(target)
}

testingRMSD <- function() {

  resultFile = "CASF14_refined_RMSD_DIG10.2_RoF-RT_elementsv2-SIFt_c12b0_SP.csv"

  IDFile = "DIG10.2_elementsv2-SIFt_c12b0_SP_process.csv"
  
  predictFile = "CASF14_refined_DIG10.2_RoF-RT_elementsv2-SIFt_c12b0_SP.csv"
  
  predict = readWekaResult(resultPath, predictFile)
  test = readWekaResult(resultPath, resultFile)
  # replace the actual with the predicted by same score
  #test$actual = predict$predicted
  test = addID(test, IDPath, IDFile)
  
  #print(test)
  #print(rank_bestpose(test, IDindex = 3))
  #x = list_top_bestpose(test, numberOfPoses = 1, IDindex = 3, decreasing = FALSE)
  #print(x)
  #x = x[order(x$ID),]
  #x = x[c(3,2)]
  #print(x)
  
  rankpose = list_top_bestpose(test, numberOfPoses = 3, IDindex = 3, decreasing = FALSE)
  print(rankpose)
  rankpose = calc_average_bestpose(list_top_bestpose(test, numberOfPoses = 3, IDindex = 3, decreasing = FALSE), scoreIndex = 1, IDindex = 3)
  rankpose[,2] = abs(rankpose[,2])
  print(rankpose)
  print(rank(-rankpose[,2]))
  print(cor(TRUE_RANK_CSAR13, rank(-rankpose[,2]), method="spearman"))
  
}
#testingRMSD()
