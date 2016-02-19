library(ROCR);
source("lib/lib.analysis.R")

DUDE_sets = c("ADA", "RENI", "FGFR1")
# calculate AUC and EF of a typical VS assertment
calcVSassertment <- function(result_pose, actives) {
  result_pose$IsActive <- as.numeric(result_pose$ID %in% actives)
  
  pred <- prediction(result_pose$predicted, result_pose$IsActive)
  perf <- performance(pred, 'tpr','fpr')
  auc_rdock <- performance(pred, "auc")
  auc.area_rdock <- slot(auc_rdock, "y.values")[[1]]
  
  EF_rdock <- perf@y.values[[1]]/perf@x.values[[1]]
  EF_rdock_1 <- EF_rdock[which(perf@x.values[[1]] > 0.01)[1]]
  EF_rdock_20 <- EF_rdock[which(perf@x.values[[1]] > 0.2)[1]]  
  
  return (list("AUC" = auc.area_rdock, "EF_1" = EF_rdock_1, "EF_20" = EF_rdock_20))
}

#
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

#"/Users/knight/Dropbox/data/DUD-E/"
activesPath = "/home/dat/WORK/DB/DUD-E/"

#resultPath = "/Users/knight/Dropbox/data/"
resultPath = "/home/dat/WORK/output/results/2015-06-23/DUD-E/"

#IDPath = "/Users/knight/Dropbox/data/"
IDPath = "/home/dat/WORK/DB/DESCRIPTORS/DUD-E/"

IDFile = "DUD-E_SIFt_c12b0_ADA.csv"
#IDFile = "DIG10.2_elementsv2-SIFt_c12b0_asp_process.csv"

resultFile = "CASF12_refined_DUD-E_RoF-REP_elementsv2-SIFt_c12b0_ADA.csv"
activesFile = "ADA_actives.txt"


#load ligands and decoys
#dec <- unique(read.table("decoys.txt")[,1]);
#analyse_VS <- function(activesFile, resultFile, IDFile) {
  actives <- unique(read.table(paste(activesPath, activesFile, sep=""))[,1])
  
  result = readWekaResult(resultPath, resultFile)
  result = addID(result, IDPath, IDFile, IDIndex = "last")
  result_bestpose = list_top_bestpose(result, numberOfPoses = 1, IDindex = 3)
  #result_bestpose = result_bestpose[order(result_bestpose$ID),]
  result_bestpose = result_bestpose[c(3,2)]
  bestpose = calcVSassertment(result_bestpose, actives)
  
  result_bestavgPose = calc_average_bestpose(list_top_bestpose(result, IDindex = 3), IDindex = 3)  
  #result_bestaveragepose$IsActive <- as.numeric(result_bestaveragepose$ID %in% actives)
  colnames(result_bestavgPose) = c("ID", "predicted")
  bestavgPose= calcVSassertment(result_bestavgPose, actives)
  
  output = cbind(do.call(rbind.data.frame, bestpose), do.call(rbind.data.frame, bestavgPose))
  rownames(output) = c("bestpose", "bestavgPose")
  print(output)
  #return (output)
#}

analyse_result_MLscoring <- function(trainingData, desc, numberOfPoses = 1, RMSD=FALSE) {
  result = data.frame(row.names=dockingMethods)
  for (trainingMethod in methods) {
    if (RMSD) {
      suffix = "refined_RMSD_DUD-E"  
    } else {
      suffix = "refined_DUD-E"  
    }
    
    resultFile = paste(trainingData, suffix, trainingMethod, desc, sep="_")
    print(resultFile)
    IDFile = paste("DUD-E_", desc, sep="")
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

runCalcVS <- function(RMSD=FALSE) {
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
run_analyse_DUDE <- function() {
  for (desc in descSets) {
    for (trainingData in dataSets) {      
      for (trainingMethod in methods) {
        for (target in DUDE_sets) {
          suffix = "refined_DUD-E"
          print(paste(trainingData, suffix, trainingMethod, desc, target, sep="_"))
          #print(paste("DUD-E_SIFt_c12b0_", target, ".csv", sep=""))
          #print(paste(target, "_actives.txt", sep=""))
        }
      }
    }
  }
}

#run_analyse_DUDE()
#x = analyse_VS(activesFile, resultFile, IDFile)

#print(cor(TRUE_RANK_CSAR13, rank(-rankpose[,2]), method="spearman"))
#poseList = read.csv(paste("/Users/knight/Dropbox/data/DUD-E/", "DUD-E_SIFt_c12b0_ADA.csv", sep=""))
#len = length(colnames(poseList))
#print(strsplit(toString(poseList[,1]), '_'))
#post = regexpr('_', poseList[, len]) - 1
#print(unique(substr(poseList[, len],1,post)))
