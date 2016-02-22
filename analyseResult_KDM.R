###############################################################
# analyse the results for MLSF Project, Feb-2016
# for target KDM
###############################################################
source("libs/libProcessData.R")
source("libs/libValidation.R")
###########################################################################
# DESC: create scatter plot from result
###########################################################################
createPlot <- function(result, evalSet="actives", data="", desc="", method="") {
  x = result[,2]
  y = result[,3]
  plot(x, y, col = "blue", xlab="RMSD", ylab=paste("Predicted Score"))
  prline = lm(y ~ x)  #Linear fit between predicted and measured binding affinity
  abline(prline)     
  new.x = data.frame( x = seq(from = range(x)[1]-2, to = range(x)[2]+2) )
  #lines(x = new.x$x, y = predict( prline, new.x, interval = "confidence" )[ ,"fit" ], col = "red" )
  lines(x = new.x$x, y = predict( prline, new.x, interval = "prediction" )[ ,"upr" ], col = "violet")
  lines(x = new.x$x, y = predict( prline, new.x, interval = "prediction" )[ ,"lwr" ], col = "violet")
  rValue = round(sqrt(calcValidationMetric(x,y)$r2.pred), digits=3)
  title(main=paste("R=", rValue, "on", evalSet, "set (", length(x), "complexes)\n", data, method, desc ,sep=" "))
}
###########################################################################
# DESC: analysis for the actives set
###########################################################################
analyseActivesSet <- function(target, classifyMethod, scorePath=TRUE_SCORES_PATH) {
  resultFile  = paste0(SPLIT_RESULT_PATH, classifyMethod, "_", target, ".csv")
  scoreFile   = paste0(scorePath, target, ".csv")
  result = mergeRowData(resultFile, scoreFile, nameIndex1 = 1, nameIndex2 = 1)  
  
  print(sqrt(calcValidationMetric(result[,2], result[,3])$r2.pred) )
  createPlot(result)
}
#analyseActivesSet("T36_JMJ_Xray_actives", "CASFv2013-refined_sampling_clusters10_RF_", scorePath=PREDICTED_SCORES_PATH)
###########################################################################
# DESC: analysis for the actives set only, by creating scatter plot
# NOTES: not really useful
###########################################################################
analysis_actives_set <- function() {
  #for (desc in descSets) {
  for (CASFset in CASF_SETS) {
    for (trainSet in TRAINING_SETS) {
      for (method in ML_METHODS) {
        #pdf(paste(path2Result,"scatterplot_",protein,"_", trainSet,"_", desc, ".pdf",sep=""), width = 12, height = 9)
        #par(mfrow = c(2,3))
        classifyMethod = paste(CASFset, trainSet, method, sep="_")
        print(classifyMethod)
        analyseActivesSet("T36_JMJ_Xray_actives", classifyMethod, scorePath=PREDICTED_SCORES_PATH)
        #dev.off()      
      }
    }
  }
}
#analysis_actives_set()
###########################################################################
# DESC: merge actives and inactives set to one set, only for KDM target
###########################################################################
createMergeList <- function(classifyMethod, target = "T36_JMJ_Xray", scorePath=PREDICTED_SCORES_PATH) {
  activesFile   = paste0(SPLIT_RESULT_PATH, classifyMethod, "_", target,"_actives.csv")
  scoreFile     = paste0(scorePath, target, "_actives.csv")
  activesList   = mergeRowData(activesFile, scoreFile, nameIndex1 = 1, nameIndex2 = 1)  
  
  inactivesFile = paste0(SPLIT_RESULT_PATH, classifyMethod, "_", target,"_inactives.csv")
  scoreFile     = paste0(scorePath, target, "_inactives.csv")
  inactivesList = mergeRowData(inactivesFile, scoreFile, nameIndex1 = 1, nameIndex2 = 1)  
  
  totalList     = rbind(activesList, inactivesList)
  colnames(totalList)[1] = "ligand"
  return (totalList)
}
#print(createMergeList("CASFv2013-refined_sampling_clusters10_RF_"))
###########################################################################
# DESC: print the actives percent in a mixing list
###########################################################################
calcFoundActivesPercent <- function(activesCutoff, classifyMethod, target = "T36_JMJ_Xray", typeOfScore = "") {
  totalList = createMergeList(classifyMethod, target=target)
  #print(totalList)
  RMSDvalues = totalList[,2]
  scores     = abs(totalList[,3])
  normalizedRMSD    = (RMSDvalues - min(RMSDvalues)) / (max(RMSDvalues) - min(RMSDvalues))
  normalizedScores  = (scores - min(scores)) / (max(scores) - min(scores))
  if (typeOfScore == "") {
    # careful with the order sign, it means sorting ascending
    sortedList = totalList[ order(RMSDvalues),]
  }
  else if (typeOfScore == "original") {
    # careful with the minus order sign, minus means sorting descending
    sortedList = totalList[ order(-scores),]
  }  
  else if (typeOfScore == "NormA") {
    # score type A
    newScores = scores - normalizedRMSD
    # careful with the minus order sign, minus means sorting descending
    sortedList = totalList[ order(-newScores),]
  } else if (typeOfScore == "NormB") {
    # score type B
    newScores = normalizedScores - normalizedRMSD
    # careful with the minus order sign, minus means sorting descending
    sortedList = totalList[ order(-newScores),]
  }
  activesFile   = paste0(SPLIT_RESULT_PATH, classifyMethod, "_", target,"_actives.csv")
  activesNameList = read.csv(activesFile)[,1]

  foundActives = sum( sortedList[1:activesCutoff, "ligand"] %in% activesNameList )
  #foundActives = sum(sortedList[1:activesCutoff, "RMSD"] < activesThreshold)
  foundActivesPercent = foundActives/activesCutoff
  
  return (foundActivesPercent)
}
#print( calcFoundActivesPercent(activesCutoff=53, classifyMethod="CASFv2013-refined_sampling_clusters10_RF_") )

###########################################################################
# DESC: analyse the test set
###########################################################################
analysis_actives_set_all <- function(typeOfScore="") {
    pdf(paste0(RESULT_PATH, "KDM_activespercent",typeOfScore,".pdf"), width = 12, height = 9)
    par(mfrow = c(2,2))
    cutoffList = c(10, 20, 35, 53)
    sapply(cutoffList, buildAnalysisPercentActives, typeOfScore)
    dev.off()
}
###########################################################################
analysis_actives_set_all()
analysis_actives_set_all(typeOfScore="NormA")
analysis_actives_set_all(typeOfScore="NormB")