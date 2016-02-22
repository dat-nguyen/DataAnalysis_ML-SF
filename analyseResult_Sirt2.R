###############################################################
# analyse the results for MLSF Project, Feb-2016
# for target HDAC
###############################################################
source("libs/libProcessData.R")
source("libs/libValidation.R")

###########################################################################
# DESC: print the actives percent in a mixing list
###########################################################################
calcFoundActivesPercent <- function(activesCutoff, classifyMethod, target=TARGET_LIST[1], typeOfScore = "") {
  splitChar = '_'
  
  ligandList = read.csv(paste0(PREDICTED_SCORES_PATH, "Sirt2.txt"), header=FALSE)
  # modify the ligand name because of splitting problem
  ligandList[,1] = gsub(' $', '', ligandList[,1])
  ligandList[,1] = gsub(' ', splitChar, ligandList[,1])
  
  scoreFile   = paste0(SPLIT_RESULT_PATH, classifyMethod, "_", target,".csv")
  scoreList   = read.csv(scoreFile, na.strings=c(".", "NA", "", "?"))
  # modify the ligand name because of splitting problem
  scoreList[,1] = gsub(pattern=paste0(splitChar,"[0-9]+$"), replacement="", scoreList[,1])
  
  ligandList    = merge(ligandList, scoreList, by.x = 1, by.y = 1)
  activesIndex  = (ligandList[,2] > 5)
  activesList   = na.omit(ligandList[activesIndex, c(1,5,3)])
  totalList     = ligandList[, c(1,5,3)]
  colnames(totalList)[1] = "ligand"
  activesNameList = activesList[,1]
    
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
  #print(totalList)
  foundActives = sum( sortedList[1:activesCutoff, "ligand"] %in% activesNameList )
  #foundActives = sum(sortedList[1:activesCutoff, "RMSD"] < activesThreshold)
  foundActivesPercent = foundActives/activesCutoff
  
  return (foundActivesPercent)
}
#print( calcFoundActivesPercent(activesCutoff=93, classifyMethod="CASFv2013-refined_sampling_clusters10_RF_") )
###########################################################################
# DESC:
###########################################################################
buildAnalysisPercentActives <- function(activesCutoff, typeOfScore="") {
  dataSets = sapply(TRAINING_SETS, paste, ML_METHODS, sep="_")
  dataSets = c(dataSets, "original scores")
  foundActives = matrix(NA, nrow = length(dataSets), ncol = length(CASF_SETS))
  rownames(foundActives) = dataSets
  colnames(foundActives) = CASF_SETS
  for (i in seq(length(dataSets))) {
    for (j in seq(length(CASF_SETS))) {
      if (i == length(dataSets)) {
        classifyMethod = paste(CASF_SETS[j], dataSets[i-1], sep="_")
        foundActives[i,j] = calcFoundActivesPercent(activesCutoff, classifyMethod, typeOfScore="original")
      }
      else {
        classifyMethod = paste(CASF_SETS[j], dataSets[i], sep="_")
        foundActives[i,j] = calcFoundActivesPercent(activesCutoff, classifyMethod, typeOfScore=typeOfScore)
      }
    }
  }
  heatmap = foundActives[nrow(foundActives):1,]
  image(1:length(CASF_SETS), 1:length(dataSets), z = t(heatmap), axes = FALSE, xlab = "", ylab = "", main = activesCutoff, cex.lab=2, cex.main = 2)
  axis(1, 1:length(CASF_SETS), colnames(heatmap))
  axis(2, 1:length(dataSets), rownames(heatmap), labels=FALSE, las=2)
  at = 1:length(dataSets)
  dataSetsShortName = sapply(TRAINING_SETS_NAME, paste, ML_METHODS, sep="_")
  dataSetsShortName = c(dataSetsShortName, "original scores")
  mtext(side = 2, text = dataSetsShortName[at], at = length(dataSets):1, line=0.5, las=2, cex=0.7)
  
  for (x in 1:ncol(heatmap))
    for (y in 1:nrow(heatmap))
      text(x, y, round(heatmap[y,x], digits = 4)*100, cex=1)
  
}
#buildAnalysisPercentActives(93, typeOfScore="")
###########################################################################
# DESC: analyse the test set
###########################################################################
analysis_test_set_all <- function(typeOfScore="") {
    #pdf(paste0(RESULT_PATH, "Sirt2_activespercent",typeOfScore,".pdf"), width = 12, height = 9)
    par(mfrow = c(1,2))
    cutoffList = c(10, 22)
    sapply(cutoffList, buildAnalysisPercentActives, typeOfScore)
    #dev.off()
}
###########################################################################
#test <- function() {
  classifyMethod = "CASFv2013-refined_sampling_clusters10_RF_"
  target = TARGET_LIST[1]
  splitChar = '_'

  ligandList = read.csv(paste0(PREDICTED_SCORES_PATH, "Sirt2.txt"), header=FALSE)
  # modify the ligand name because of splitting problem
  ligandList[,1] = gsub(' $', '', ligandList[,1])
  ligandList[,1] = gsub(' ', splitChar, ligandList[,1])
  
  scoreFile   = paste0(SPLIT_RESULT_PATH, classifyMethod, "_", target,".csv")
  scoreList   = read.csv(scoreFile, na.strings=c(".", "NA", "", "?"))
  # modify the ligand name because of splitting problem
  scoreList[,1] = gsub(pattern=paste0(splitChar,"[0-9]+$"), replacement="", scoreList[,1])
  
  ligandList    = merge(ligandList, scoreList, by.x = 1, by.y = 1)
  activesIndex  = (ligandList[,2] > 5)
  activesList   = na.omit(ligandList[activesIndex, c(1,5,3)])
  totalList     = ligandList[, c(1,5,3)]
  colnames(totalList)[1] = "ligand"
  activesNameList = activesList[,1]
  
  #print(ligandList[,1] %in% scoreList[,1])
  
#}
###########################################################################
analysis_test_set_all()
analysis_test_set_all(typeOfScore="NormA")
analysis_test_set_all(typeOfScore="NormB")