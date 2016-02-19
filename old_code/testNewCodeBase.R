# NEW CODE BASE FOR MLSF Project
library(randomForest)
library(e1071)
library(gbm)
library(stats)

# return a list of identical cols from data
checkIdenticalData <- function(data) {
  removeIndex = c()
  for (index in seq(1:ncol(data))) {
    if (length(unique(data[,index]))==1) {
      removeIndex = c(removeIndex, index)
    }
  }  
  return (-(removeIndex))
}

# read the result from Weka predict output file (CSV)
readWekaResult <- function(wekaFile)  {
  result = read.csv(wekaFile, skip = 4)
  return (result[,3])  
}

methods = c("RF", "BRT", "SVM", "MLR")
resultPath  = "/home/dat/WORK/output/results/2015-08-26/"
dataPath = "/home/dat/WORK/DB/DESCRIPTORS/August-2015/"
scorePath = "/home/dat/WORK/output/scores/"
#scorePath = "/Users/knight/MyClouds/scores/"

CASF      = c("CASF07", "CASF12", "CASF13", "CASF14") 
descList  = c("elementsv2-SIFt", "elementsv2-SIFt-xscore")
metaClass = c("RoF", "Bagging", "RandSS")
baseClass = c("RoT", "REPT")
scoreTypeList = c("ELscores", "MLscores", "scores")

testPath = paste(dataPath, "CASF07_test_elementsv2-SIFt-xscore.csv", sep="")
trainPath = paste(dataPath, "CASF07_training_elementsv2-SIFt-xscore.csv", sep="")

trainData = read.csv(trainPath, na.strings=c(".", "NA", "", "?"))
testData = read.csv(testPath, na.strings=c(".", "NA", "", "?"))
method = "MLR"
nameIndex = 1
targetIndex = 2

##################################################
# predict the scoring with choosen ML method     #
##################################################
predictMLscore <- function(method = "BRT", trainData, testData, targetIndex = 2, nameIndex = 1) {
  # remove identical col data because SVM can't handle it
  removeIndex = checkIdenticalData(trainData)
  trainData = trainData[,removeIndex]
  testData = testData[,removeIndex]

  # check if the training data intersect with the test data
  if (length(intersect(testData[,nameIndex], trainData[,nameIndex])) > 0) {
    print("##################################################")
    print("## WARNING: intersect between training and test ##")
    print("##################################################")  
    print(intersect(testData[,nameIndex], trainData[,nameIndex]))
  }
  
  # targetIndex: column of target in table
  nTrainData 	= nrow(trainData)		# number of pdb complexes for training
  nTestData 	= nrow(testData) 		# number of pdb complexes for testing
  # get random seed
  seed		= as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
  
  iTrain		= seq(1,nTrainData, 1)
  nSample		= nTrainData
  set.seed(seed)
  selectDesc	= setdiff(c(1:ncol(trainData)), c(nameIndex, targetIndex))
  
  iTrain		= sample(iTrain)[1:nSample]		# shuffle selected complexes
  trainTarget = trainData[ iTrain, targetIndex]
  trainDesc  	= trainData[ iTrain, selectDesc]
  rownames(trainDesc)[1:nrow(trainDesc)] = as.vector(trainData[ iTrain, 1])
  
  iTest		= seq(1,nTestData, 1)
  testTarget  = testData[ iTest, targetIndex] 
  testDesc   	= testData[ iTest, selectDesc]
  rownames(testDesc)[1:nrow(testDesc)]   = as.vector(testData[ iTest, 1])
  
  ##################################################
  # SELECTING ML WITH BEST INTERNAL VALIDATION     #
  ##################################################
  if (method == "SVM") {
    ML_Score 	= svm(trainTarget ~ ., data=trainDesc, na.action=na.omit) 
  } else if (method == "RF") {
    set.seed(seed)
    ML_Score = randomForest(trainTarget ~ ., data=trainDesc, ntree=500, na.action=na.omit)    
  } else if (method == "BRT") {
    ML_Score   = gbm(trainTarget ~ ., data=trainDesc, n.trees=1000)
  } else if (method == "MARS") {
    # \TODO: still doesn't work
    ML_Score   = earth(trainTarget ~ ., data=trainDesc)
  } else if (method == "MLR") {
    ML_Score   = lm(trainTarget ~ ., data=trainDesc)
  } else if (method == "NLS") {
    # \TODO: still doesn't work    
    ML_Score   = nls(trainTarget ~ ., data=trainDesc)
  }

  #######################################################
  # ML-SCORE PREDICTION OF TEST DATASET	       		#
  #######################################################
  
  if (method == "BRT") {
    testPred   = predict(ML_Score, newdata = testDesc, n.trees=1000) 
  } else {
    testPred   = predict(ML_Score, newdata = testDesc) 
  }

#  rmse 		= format(sqrt(mean(( testTarget - testPred)^2)), digits=3) 
#  sdev 		= format(sd( testTarget - testPred ), digits=3)
#  fitpoints 	= na.omit(cbind(testTarget, testPred))
#  fitcorr  	= format(round(cor(fitpoints[,1], fitpoints[,2]), digits=3))   	#Pearson correlation coefficient (R)
#  sprcorr     	= format(cor(rank(fitpoints[,1]), rank(fitpoints[,2])), digits=3) #Spearman rank correlation   
   return (testPred)
}
#predictMLscore(method, trainData, testData)

runML <- function() {
  for (CASFset in CASF ) {
    for (desc in descList) {
      testPath = paste(dataPath, CASFset, "_test_", desc, ".csv", sep="")
      trainPath = paste(dataPath, CASFset, "_training_", desc, ".csv", sep="")
    
      trainData = read.csv(trainPath, na.strings=c(".", "NA", "", "?"))
      testData = read.csv(testPath, na.strings=c(".", "NA", "", "?"))
      
      nTestData   = nrow(testData)   	# number of pdb complexes for testing
      iTest  	= seq(1,nTestData, 1)
      testTarget  = testData[ iTest, 1:2] 
      colnames(testTarget) = c("PDB","experimental")
      
      scoreData = sapply(methods, predictMLscore, trainData, testData)    
      scoreData = cbind(testTarget, scoreData)
      
      scorePath = paste("/home/dat/WORK/output/scores/", CASFset, "_core_MLscores_", desc, ".csv", sep="")
      write.table(scoreData, file=scorePath, sep=",", quote=FALSE, row.names=FALSE)
      
    }  
  }
}
#runML()

runScoreSummaryEL <- function() {
  for (CASFset in CASF) {
    for (desc in descList) {
      testPath = paste(dataPath, CASFset, "_test_", desc, ".csv", sep="")
      testData = read.csv(testPath, na.strings=c(".", "NA", "", "?"))
      
      nTestData   = nrow(testData)     # number of pdb complexes for testing
      iTest    = seq(1,nTestData, 1)
      scoreData  = testData[ iTest, 1:2] 
      colnames(scoreData) = c("PDB","experimental")
      
      for (meta in metaClass) {
        for (base in baseClass) {
          scorePath = paste(resultPath,CASFset,"_test_",meta,"-",base,"_",desc,".csv", sep="")
          #print(scorePath)
          scoreData = cbind(scoreData, readWekaResult(scorePath))
          colnames(scoreData)[ncol(scoreData)] = paste(meta,"-",base,sep="")
        }
      }
      
      scorePath = paste("/home/dat/WORK/output/scores/", CASFset, "_core_ELscores_", desc, ".csv", sep="")
      write.table(scoreData, file=scorePath, sep=",", quote=FALSE, row.names=FALSE)
    }
  }
}

calcCorrelation <- function(scoreName, scoreData, metric="r2.pred") {
  x = scoreData[,"experimental"]
  y = scoreData[,scoreName]
  #print(calcValidationMetric(x,y)[metric])
  rValue = round(sqrt(unlist(calcValidationMetric(x,y)$r2.pred)), digits=3)  
  return (rValue)  
}

calcCorrelationAll <- function(scoreData) {
  scoreList = colnames(scoreData)[c(-1,-2)]
  return (sapply(scoreList, calcCorrelation, scoreData))
  
}

createCorrelList <- function(CASFset) {
  correlList = c()
  for (scoreType in scoreTypeList) {
    if (scoreType == "scores") {
      scoreFile = paste(scorePath,CASFset,"_core_",scoreType,".csv", sep="")
    } else {      
      scoreFile = paste(scorePath,CASFset,"_core_",scoreType,"_",descList[2],".csv", sep="")
    }
    scoreData = read.csv(scoreFile, na.strings=c(".", "NA", "", "?"))
    correl = calcCorrelationAll(scoreData)
    correlList = c(correlList, correl)
    
  }  
  names(correlList) = gsub("REPT", "C4.5", names(correlList))
  
  return (correlList)
}

source("lib/lib.validation_extern.R")

correlData = sapply(CASF, createCorrelList)
correlData = t(correlData)
write.table(correlData, file=paste(scorePath,"Correl_Pearson_core.csv",sep=""), sep=",", quote=FALSE, row.names=TRUE)

createHeatmap <- function(heatmap, title="") {
  op <- par(mar = c(5,9,4,2) + 0.1)
  heatmap = t(heatmap)
  methods = rownames(heatmap)
  methods = methods[c(length(methods):1)]
  targets = colnames(heatmap)
  heatmap = heatmap[nrow(heatmap):1,]
  
  textDesc = ""
  mainLabel = paste(title)
  image(1:length(targets), 1:length(methods), z = t(heatmap), axes = FALSE, xlab = "", ylab = "", main = mainLabel, cex.lab=3, cex.main = 2)
  axis(1, 1:length(targets), targets)
  axis(2, 1:length(methods), methods, labels = FALSE, las = 2)
  at = 1:11
  mtext(side = 2, text = methods[at], at = at, col = "blue", line = 0.8, las = 2)
  at = 12:14
  mtext(side = 2, text = methods[at], at = at, col = "green", line = 0.8, las = 2)
  at = 15:length(methods)
  mtext(side = 2, text = methods[at], at = at, col = "red", line = 0.8, las = 2)
  #legend(0,0, c('Ensemble Learning','Machine Learning','u'), lty=1, col=c('red', 'blue', 'black'), bty='n', cex=.75)
  
  for (x in 1:ncol(heatmap))
    for (y in 1:nrow(heatmap)) {
      text(x, y, round(heatmap[y,x], digits = 3), cex = 1)
    }
  
}

#pdf(paste(scorePath,"scoring.pdf",sep=""))
#createHeatmap(correlData, "Scoring power")
#dev.off()

#rankData = read.csv(paste(scorePath,"Ranking_1xx_core.csv",sep=""))
rankData = read.csv(paste(scorePath,"Ranking_123_core.csv",sep=""))
names(rankData) = gsub("REPT", "C4.5", names(rankData))
#pdf(paste(scorePath,"ranking.pdf",sep=""))
createHeatmap(rankData, "Ranking power")
#dev.off()

testHeatmap <- function() {
for (desc in c("_elementsv2-SIFt_c12b0-xscore", "_elementsv2-SIFt_c12b0")) {
  textDesc = ""
  if (desc == "_elementsv2-SIFt_c12b0-xscore") {
    textDesc = "ESX"
  } else if (desc == "_elementsv2-SIFt_c12b0") {
    textDesc = "ES" }
  pdf(paste(path,"heatmap_spearman_refined_RMSD_", textDesc, ".pdf",sep=""), width = 16, height = 8)
  par(mfrow = c(1,2))
  for (pose in c(1,3)) {
    #pdf(paste(path,"heatmap_spearman_", textDesc, "_Pose",pose, ".pdf",sep=""), width = 14, height = 12)
    createHeatmap(desc, pose)
    #dev.off()
  }
  dev.off()
}
}