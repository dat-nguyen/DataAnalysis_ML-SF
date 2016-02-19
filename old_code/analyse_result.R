source("lib/lib.validation_extern.R")
source("lib/lib.process_data.R")
source("lib/lib.analysis.R")
source("lib/configuration.R")

###########################################################################
# DESC: create scatter plot from result
###########################################################################
createPlot <- function(result, evalSet="actives", data="", desc="", method="") {
  x = result[,1]
  y = result[,2]
  plot(x, y, col = "blue", xlab="Measured binding affinity", ylab=paste("Predicted binding affinity"))
  prline = lm(y ~ x)  #Linear fit between predicted and measured binding affinity
  abline(prline)     
  new.x = data.frame( x = seq(from = range(x)[1]-2, to = range(x)[2]+2) )
  #lines(x = new.x$x, y = predict( prline, new.x, interval = "confidence" )[ ,"fit" ], col = "red" )
  lines(x = new.x$x, y = predict( prline, new.x, interval = "prediction" )[ ,"upr" ], col = "violet")
  lines(x = new.x$x, y = predict( prline, new.x, interval = "prediction" )[ ,"lwr" ], col = "violet")
  rValue = round(sqrt(calcValidationMetric(x,y)$r2.pred), digits=3)
  title(main=paste("R=", rValue, "on", evalSet, "set (", length(x), "complexes)\n", data, method, desc ,sep=" "))
}

# DESC: reanalysis the CASF 
analysis_CASF <- function() {
  for (desc in descSets) {
    for (evalSet in evalSets) {
      pdf(paste(path2Result,"heatmap_",desc, "_",evalSet, ".pdf",sep=""), width = 14, height = 8)
      par(mfrow = c(1,3))
      sapply(c("r2.pearson", "r2m", "rse"), createHeatmap, dataSets, methods, resultPath, evalSet, desc)
      dev.off()
      pdf(paste(path2Result,"scatterplot_",desc, "_",evalSet, ".pdf",sep=""), width = 12, height = 9)
      par(mfrow = c(2,3))
      for (method in methods) {
        for (data in dataSets) {      
          result = readClassifyResult(resultPath, data, evalSet, method, desc)
          createPlot(result, evalSet, data, desc, method)
        }
      }
      dev.off()
    }
  }
}

#analysis_CASF()

#createHeatmap(metric = "r2.pearson", dataSets, methods, path2Result, evalSet=trainingSets[1], protein=proteinList[1], desc, testSet="actives")
###########################################################################
# DESC: analysis for the actives set
###########################################################################

analyseActivesSet <- function(path2Result, dataSet, trainSet, protein, method, desc, cutoff = 12, binsize = 0) {
  activesName = paste(path2Result, dataSet, "_", trainSet, "_", protein, "_", method, "_", concatCutoffAndBinsize(desc, cutoff, binsize), "_actives.csv", sep="")  
  activesList = read.csv(activesName, skip = 4)
  activesList = activesList[,2:3]
  
  #print(calcValidationMetric(activesList[,1], activesList[,2])$r2.pred)
  createPlot(activesList, data=paste(dataSet, trainSet), desc=desc, method=method)
}

# DESC: analysis for the actives set only
analysis_actives_set <- function() {

  for (protein in proteinList) {
    for (trainSet in trainingSets) {
      for (desc in descSets) {        
        pdf(paste(path2Result,"scatterplot_",protein,"_", trainSet,"_", desc, ".pdf",sep=""), width = 12, height = 9)
        par(mfrow = c(2,3))
        for (method in methods) {  
          for (data in dataSets) {
             analyseActivesSet(path2Result, data, trainSet, protein, method, desc)
          }        
        }
        dev.off()      
      }
    }
  }
}

#analysis_actives_set()

###########################################################################

# DESC: merge actives and inactives set to one set
createMergeList <- function(path2Result, dataSet, trainSet, protein, method, desc, cutoff = 12, binsize = 0) {
  activesName = paste(path2Result, dataSet, "_", trainSet, "_", protein, "_", method, "_", concatCutoffAndBinsize(desc, cutoff, binsize), "_actives.csv", sep="")  
  activesList = read.csv(activesName, skip = 4) 
  
  inactivesName = paste(path2Result, dataSet, "_", trainSet, "_", protein, "_", method, "_", concatCutoffAndBinsize(desc, cutoff, binsize), "_inactives.csv", sep="")  
  inactivesList = read.csv(inactivesName, skip = 4) 
  
  totalList     = rbind(activesList, inactivesList)
  totalList     = makeUniqueList(totalList, col = 2)
  return (totalList)
}

# DESC: merge the ligand ID in the result
mergeLigandID <- function(path2Result, dataSet, trainSet, protein, method, desc, cutoff = 12, binsize = 0, testSet = "") {
  if (testSet == "") {
    input = createMergeList(path2Result, dataSet, trainSet, protein, method, desc)
    activesID = read.csv(paste(path2Name,"JMJ_actives.csv",sep=""))
    activesID = makeUniqueList(activesID, col = 2)
    activesID = activesID[,c(1,2)]
    colnames(activesID)[2] = "goldscore"
    
    inactivesID = read.csv(paste(path2Name,"JMJ_inactives.csv",sep=""))
    inactivesID = makeUniqueList(inactivesID, col = 2)
    
    testSetID = rbind(activesID, inactivesID)
  }
  else {
    input = read.csv(paste(path2Result, dataSet, "_", trainSet, "_", protein, "_", method, "_", concatCutoffAndBinsize(desc, cutoff, binsize), "_", testSet, ".csv", sep=""), skip = 4)
    input = makeUniqueList(input, col = 2)    
    testSetID = read.csv(paste(path2Name,"JMJ_",testSet,".csv",sep=""))    
    testSetID = makeUniqueList(testSetID, col = 2)    
  }
  
  mergeList    = merge(input, testSetID, by.x=2 , by.y=2)
  return (mergeList[,c(5,1,3)])
}

# DESC:
calcFoundActivesPercent <- function(activesCutoff, path2Result, dataSet, trainSet, protein, method, desc, activesThreshold = 10, cutoff = 12, binsize = 0, testSet = "") {
  totalList = mergeLigandID(path2Result, dataSet, trainSet, protein, method, desc, cutoff = 12, binsize = 0, testSet)
  #print(totalList)
  sortedList = totalList[ order(-totalList[,3]),]    
  #print(sortedList)
  foundActives = sum(sortedList[1:activesCutoff, "actual"] < activesThreshold)
  #print(foundActives)
  #print(sortedList[1:activesCutoff, "actual"] < activesThreshold)
  foundActivesPercent = foundActives/activesCutoff
  
  return (foundActivesPercent)
}

# DESC:
buildAnalysis <- function(activesCutoff, descSets, activesThreshold = 10, trainSet = "refined", protein = "3PDQ_new") {
  foundActives = list()
  for (desc in descSets) {
    foundActives[[desc]] = matrix(NA, nrow = length(dataSets), ncol = length(methods))
    rownames(foundActives[[desc]]) = dataSets
    colnames(foundActives[[desc]]) = methods
    for (i in seq(length(dataSets)))
      for (j in seq(length(methods))) {
        dataSet = dataSets[i]
        method = methods[j]
        foundActives[[desc]][i,j] = calcFoundActivesPercent(activesCutoff, path2Result, dataSet, trainSet, protein, method, desc, activesThreshold=activesThreshold)
      }
    
    heatmap = foundActives[[desc]][nrow(foundActives[[desc]]):1,]
    image(1:length(methods), 1:length(dataSets), z = t(heatmap), axes = FALSE, xlab = desc, ylab = "", main = activesCutoff, cex.lab=2, cex.main = 2)
    axis(1, 1:length(methods), colnames(heatmap))
    axis(2, 1:length(dataSets), rownames(heatmap))
    for (x in 1:ncol(heatmap))
      for (y in 1:nrow(heatmap))
        text(x, y, round(heatmap[y,x], digits = 4)*100, cex = 2)
  }
}

# DESC: analyse the test set with all new features (SIFtsv2 + metallo)
analysis_test_set_all <- function() {

  for (protein in proteinList) {
      trainSet=trainingSets[1]
      pdf(paste(path2Result,"activespercent_",protein,"_", trainSet, "_", descSets[1], ".pdf",sep=""), width = 12, height = 9)
      par(mfrow = c(2,3))
      cutoffList = c(10, 20, 30, 40, 53)
      sapply(cutoffList, buildAnalysis, descSets=descSets[1], activesThreshold = 10, trainSet, protein)
      dev.off()
      for (trainSet in trainingSets) {
        pdf(paste(path2Result,"activespercent_",protein,"_", trainSet, ".pdf",sep=""), width = 12, height = 9)
        par(mfrow = c(2,5))
        cutoffList = c(10, 20, 30, 40, 53)
        sapply(cutoffList, buildAnalysis, descSets, activesThreshold = 10, trainSet, protein)
        dev.off()
      }
  }
}

#analysis_test_set_all()
###########################################################################
dataSet = dataSets[3] # 2014
trainSet = trainingSets[2] # metallo
desc = descSets[1] #SIFts
protein = proteinList[2] #T36_JMJ
method = methods[2] #RT

test = readClassifyResult(resultPath, data="CASF13_refined", evalSet = "DIG10.2", method, desc, binsize="0_XP")
#test = makeUniqueList(test, col = 1)  
test[,1] = round(test[,1],3)
activesID = read.csv(paste(path2Name,"DIG10.2_XP.csv",sep=""))
activesID = makeUniqueList(activesID, col = 2)
activesID[,2] = round(activesID[,2],3)
activesID = activesID[,c(1,2)]

mergeList    = merge(test, activesID, by.x=1 , by.y=2, sort = FALSE, all = TRUE)
sortedList = mergeList[ order(-mergeList[,2]),]
sortedList[,3] = substr(sortedList[,3],1,5)
print(unique(sortedList[,3]))


#test = mergeLigandID(path2Result, dataSet, trainSet, protein, method, desc, testSet = "test_gold")
#test = mergeLigandID(path2Result, dataSet, trainSet, protein, method, desc, testSet = "")
#sortedList = test[ order(-test[,3]),]
#colnames(sortedList)[2] = "goldscore"
#write.table(sortedList, file = "/home/dat/WORK/output/results/CASF14_metallo_T36_JMJ_Xray_RoF-RT_elementsv2-SIFt_c12b0_test_gold.csv", sep = ",", row.names = FALSE)
#for (desc in descSets) {
#  for (data in dataSets) {
#    for (method in methods) {
#      print (calcFoundActivesPercent(activesCutoff = 10, activesThreshold = 10, path2Result, data, trainSet, protein, method, desc, testSet=""))
#    }
#  }
#}
