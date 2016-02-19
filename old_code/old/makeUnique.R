#source("lib/lib.data.R")

# 
concatName <- function(path, desc, method, trainset, state) {
  if (desc == "credo") 
    return (paste(path, "JMJD2A_", trainset, "_", desc, "_", state, "_c12_RoF_RC_", method,".csv", sep=""))
  else
    return (paste(path, "JMJD2A_", trainset, "_", desc, "_", state, "_c12b12_RoF_RC_", method,".csv", sep=""))
}

path2Result = "/Users/knight/dev/classification_R/JMJD/Result/"
path2Name   = "/Users/knight/dev/classification_R/JMJD/"
#path2Result = "/home/dat/WORK/DB/JMJD/Result/"
#path2Name   = "/home/dat/WORK/DB/DESCRIPTORS/JMJD/"

descList       = c("credo", "elements", "sybyl")
methodList     = c("REP", "RT")
trainsetList   = c("CASF12", "CASF13", "CSAR12")
stateList      = c("actives", "inactives")

# 
createMergeList <- function(desc, method, trainset, state) {
  if (state == "actives") {
    name     = read.csv(paste(path2Name,"JMJD2A_dock_actives.csv",sep=""), sep = "\t")
  }
  else {
    name     = read.csv(paste(path2Name,"JMJD2A_dock_inactives_gold.csv",sep=""), sep = "\t")
  }
  result   = read.csv(concatName(path2Result, desc, method, trainset, state)) 
  
  name     = makeUniqueList(name)
  result   = makeUniqueList(result[,2:3])
  mergeList    = merge(name, result, by.x=1 , by.y=1)
  if (state == "actives") {  
    mergeList    = data.frame(goldscore = rep(0, length(mergeList[,1])), mergeList)
  }
  else {
    mergeList    = data.frame(pIC50 = rep(0, length(mergeList[,1])), mergeList)
  }
  return (mergeList)
}

createSortedList <- function(desc, method, trainset) {
  #activesList = createMergeList(descList[1], methodList[1], trainsetList[1], stateList[1])
  activesList = createMergeList(desc, method, trainset, stateList[1])
  inactivesList = createMergeList(desc, method, trainset, stateList[2])
  totalList = rbind(activesList, inactivesList)
  sortedList = totalList[ order(-totalList[,4]),]
  
  return (sortedList)
}

calcFoundActivesPercent <- function(desc, method, trainset, cutoff) {
  sortedList = createSortedList(desc, method, trainset)
  foundActives = sum(sortedList[1:cutoff, "pIC50"] > 0)
  foundActivesPercent = foundActives/cutoff
  
  return (foundActivesPercent)
}

buildAnalysis <- function(cutoff) {
  foundActives = list()
  for (method in methodList) {
    foundActives[[method]] = matrix(NA, nrow = length(trainsetList), ncol = length(descList))
    rownames(foundActives[[method]]) = trainsetList
    colnames(foundActives[[method]]) = descList
    # for (trainset in trainsetList)
    for (i in seq(length(trainsetList)))
      # for (desc in descList)
      for (j in seq(length(descList))) {
        trainset = trainsetList[i]
        desc = descList[j]
        foundActives[[method]][i,j] = calcFoundActivesPercent(desc, method, trainset, cutoff)
      }
    # \TODO: fixing bugs, trainsetList verwechselt mit descList?
    image(1:length(trainsetList), 1:length(descList), z = t(foundActives[[method]]), axes = FALSE, xlab = method, ylab = "", main = cutoff)
    axis(1, 1:length(trainsetList), colnames(foundActives[[method]]))
    axis(2, 1:length(descList), rownames(foundActives[[method]]))
  }
}

path = "/Users/knight/"
pdf(paste(path,"JMJD2A_heatmap.pdf",sep=""), width = 8, height = 12)
par(mfrow = c(4,2))
cutoffList = c(10, 20, 30, 46)
sapply(cutoffList, buildAnalysis)
dev.off()

source("lib/lib.validation.extern.R")

plotActives <- function() {
  for (trainset in trainsetList)
    for (desc in descList) 
      for (method in methodList) {
        result   = read.csv(concatName(path2Result, desc, method, trainset, "actives")) 
        plot(result[,2], result[,3], xlab = "measured pIC50", ylab = "predicted pKd", main = paste(method,"_",trainset,"_",desc,sep=""))
        legend("topleft", paste("R=",round(sqrt(calcValidationMetric(result[,2], result[,3])$r2.pearson), 3), sep=""))
      }
}

pdf(paste(path,"JMJD2A_plot.pdf",sep=""), width = 16, height = 8)
par(mfrow = c(3,6))
plotActives()
dev.off()


#cutoff = 40
#activesList = createMergeList(desc, method, trainset, stateList[1])
#inactivesList = createMergeList(desc, method, trainset, stateList[2])
#totalList = rbind(activesList, inactivesList)
#sortedList = totalList[ order(-totalList[,4]),]
#foundActives = sum(sortedList[1:cutoff, "pIC50"] > 0)
#foundActivesPercent = foundActives/cutoff