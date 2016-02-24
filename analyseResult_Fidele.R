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
#  hist(ligandList[,2], xlim=c(-15,10), col=rgb(0,0,1,1/4), border=F)
#  hist(ligandList[,3], col=rgb(1,0,0,1/4), add=T)

readPredictScore <- function(target, classifyMethod="CASFv2014-refined_sampling_clusters10_RF_") {  
  ligandList = read.csv(paste0(PREDICTED_SCORES_PATH, target, ".txt"), header=FALSE)
  # modify the ligand name because of splitting problem
  ligandList[,1] = gsub(' $', '', ligandList[,1])
  ligandList[,1] = gsub(' ', splitChar, ligandList[,1])
  
  scoreFile   = paste0(SPLIT_RESULT_PATH, "Fidele/", classifyMethod, "_", target,".csv")
  #print(classifyMethod)
  scoreList   = read.csv(scoreFile, na.strings=c(".", "NA", "", "?"))
  # modify the ligand name because of splitting problem
  scoreList[,1] = gsub(pattern=paste0(splitChar,"[0-9]+$"), replacement="", scoreList[,1])
  
  ligandList    = merge(ligandList, scoreList, by.x = 1, by.y = 1)
  colnames(ligandList)[1] = "ligand"
  
  return (ligandList)
}

plotHistScore <- function(index) {
  africa = readPredictScore(FideleList[index])
  npact = readPredictScore(FideleList[index+2])
  mycolor = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4))
  hist(africa[,2], breaks= 30, ylim=c(0, 2000), col=mycolor[1], border=F, main=substr(FideleList[index],1,4), xlab = "XP score")
  hist(npact[,2], breaks= 30, col=mycolor[2], add=T)
  legend('topleft',c('Africa','NPACT'),
         fill = c(mycolor[1], mycolor[2]), bty = 'n',
         border = NA)
}
#pdf(paste0(RESULT_PATH, "Fidele_3E37_5P21_XP.pdf"), width = 12, height = 9)
#par(mfrow = c(1,2))
#plotHist(index = 2)
#plotHist(index = 22 )
#dev.off()

plotHistRMSD <- function(index, classifyMethod) {
  africa = readPredictScore(FideleList[index], classifyMethod)
  npact = readPredictScore(FideleList[index+2], classifyMethod)
  mycolor = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4))
  hist(africa[,3], breaks= 30, ylim=c(0, 2000), col=mycolor[1], border=F, main=substr(FideleList[index],1,4), xlab = "RMSD")
  hist(npact[,3], breaks= 30, col=mycolor[2], add=T)
  legend('topleft',c('Africa','NPACT'),
         fill = c(mycolor[1], mycolor[2]), bty = 'n',
         border = NA)
}
for (CASFset in CASF_SETS) {
  pdf(paste0(RESULT_PATH, "Fidele_3E37_5P21_RMSD",CASFset, ".pdf"), width = 12, height = 9)
  par(mfrow = c(2,6))
  for (trainSet in TRAINING_SETS)
     for (method in ML_METHODS) {
        #par(mfrow = c(1,2), new=TRUE)
        classifyMethod = paste(CASFset, trainSet, method, sep="_")
        plotHistRMSD(index = 2, classifyMethod )
        plotHistRMSD(index = 22, classifyMethod )
     }
  dev.off()
}


###########################################################################
#analysis_test_set_all()
#analysis_test_set_all(typeOfScore="NormA")
#analysis_test_set_all(typeOfScore="NormB")