
#proteinList = read.csv("/Users/knight/Dropbox/data/CASF/CASF14_refined_SIFt_c12b0.csv", na.strings=c(".", "NA", "", "?"))
#print(proteinList[,length(colnames(proteinList))])

elements <- function() {
  misc = read.csv("/Users/knight/Dropbox/data/RMSD/CASF14_refined_RMSD_elementsv2_misc.csv", na.strings=c(".", "NA", "", "?"))
  gold = read.csv("/Users/knight/Dropbox/data/RMSD/CASF14_refined_RMSD_elementsv2_gold.csv", na.strings=c(".", "NA", "", "?"))
  
  total = rbind(misc, gold)
  #lengthCol = length(colnames(total))
  total = total[order(total$ligandID),]
  write.table(total, file = "/Users/knight/Dropbox/data/RMSD/CASF14_refined_RMSD_elementsv2.csv", sep = ",", row.names = FALSE)  
}

credo <- function() {
  misc = read.csv("/Users/knight/Dropbox/data/RMSD/CASF14_refined_RMSD_SIFt_c12b0_misc.csv", na.strings=c(".", "NA", "", "?"))
  gold = read.csv("/Users/knight/Dropbox/data/RMSD/CASF14_refined_RMSD_SIFt_c12b0_gold.csv", na.strings=c(".", "NA", "", "?"))
  
  total = rbind(misc, gold)
  total = total[order(total$ligandID),]
  write.table(total, file = "/Users/knight/Dropbox/data/RMSD/CASF14_refined_RMSD_SIFt_c12b0.csv", sep = ",", row.names = FALSE)  
}

xscore <- function() {
  asp = read.csv("/Users/knight/Dropbox/data/RMSD/CASF14_refined_RMSD_xscore_asp.csv", na.strings=c(".", "NA", "", "?"))
  plp = read.csv("/Users/knight/Dropbox/data/RMSD/CASF14_refined_RMSD_xscore_plp.csv", na.strings=c(".", "NA", "", "?"))
  gold = read.csv("/Users/knight/Dropbox/data/RMSD/CASF14_refined_RMSD_xscore_goldscore.csv", na.strings=c(".", "NA", "", "?"))
  chem = read.csv("/Users/knight/Dropbox/data/RMSD/CASF14_refined_RMSD_xscore_chemscore.csv", na.strings=c(".", "NA", "", "?"))
  
  totalxscore = rbind(asp, plp, gold, chem)
  totalxscore = totalxscore[order(totalxscore$ID),]
  totalxscore[,1] = gsub(".mol2", "", totalxscore[,1])
  
  write.table(totalxscore, file = "/Users/knight/Dropbox/data/RMSD/CASF14_refined_RMSD_xscore.csv", sep = ",", row.names = FALSE)  
}


#all = merge(total, totalxscore, by.x = length(colnames(total)), by.y = 1)
#write.table(all, file = "/Users/knight/Dropbox/data/RMSD/CASF14_refined_RMSD_elementsv2_c12b0-xscore_process.csv", sep = ",", row.names = FALSE)  

#desc = "elementsv2-SIFt_c12b0"
desc = "elementsv2-SIFt_c12b0-xscore"

#splitData <- function(desc) {
  total = read.csv(paste("/Users/knight/Dropbox/data/RMSD/CASF14_refined_RMSD_", desc, "_process.csv", sep=""), na.strings=c(".", "NA", "", "?"))
  nTrainData 	= nrow(total)
  #siz		= round((10*nTrainData)/100)
  #newDataSize = round(nTrainData/8)
  newDataSize = 86150
  # set seed to "randomness"
  as.numeric(Sys.time())-> t; set.seed((t - floor(t)) * 1e8 -> seed);
  choosenPDB	= sample(nTrainData, size=newDataSize)
  totalNew	= total[ choosenPDB,]  
  
  write.table(totalNew, file = paste("/Users/knight/Dropbox/data/RMSD/CASF14_reduced_RMSD_", desc, ".csv", sep=""), sep = ",", row.names = FALSE)  
#}

