#source("lib/configuration.R")

###########################################################################
# add ligand/protein ID to result of weka
###########################################################################
addID <- function(result, IDPath, IDFile, IDIndex = 1) {
  IDdata        = read.csv(paste(IDPath, IDFile, sep=""), na.strings=c(".", "NA", "", "?"))  
  if (IDIndex == "last") {
    IDIndex = length(colnames(IDdata))
  }
  result["ID"]  = IDdata[,IDIndex]    
  return (result)
}

###########################################################################
# preprocess the scorelist by sorting it by the absolute value
###########################################################################
preprocess_scorelist <- function(scorelist, scoreIndex = 2, IDindex, splitChar = '_', decreasing = TRUE) {
  # be sure that the scores are in absolute value
  scorelist[,scoreIndex] = abs(scorelist[,scoreIndex])
  # sorting the list
  if (decreasing) {
    scorelist = scorelist[ order(-scorelist[,scoreIndex]), ]
  }
  else {
    scorelist = scorelist[ order(scorelist[,scoreIndex]), ]
  }
  # only getting the cs-ligand name
  if (!missing(IDindex)) {
    charPos = regexpr(splitChar, scorelist[, IDindex]) - 1
    scorelist[,IDindex] = substr(scorelist[,IDindex], 1, charPos)
  }
  return (scorelist)
}

 
###########################################################################
# list top numberOfPoses (standard = 1) best poses
###########################################################################
list_top_bestpose <- function(scorelist, numberOfPoses = 1, IDindex, decreasing = TRUE) {
  poselist = preprocess_scorelist(scorelist, IDindex = IDindex, decreasing = decreasing)
  #print(length(poselist[,1]))
  #print(poselist)
  library(hash)
  keep_pose = hash(unique(poselist[,IDindex]), numberOfPoses) 
  index = 0
  for (pose in poselist[,IDindex]) {
    index = index + 1
    if (keep_pose[[pose]] > 0) {
      keep_pose[[pose]] = keep_pose[[pose]] - 1
    }
    else {
      poselist = poselist[-c(index),]     
      index = index - 1
    }
  }
  return (poselist)
} 

###########################################################################
# 
###########################################################################
calc_average_bestpose <- function(ranklist, scoreIndex = 2, IDindex) {
  library(hash)
  numberOfPoselist = hash(unique(ranklist[,IDindex]), 0) 
  totalScoreList = hash(unique(ranklist[,IDindex]), 0) 
  index = 0
  for (pose in ranklist[,IDindex]) {
    index = index + 1    
    numberOfPoselist[[pose]] = numberOfPoselist[[pose]] + 1
    totalScoreList[[pose]] = totalScoreList[[pose]] + ranklist[index, scoreIndex]
  }  
  temp = data.frame(data.frame(cbind(Name = NA, score = 1:length(keys(totalScoreList)))))
  index = 0
  for (pose in keys(totalScoreList)) {
    index = index + 1    
    totalScoreList[[pose]] = totalScoreList[[pose]] / numberOfPoselist[[pose]]
    temp[index,] = c(pose, totalScoreList[[pose]] )
  }
  # convert to float
  temp[,2] = as.numeric(temp[,2])
  # sort in the descending order
  #temp = temp[ order(-temp[,2]), ] 
  colnames(temp) = c("ID", "predicted")
  return (temp)
}