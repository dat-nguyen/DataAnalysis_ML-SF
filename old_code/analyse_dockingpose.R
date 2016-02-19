#SP = read.csv("/home/dat/WORK/DB/DESCRIPTORS/DIG10.2/confgen/cs-confgen_all_SP.csv")
#XP = read.csv("/home/dat/WORK/DB/DESCRIPTORS/DIG10.2/confgen/cs-confgen_all_XP.csv")
ASP = read.csv("/home/dat/WORK/DB/DESCRIPTORS/DIG10.2/confgen/cs-confgen_all_asp.csv")
#PLP = read.csv("/home/dat/WORK/DB/DESCRIPTORS/DIG10.2/confgen/cs-confgen_all_plp.csv")
#chem = read.csv("/home/dat/WORK/DB/DESCRIPTORS/DIG10.2/confgen/cs-confgen_all_chemscore.csv")
gold = read.csv("/home/dat/WORK/DB/DESCRIPTORS/DIG10.2/confgen/cs-confgen_all_goldscore.csv")

true_rank = c(7, 8, 6, 4, 3, 5, 2, 1, 10, 9)

preprocess_scorelist <- function(scorelist) {
  # be sure that the scores are in absolute value
  scorelist[,2] = abs(scorelist[,2])
  # sorting the list
  scorelist = scorelist[ order(-scorelist[,2]), ]
  # 
  scorelist[,1] = substr(scorelist[,1],1,5)
  
  return (scorelist)
}

rank_bestpose <- function(scorelist) {
  scorelist = preprocess_scorelist(scorelist)
  return (unique(scorelist[,1]))
}

rank_average_bestpose <- function(scorelist, numberOfPoses = 3) {
  poselist = preprocess_scorelist(scorelist)
  library(hash)
  keep_pose = hash(unique(poselist[,1]), numberOfPoses) 
  index = 0
  for (pose in poselist[,1]) {
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

calc_average_bestpose <- function(ranklist) {
  library(hash)
  numberOfPoselist = hash(unique(ranklist[,1]), 0) 
  totalScoreList = hash(unique(ranklist[,1]), 0) 
  index = 0
  for (pose in ranklist[,1]) {
    index = index + 1    
    numberOfPoselist[[pose]] = numberOfPoselist[[pose]] + 1
    totalScoreList[[pose]] = totalScoreList[[pose]] + ranklist[index, 2]
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
  return (temp)
}
#print(rank_bestpose(SP))
#print(rank_bestpose(XP))
print(rank_bestpose(ASP))
#print(rank_bestpose(PLP))
#print(rank_bestpose(chem))
print(rank_bestpose(gold))
#print(rank_average_bestpose(PLP))
test = (calc_average_bestpose(rank_average_bestpose(gold)))
print(test)
print(rank(-test[,2]))
print(cor(true_rank, rank(-test[,2]), method="spearman"))

# test
x = c(9,7,10,4,3,6,2,1,8,5)
print(cor(true_rank, x, method="spearman"))