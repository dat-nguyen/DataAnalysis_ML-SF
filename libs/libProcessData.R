###############################################################
# library for processing data from MLSF Project, Feb-2016
###############################################################

# merge row data from path1 and path2, both in csv format
# if removeIndex != 0, the col with removeIndex in 2. dataset will be removed
mergeRowData <- function(path1, path2, nameIndex1, nameIndex2, removeIndex = 0, createFakeValue = FALSE) {
  data1 = read.csv(path1, na.strings=c(".", "NA", "", "?"))
  data2 = read.csv(path2, na.strings=c(".", "NA", "", "?"))
  if (missing(nameIndex1)) # if the index is missing, then it's the last col of the set
    nameIndex1 = length(data1[1,])
  if (missing(nameIndex2)) # if the index is missing, then it's the last col of the set
    nameIndex2 = length(data2[1,])
  if (removeIndex > 0) { # sometimes we need to remove a col @ removeIndex @ the second dataset
    data2 = data2[,-removeIndex]
    if (nameIndex2 > removeIndex) # if the col to be merged from second set is at right position from the removed col, need to adjust its position as well
      nameIndex2 = nameIndex2-1
  }
  else if (removeIndex == -1) { # remove the last col of second dataset
    data2 = data2[,-(length(data2[1,]))]    
    nameIndex2 = nameIndex2-1
  }
  mergeData = merge(data1, data2, by.x = nameIndex1, by.y = nameIndex2)
  # create a fake col with fake measurement data for weka, all the fake data are unique
  if (createFakeValue) {
    firstCol  = mergeData[,1]
    firstColName = colnames(mergeData)[1]
    fakeCol     = seq(from = 0, by=0.001, length.out = nrow(mergeData))
    mergeData = cbind(firstCol, "RMSD" = fakeCol, mergeData[,-1])
    colnames(mergeData)[1] = firstColName
  }
  
  return (mergeData)
}

# merge col data from path1 and path2, both in csv format
# 
mergeColData <- function(path1, path2) {
  data1 = read.csv(path1, na.strings=c(".", "NA", "", "?"))
  data2 = read.csv(path2, na.strings=c(".", "NA", "", "?"))  
  mergeData = rbind(data1, data2)
  
  return (mergeData)
}

# merge description data from 2 sets to 1 set (in other words, create combi desc data)
mergeDesc <- function(dataPath, CASFset, trainingSet, createFakeValue = FALSE) {
  print(CASFset)
  path1 = paste(dataPath, CASFset, "_RMSD_", trainingSet, "-", DESCRIPTORS[1], ".csv", sep="")
  path2 = paste(dataPath, CASFset, "_RMSD_", trainingSet, "-", DESCRIPTORS[2], ".csv", sep="")
  print(path1)
  print(path2)
  # all training data always has the first col as the pKd/RMSD value, so after merging it needs to be removed to avoid redundant, that's why removeIndex = 1
  mergeData12 = mergeRowData(path1, path2, removeIndex = 1, createFakeValue=createFakeValue)
  # write to file
  path12 = paste(dataPath, "processed/", CASFset, "_RMSD_", trainingSet, "-", DESCRIPTORS[1], "-", DESCRIPTORS[2], ".csv", sep="")
  print(path12)
  write.table(mergeData12, file = path12, sep = ",", quote=FALSE, row.names = FALSE)
  
  # for merging the third desc set into one, but for later
  # \TODO: move into a new function
  #  path3 = paste(dataPath, CASFset, "_xscore.csv", sep="")
  #  mergeData123 = mergeDesc(path12, path3, nameIndex1=1,nameIndex2=1)
  #  path123 = paste(dataPath, CASFset,"_elementsv2-SIFt-xscore.csv", sep="")
  #  write.table(mergeData123, file = path123, sep = ",", quote=FALSE, row.names = FALSE)  
}

