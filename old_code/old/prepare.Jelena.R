descList       = c("credo", "elements", "sybyl")
binsizeList    = c(0, 12)
cutoff         = 12
trainingList   = c("CASF12", "CASF13")
# 
addID2DescData <- function(dataPath = "/home/dat/WORK/DB/DESCRIPTORS/Jelena/", DBsetPrefix = "140722_set1_moe_") {
    for (desc in descList) 
      for (binsize in binsizeList) {    
        dataName = paste(DBsetPrefix, concatCutoffAndBinsize(desc, cutoff, binsize), sep="")
        data = read.csv(paste(dataPath, dataName, ".csv", sep=""))
        # create new table with ligand name and numeric ID
        ligandID = t(rbind(levels(data[,ncol(data)]), c(data[,ncol(data)])))
        colnames(ligandID) = c("ligand", "ID")
        # create new data with extra numeric ID
        newData = merge(ligandID, data, by.y = ncol(data), by.x = 1)
    
        # write to file
        dataNameWithID = paste(DBsetPrefix, concatCutoffAndBinsize(desc, cutoff, binsize),"_ID", sep="")
        ligandNameID   = paste(DBsetPrefix, "ID", sep="")

        write.table(newData[,2:ncol(newData)], file = paste(dataPath, "WithID/", dataNameWithID, ".csv", sep=""), sep = ",", row.names = FALSE)
      }
  write.table(ligandID, file = paste(dataPath, "WithID/", ligandNameID, ".csv", sep=""), sep = ",", row.names = FALSE)
}

combineResult <- function(resultPath = "/home/dat/WORK/dev/weka/results/", DBsetPrefix = "140722_set1_moe_") {
  methodList = c("RT", "REP")
  # read the ID file first
  dataWithFullName = read.csv(paste("/home/dat/WORK/DB/DESCRIPTORS/Jelena/WithID/",DBsetPrefix,"ID.csv", sep=""))
  # exchange the cols
  dataWithFullName = dataWithFullName[,2:1]
  for (trainingSet in trainingList)           
    for (desc in descList) 
      for (binsize in binsizeList)  
        for (method in methodList) {
          resultName = paste(trainingSet, "-", method, "-", DBsetPrefix, concatCutoffAndBinsize(desc, cutoff, binsize), sep="")
          # read the classify data, skipping the first 4 lines
          data = read.csv(paste(resultPath, resultName, ".csv", sep=""), skip=4)
          # merge the full name data with classify result
          dataWithFullName = merge(dataWithFullName, data[,2:3], by.x = 1, by.y = 1)
          len = ncol(dataWithFullName)
          colnames(dataWithFullName)[len] = paste(trainingSet,"_",method,"_",concatCutoffAndBinsize(desc, cutoff, binsize),sep="")
          totalResultsName = paste(resultPath,"Result_",DBsetPrefix, ".csv", sep="")
          write.table(dataWithFullName, file = totalResultsName, sep = ",", row.names = FALSE)        
        }  
  totalResultsName = paste(resultPath,"Result_",DBsetPrefix, ".csv", sep="")
  write.table(dataWithFullName, file = totalResultsName, sep = ",", row.names = FALSE)
}

resultPath = "/home/dat/WORK/dev/weka/results/"
#DBsetPrefix = "140722_set1_moe_"
#method = "REP"
#desc = "credo"
#binsize = 12
JelenaAllSets = c("140722_set1_moe_", "140722_set2_moe_", "140722_set3_moe_", "140722_set2_maestro_", "140722_set3_maestro_")
for (eachset in JelenaAllSets) {
#  combineResult(DBsetPrefix = eachset)
}
combineResult(DBsetPrefix = "140722_set1_maestro_")
#addID2DescData(DBsetPrefix = "140722_set1_moe_")
#addID2DescData(DBsetPrefix = "140722_set1_maestro_")
#addID2DescData(DBsetPrefix = "140722_set2_moe_")
#addID2DescData(DBsetPrefix = "140722_set2_maestro_")
#addID2DescData(DBsetPrefix = "140722_set3_moe_")
#addID2DescData(DBsetPrefix = "140722_set3_maestro_")
