#lig <- unique(read.table("/home/dat/Downloads/ligands.txt")[,1]);
#uniqRes <- read.table("/home/dat/Downloads/dataforR_uq.txt",header=T);
#colnames(uniqRes)[1]="LigandName";
#uniqRes$IsActive <- as.numeric(uniqRes$LigandName %in% lig)
#predINTERuq <- prediction(uniqRes$INTER*-1, uniqRes$IsActive)

#data1 = read.csv("/home/dat/WORK/DB/DESCRIPTORS/DIG10.2/confgen/DIG10.2_elementsv2_c12b0_SP.csv", na.strings=c(".", "NA", "", "?"))
#data2 = read.csv("/home/dat/WORK/DB/DESCRIPTORS/DIG10.2/confgen/DIG10.2_SIFt_c12b0_SP.csv", na.strings=c(".", "NA", "", "?"))
#data1 = read.csv("/home/dat/WORK/DB/DESCRIPTORS/CASF/CASF14_refined_SIFt_c12b0.csv", na.strings=c(".", "NA", "", "?"))
#data2 = read.csv("/home/dat/WORK/DB/DESCRIPTORS/CASF/CASF14_refined_elementsv2_c12b0.csv", na.strings=c(".", "NA", "", "?"))
#nameIndex1 = length(data1[1,])
#nameIndex2 = length(data2[1,])
#score = read.csv("/home/dat/WORK/DB/DESCRIPTORS/DIG10.2/cs-confgen_all_SP.csv", na.strings=c(".", "NA", "", "?"))
#score = read.csv("/home/dat/WORK/DB/DESCRIPTORS/CASF/CASF14_refined_xscore.csv", na.strings=c(".", "NA", "", "?"))
#mergeData = merge(mergeData12, score, by.x = 1, by.y = 1)

path = "/home/dat/MyClouds/reports/pr_2015-06-23/"
#file = "Spearman_DIG10.2_elementsv2-SIFt_c12b0-xscore_Pose_1.csv"
file = "Spearman_DIG10.2_elementsv2-SIFt_c12b0-xscore_Pose_1.csv"
filePrefix = "Spearman_DIG10.2"
#filePrefix = "Spearman_reduced_RMSD_DIG10.2"
filePrefix = "Spearman_refined_RMSD_DIG10.2"

createHeatmap <- function(desc, numberOfPoses) {
  result = read.csv(paste(path, filePrefix, desc, '_Pose_', numberOfPoses, ".csv", sep=""))
  
  heatmap = result
  
  methods = gsub(desc, "", rownames(heatmap))
  methods = gsub("_RoF", "", methods)
  methods = methods[c(length(methods):1)]
  targets = colnames(heatmap)
  heatmap = heatmap[nrow(heatmap):1,]
  
  textDesc = ""
  if (desc == "_elementsv2-SIFt_c12b0-xscore") {
    textDesc = "ESX"
  } else if (desc == "_elementsv2-SIFt_c12b0") {
    textDesc = "ES" }
  
  mainLabel = paste("Spearman - RMSD - descriptors", textDesc, "pose", numberOfPoses)
  image(1:length(targets), 1:length(methods), z = t(heatmap), axes = FALSE, xlab = "Poses gen. by docking method", ylab = "", main = mainLabel, cex.lab=2, cex.main = 2)
  axis(1, 1:length(targets), targets)
  axis(2, 1:length(methods), methods)
  for (x in 1:ncol(heatmap))
    for (y in 1:nrow(heatmap))
      text(x, y, round(heatmap[y,x], digits = 3), cex = 1.8)
  
}

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