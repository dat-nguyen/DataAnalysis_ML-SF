# configuration file
# contains all constants

#resultPath  = "/home/dat/WORK/output/results/2015-01-13/PDBbind/"
resultPath  = "/home/dat/WORK/output/results/2015-06-23/"
path2Result = "/home/dat/WORK/output/results/2015-06-23/"
#path2Result = "results/2015-01-13/"
path2Name   = "/home/dat/WORK/DB/DESCRIPTORS/CASF/"
#path2Name   = "New/"

dataSets = c("CASF12", "CASF13", "CASF14")
descSets = c("elementsv2-SIFt_c12b0", "elementsv2-SIFt_c12b0-xscore")
methods  = c("RoF-REP", "RoF-RT")
evalSets = c("training", "test")
trainingSets = c("refined", "metallo")

proteinList         = c("3PDQ_new", "T36_JMJ_Xray")
stateList           = c("actives", "inactives")

TRUE_RANK_CSAR13 = c(7, 8, 6, 4, 3, 5, 2, 1, 10, 9)
dockingMethods = c("asp", "plp", "chemscore", "goldscore", "SP", "XP")