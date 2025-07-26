#install and load dependencies
source("qPCRfunctions.R")
source("dependencies.R")

# enter unique experiment identifier
expID <- "dsCHS1"
# enter your reference gene
RefG <- "RpS18" 
#experimental group
expGroup <- "dsCHS1"
#Reference group
refGroup <- "dsVer"
# plot data sorted by primer pairs
level_list <- c("CHS1")#c("Sec23", "E75", "EcRA", "G8A", "CHS1")
# enter path to prepared data in *.xlsx file
loadPath <-  "C:/Users/beckm/OneDrive/Desktop/Kim_RNAi/2507dsCHS1.xlsx"
# enter path where data shall be stored
storePath <- "C:/Users/beckm/OneDrive/Desktop/Kim_RNAi"

# if you want to show stats, set TRUE
add_statistics <-F

# load data from different sheets if necassary
DataTable <- load_Data(loadPath)


# prepare your data and calc Cts
dataPool <- rel_data_prep(DataTable, RefG, genes, refGroup)

dataPool2 <- dataPool 

rel_plot(dataPool,bars = "se") #optional add: bars = "se/sd/ci" se is standard

ggsave(paste0(expGroup, "_rel_plot.png"), path = storePath)
ggsave(paste0(expGroup, "_rel_plot.pdf"), path = storePath)

# save data
saveData(dataPool, extractedData, storePath, expID)

# clean up only if you're certain to remove all data from the workspace
# saved data will NOT be deleted from the disk
#rm(list = ls())








