#install and load dependencies
source("qPCRfunctions.R")
source("dependencies.R")

# enter unique experiment identifier
expID <- "Coreg"
# enter your reference gene
RefG <- "RpS18" 
#experimental group
expGroup <- "dsKu70"
#Reference group
refGroup <- "dsVer"
# enter path to prepared data in *.xlsx file
loadPath <- "C:/Users/Marius Beck/OneDrive/Desktop/RNAi/LSRNAi/Sammelmappe_extra_Rebuf.xlsx"
# enter path where data shall be stored
storePath <- "C:/Users/Marius Beck/OneDrive/Desktop/RNAi/LSRNAi/"

# if you want to show stats, set TRUE
add_statistics <-F

# load data from different sheets if necassary
DataTable <- load_Data(loadPath)

# prepare your data and calc Cts
dataPool <- rel_data_prep(DataTable, RefG, genes, refGroup)


# plot data sorted by primer pairs
level_list <- c("FABP15275", # 1st gene should always be the KnockDown!!!
                "FABP1310",
                "FABP12473",
                "Reb",
                "CHS1",
                "Ctl2", 
                "Tsr", 
                "Osi18",
                "SERCA",
                "Lac2",
                "AGM",
                "EcRA"
                )

rel_plot(dataPool) #optional add: bars = "se/sd/ci" se is standard

ggsave(paste0(expGroup, "_coreg_plot.png"), path = storePath)
ggsave(paste0(expGroup, "_coreg_plot.pdf"), path = storePath)

# save data
saveData(dataPool, extractedData, storePath, expID)

# clean up only if you're certain to remove all data from the workspace
# saved data will NOT be deleted from the disk
rm(list = ls())








