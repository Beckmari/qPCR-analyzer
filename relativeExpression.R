

#install and load dependencies
source("qPCRfunctions.R")
source("dependencies.R")

# enter unique experiment identifier
expID <- "dsFABP3_Lac2_pool"
# enter your reference gene
RefG <- "RpS18" 
# enter all your genes w/o reference gene
genes <- list("Lac2") 
# enter path to prepared data in *.xlsx file
loadPath <- "C:/Users/Marius Beck/OneDrive/Desktop/FABP/dsFAPB/dsFABP-KDs/Coreg/Lac2/dsFABP_Lac2_coreg.xlsx"

# enter path where data shall be stored
storePath <- "C:/Users/Marius Beck/OneDrive/Desktop/FABP/dsFAPB/dsFABP-KDs/Coreg/Lac2/"
# if you want to show stats, set TRUE
add_statistics <-T

# enter your empirical primerefficancies in the following list separated by comma
# if you don't know your empirical primerefficancies enter a theoretical value of 2
efficancies <- list(
  RpS6 = 1.99, #empiric
  RpS18 = 1.93, #empiric
  FABP15275 = 1.91, #empiric
  FABP12473 = 2.00, #theoretic
  FABP1310 = 2.00, #theoretic
  CHS1 = 2.00, #empiric
  Exp = 1.95, #empiric
  Osiris18 = 2.00, #theoretic
  Tsr = 2.00, #theoretic
  Ctl2_1323 = 2.00, #theoretic
  SERCA3 = 2.049, #empiric
  SERCA2 = 2.0, #theoretic
  Tc6639 = 2.00, #theoretic
  G8A = 2.04, #theoretic
  E75 = 2.02, #theoretic
  KU70 = 2.00, #theoretic
  KU80 = 2.00, #theoretic
  Lac2 = 2.00, #theoretic
  EcRA = 1.95 #theoretic
)

# load data

DataTable1 <- read_excel(loadPath, sheet = 1)


DataTable1$Cq <- gsub(",", ".", DataTable1$Cq)
DataTable1$Cq <- gsub("-1", NA, DataTable1$Cq)
DataTable1$Cq <- as.numeric(DataTable1$Cq)
DataTable1 <- na.omit(DataTable1)

DataTable2 <- read_excel(loadPath, sheet = 2)



DataTable2$Cq <- gsub(",", ".", DataTable2$Cq)
DataTable2$Cq <- gsub("-1", NA, DataTable2$Cq)
DataTable2$Cq <- as.numeric(DataTable2$Cq)
DataTable2 <- na.omit(DataTable2)








DataTable1 <- aggregate(DataTable1, list(DataTable1$Sample, DataTable1$Gene), mean)
DataTable2 <- aggregate(DataTable2, list(DataTable2$Sample, DataTable2$Gene), mean)

DataTable1$Sample <- DataTable1$Group.1
DataTable1$Gene <- DataTable1$Group.2
DataTable2$Sample <- DataTable2$Group.1
DataTable2$Gene <- DataTable2$Group.2
DataTable <- rbind(DataTable1, DataTable2)


# prepare your data
ValueTable <- dataPrep(DataTable2, RefG, genes)

# normalize Ct values with Primerefficancies against reference gene and calculate deltaCt
uni_samples <- unique(ValueTable$Sample)

dataPool <- c()
tempData <- c()
for (sample_name in uni_samples) {
  var_name <- paste0(sample_name, "_plot")
  for (gene_name in genes) {
    tempData <- plotPrep(efficancies, RefG, gene_name, ValueTable)# get(sample_name))
    dataPool <- rbind(dataPool, tempData)
  }
}
dataPool <- as.data.frame(dataPool)
# check for n < 2 in dataPool
duplicates <- duplicated(dataPool$Sample) | duplicated(dataPool$Sample, fromLast = TRUE)
extractedData <- dataPool[!duplicates, ]
dataPool <- dataPool[duplicates, ]

#dataPool$Sample <- factor(dataPool$Sample, levels = c("L1", "L2", "L3", "L4", "L5", "PP", "P0", "P1", "P2", "P3", "P4", "P5", "A"))


temp <- subset(dataPool1, grepl("dsVer", dataPool1$Sample))
temp <- 1 / median(temp$normCt)

dataPool1$relExp <- dataPool1$normCt * temp

temp <- subset(dataPool2, grepl("dsVer", dataPool2$Sample))
temp <- 1 / median(temp$normCt)

dataPool2$relExp <- dataPool2$normCt * temp
#
dataPool <- rbind(dataPool1, dataPool2)








# plot data sorted by primer pairs

dataPool$Sample <- factor(dataPool$Sample, levels = c("dsVer", "dsFABP"))
for (gene_name in genes) {
  tempPlotData <- subset(dataPool, grepl(gene_name, dataPool$Gene))
  plotData(tempPlotData, gene_name, add_statistics)
  # save plot and data
  ggsave(paste0(gene_name, "_normalized_boxplot1.tiff"), path = storePath)
  ggsave(paste0(gene_name, "_normalized_boxplot1.pdf"), path = storePath)
}


# plot data sorted by RNAi
rnai <- unique(dataPool$Sample)
rnai <- subset(rnai, !grepl("dsVer", rnai))
ctrlKD <- subset(dataPool, grepl("dsVer", dataPool$Gene))
for (kd in rnai) {
  tempPlotData <- subset(dataPool, grepl(kd, dataPool$Sample))
  #tempPlotData <- rbind(tempPlotData, ctrlKD)
  plotData(tempPlotData, kd, add_statistics)
  # save plot and data
  ggsave(paste0(kd, "_normalized_boxplot1.tiff"), path = storePath)
  ggsave(paste0(kd, "_normalized_boxplot1.pdf"), path = storePath)
}

# save data
saveData(dataPool, extractedData, storePath, expID)



# clean up only if you're certain to remove all data from the workspace
# saved data will NOT be deleted from the disk
rm(list = ls())

