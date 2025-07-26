#install and load dependencies
source("qPCRfunctions.R")
source("dependencies.R")

# enter unique experiment identifier
expID <- "dsFABP8"
# enter your reference gene
RefG <- "RpS18" 
# enter all your genes w/o reference gene
genes <- list( "Dop2R", "Dop1R2", "TH", "CHS1", "Lac2") 
# enter path to prepared data in *.xlsx file
loadPath <-  "C:/Users/beckm/Desktop/dsFABP8_coreg.xlsx"
# enter path where data shall be stored
storePath <- "C:/Users/beckm/Desktop/"
# if you want to show stats, set TRUE
add_statistics <-F

# enter your empirical primerefficancies in the following list separated by comma
# if you don't know your empirical primerefficancies enter a theoretical value of 2
efficancies <- list(
  RpS6 = 1.99, #empiric
  RpS18 = 1.93, #empiric
  FABP8 = 1.91, #empiric 15275
  FABP2 = 2.00, #theoretic 12473
  FABP4 = 2.00, #theoretic 1310
  CHS1 = 2.00, #empiric
  Reb = 1.95, #empiric
  Osi18 = 2.00, #theoretic
  Tsr = 2.00, #theoretic
  Ctl2_1323 = 2.00, #theoretic
  SERCA12671 = 2.049, #empiric
  SERCA2 = 2.0, #theoretic
  Tc6639 = 2.00, #theoretic
  G8A = 2.04, #theoretic
  E75 = 2.02, #theoretic
  KU70 = 2.00, #theoretic
  KU80 = 2.00, #theoretic
  Lac2 = 2.00, #theoretic
  TH = 2.00, #theoretic
  Dop1R2 = 2.00, #theoretic
  Dop2R = 2.00, #theoretic
  Lac2 = 2.00, #theoretic
  EcRA = 1.95 #theoretic
)

# load data

DataTable <- read_excel(loadPath, sheet = 1)

genes <- unique(DataTable$Gene)
genes <- subset(genes, !grepl(RefG, genes))

DataTable$Cq <- gsub(",", ".", DataTable$Cq)
DataTable$Cq <- gsub("-1", NA, DataTable$Cq)
DataTable$Cq <- as.numeric(DataTable$Cq)
DataTable <- na.omit(DataTable)

# prepare your data
ValueTable <- dataPrep(DataTable, RefG, genes)

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
dataPool <- distinct(dataPool)
#dataPool$Sample <- factor(dataPool$Sample, levels = c("L1", "L2", "L3", "L4", "L5", "PP", "P0", "P1", "P2", "P3", "P4", "P5", "A"))
dataPool$Sample <- factor(dataPool$Sample, levels = c("Elytra P0", "Elytra P1", "Elytra P2", "Elytra P3", "Elytra P4", "Elytra P5"))

# plot data sorted by primer pairs

#dataPool$Sample <- factor(dataPool$Sample, levels = c("dsVer", "dsTsr"))
for (gene_name in genes) {
  tempPlotData <- subset(dataPool, grepl(gene_name, dataPool$Gene))
  plotData(tempPlotData, gene_name, add_statistics)
  # save plot and data
  ggsave(paste0(gene_name, "_Elytra.png"), path = storePath)
  ggsave(paste0(gene_name, "_Elytra.pdf"), path = storePath)
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

