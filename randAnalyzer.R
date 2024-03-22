#install and load dependencies
source("qPCRfunctions.R")
source("dependencies.R")

# enter unique experiment identifier
expID <- "experiment name"
# enter your reference gene
RefG <- "RpS18" 
# enter all your genes w/o reference gene
genes <- list("genes", "of", "interest") 
# enter path to prepared data in *.xlsx file
loadPath <- "C:/path/to/data.xlsx"
# enter path where data shall be stored
storePath <- "C:/path/to/store/images"
# if you want to show stats, set TRUE
add_statistics <- TRUE

# enter your empirical primerefficancies in the following list separated by comma
# if you don't know your empirical primerefficancies enter a theoretical value of 2
efficancies <- list(
  RpS6 = 1.99, #empiric
  RpS18 = 1.93, #empiric
  FABP = 1.91, #empiric
  CHS1 = 2.00, #empiric
  Exp = 1.95, #empiric
  CDA5 = 2.00, #theoretic
  CDA6 = 2.026, #empiric
  CDA7 = 2.075, #empiric
  CDA8 = 2.00, #theoretic
  CDA9 = 2.00 #theoretic
)

# load data
DataTable <- read_excel(loadPath)

# prepare your data
ValueTable <- dataPrep(DataTable, RefG, genes)

# normalize Ct values with Primerefficancies against reference gene and calculate deltaCt
uni_samples <- unique(ValueTable$Sample)

dataPool <- c()
for (sample_name in uni_samples) {
  var_name <- paste0(sample_name, "_plot")
  for (gene_name in genes) {
    tempData <- plotPrep(efficancies, RefG, gene_name, get(sample_name))
    dataPool <- rbind(dataPool, tempData)
  }
}
dataPool <- as.data.frame(dataPool)
# check for n < 2 in dataPool
duplicates <- duplicated(dataPool$Sample) | duplicated(dataPool$Sample, fromLast = TRUE)
extractedData <- dataPool[!duplicates, ]
dataPool <- dataPool[duplicates, ]


# plot data
for (gene_name in genes) {
  tempPlotData <- subset(dataPool, grepl(gene_name, dataPool$Gene))
  plotData(tempPlotData, gene_name, add_statistics)
  # save plot and data
  ggsave(paste0(gene_name, "_normalized_boxplot.tiff"), path = storePath)
  ggsave(paste0(gene_name, "_normalized_boxplot.svg"), path = storePath)
}

# save data
saveData(dataPool, extractedData, storePath, expID)



# clean up only if you're certain to remove all data from the workspace
# saved data will NOT be deleted from the disk
rm(list = ls())

