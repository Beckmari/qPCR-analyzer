#install and load dependencies
source("qPCRfunctions.R")
source("dependencies.R")

# enter unique experiment identifier
expID <- "experimental ID"
# enter your reference gene
RefG <- "RpS18" 
# enter all your genes w/o reference gene
genes <- list("FABP") 
# enter path to prepared data in *.xlsx file
loadPath <- "C:/Users/beckm/OneDrive/Desktop/micPCR/FABP15275/FABP15275_L1-P0.xlsx"
# enter path where data shall be stored
storePath <- "C:/Users/beckm/OneDrive/Desktop/micPCR"
# if you want to show stats, set TRUE
add_statistics <- FALSE

# enter your empirical primerefficancies in the following list separated by comma
# if you don't know your empirical primerefficancies enter a theoretical value of 2
efficancies <- list(
  RpS6 = 1.99, #empiric
  RpS18 = 1.93, #empiric
  FABP = 1.91, #empiric
  CHS1 = 2.00, #empiric
  Exp = 1.95 #empiric
)

#load data
DataTable <- read_excel(loadPath)

#prepare your data
ValueTable <- dataPrep(DataTable, RefG, genes)

#normalize Ct values with Primerefficancies against reference gene and calculate deltaCt
dataPool <- c()
for (gene_name in genes) {
  var_name <- paste0(gene_name, "_plot")
  assign(var_name, plotPrep(efficancies, RefG, gene_name, get(gene_name)))
  dataPool <- c(dataPool, get(var_name))
}
dataPool <- as.data.frame(dataPool)
#check for n < 2 in dataPool
duplicates <- duplicated(dataPool$Sample) | duplicated(dataPool$Sample, fromLast = TRUE)
extractedData <- dataPool[!duplicates, ]
dataPool <- dataPool[duplicates, ]


#plot data
plotData(dataPool, expID, add_statistics)


# save plot and data
ggsave(paste0(expID, "_normalized_boxplot.tiff"), path = storePath)
ggsave(paste0(expID, "_normalized_boxplot.svg"), path = storePath)
write.xlsx(dataPool, paste0(storePath, "/", expID, "_boxplot_Data.xlsx"))
write.xlsx(extractedData, paste0(storePath, "/", expID, "_extractedData.xlsx"))

# clean up only if you're certain to remove all data from the workspace
# saved data will NOT be deleted from the disk
rm(list = ls())

