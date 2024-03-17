#qPCR functions

normqPCR <- function(effRef, effGOI, GOI, Ref) {
  GOI_eff <- effGOI^GOI$Ct
  Ref_eff <- effRef^Ref$Ct
  normCtVal <- Ref_eff/GOI_eff
  GOI <- cbind(GOI, normCtVal)
  colnames(GOI)[ncol(GOI)] <- "normCt"
  return(GOI)
}

deltaCtqPCR <- function(GoI, RefG_CT, GoI_Ct, norm = FALSE) {
  deltaCtVal <- RefG_CT - GoI_Ct
  GoI <- cbind(GoI, deltaCtVal)
  if (norm == FALSE) {
    colnames(GoI)[ncol(GoI)] <- "deltaCt"
  } else if (norm == TRUE) {
    colnames(GoI)[ncol(GoI)] <- "normdeltaCt"
  }
  return(GoI)
}



dataPrep <- function(df, RefG, genes) { 
  #daten sollten vor na.omit noch entsprechend sortiert werden, so dass referenz und sample immer zusammengehören
  num_genes <- length(genes)
  selCols <- c()
  columnNames <- colnames(df)
  for (col_name in colnames(df)) {
    if (grepl("Ct", col_name)) {
      selCols <- c(selCols, col_name)
    } else if (grepl("Cq", col_name)) {
      colnames(df)[colnames(df) == "Cq"] <- "Ct"
      selCols <- selCols <- c(selCols, "Ct")
    } else if (grepl("Sample", col_name)) {
      colnames(df)[colnames(df) == "Sample"] <- "Sample"      
      selCols <- c(selCols, col_name)
      split_txt <- strsplit(df$Sample, " ")
      split_txt <- matrix(unlist(split_txt), ncol = 2, byrow = TRUE)
      df <- cbind(df, split_txt)
      colnames(df)[c(ncol(df)-1, ncol(df))] <- c("Gene", "SampleName")
      selCols <- c(selCols, "Gene", "SampleName")
    } else if (grepl("Stage", col_name)) {
      selCols <- c(selCols, col_name)
    }
  }
  num_selCols <- length(selCols)
  num_selCols2 <- num_selCols+1
  num_selCols3 <- num_selCols+num_selCols
  ValueTable <- df[,selCols]
  ValueTable$Ct <- as.numeric(gsub(",", ".", ValueTable$Ct))
  assign(RefG, subset(ValueTable, grepl(RefG, ValueTable$Sample)), envir = .GlobalEnv)
  for (gene_name in genes) {
    if (!grepl(RefG, gene_name)) {
      var_name <- gene_name
      datVar <- subset(ValueTable, grepl(var_name, ValueTable$Sample))
      datVar <- cbind(datVar, get(RefG))
      datVar[datVar < 0] <- NA
      datVar <- na.omit(datVar)
      tempVar <- datVar[,1:num_selCols]
      tempVar <- rbind(tempVar, datVar[,num_selCols2:num_selCols3])
      assign(var_name, tempVar, envir = .GlobalEnv)
    }
  }
  return(ValueTable)
}
  
plotPrep <- function(eff, RefG, gene, gene_values) {
  
  plotData <- subset(gene_values, grepl(gene, gene_values$Sample))
  Ref <- subset(gene_values, grepl(RefG, gene_values$Sample))
  effRef <- eff[[RefG]]
  effGOI <- eff[[gene]]
  
  plotData <- normqPCR(effRef, effGOI, plotData, Ref)
  # standard delta Ct
  plotData <- deltaCtqPCR(plotData, Ref$Ct, plotData$Ct)
  # normalized delta Ct
  #plotData <- deltaCtqPCR(plotData, Ref$normCt, plotData$normCt, norm = TRUE)
  
  
  return(plotData)
}

plotData <- function(dataPool, expID, add_statistics = FALSE) {
  stat.test<- dataPool %>% t_test(normCt ~ Sample)
  bxp <- ggplot(data = dataPool, aes(x = Sample, y = normCt)) + geom_boxplot() + labs(title = expID) + 
    theme(plot.title = element_text(hjust = 0.5))
  if (add_statistics == TRUE) {
    stat.test <- stat.test %>% add_xy_position(x = "Sample")
    bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01)
  } else {
    bxp
  }
}


