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
      selCols <- c(selCols, "Sample")
      #split_txt <- strsplit(df$Sample, " ")
      #nsplt_text <- length(split_txt[[1]])
      #split_txt <- matrix(unlist(split_txt), ncol = nsplt_text, byrow = TRUE)
      #df <- cbind(df, split_txt)
      #colnames(df)[c(ncol(df)-1, ncol(df))] <- c("Gene", "SampleName")
      #selCols <- c(selCols, "Gene", "SampleName")
    } else if (grepl("Gene", col_name)) {
      selCols <- c(selCols, col_name)
    } 
  }
  df <- cbind(df, paste(df$Sample, df$Gene))
  colnames(df)[ncol(df)] <- "Sample_Gene"
  selCols <- c(selCols, "Sample_Gene")
  num_selCols <- length(selCols)
  num_selCols2 <- num_selCols+1
  num_selCols3 <- num_selCols+num_selCols
  ValueTable <- df[,selCols]
  ValueTable$Ct <- as.numeric(gsub(",", ".", ValueTable$Ct))
  
  uni_samples <- unique(ValueTable$Sample)
  
  
  for (sample_name in uni_samples) {
    var_name <- sample_name
    assign(var_name, subset(ValueTable, grepl(sample_name, ValueTable$Sample)))
    assign(paste0(RefG, "-", sample_name), subset(get(var_name), grepl(RefG, get(var_name)$Gene)))
  }
  for (sample_name in uni_samples) {
    tempStore <- get(sample_name)
    tempVar <- c()
    for (gene_name in genes) {
      if (!grepl(RefG, gene_name)) {
        if (any(grepl(gene_name, tempStore$Gene))) {
          var_name2 <- gene_name
          datVar <- subset(tempStore, grepl(var_name2, tempStore$Gene))
          datVar <- cbind(datVar, get(paste0(RefG, "-", sample_name)))
          datVar[datVar < 0] <- NA
          datVar <- na.omit(datVar)
          tempVar1 <- datVar[,1:num_selCols]
          #tempVar2 <- rbind(tempVar1, datVar[,num_selCols2:num_selCols3])
          tempVar <- rbind(tempVar, tempVar1)
        }
      }
    }
    tempVar <- rbind(tempVar, datVar[,num_selCols2:num_selCols3])
    assign(sample_name, tempVar, envir = .GlobalEnv)
  }
  
  return(ValueTable)
}
  
plotPrep <- function(eff, RefG, gene, sample_values) {
  
  plotData <- subset(sample_values, grepl(gene, sample_values$Gene))
  Ref <- subset(sample_values, grepl(RefG, sample_values$Gene))
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
  bxp <- ggplot(data = dataPool, aes(x = reorder(Sample_Gene, -normCt), y = normCt)) + geom_boxplot() 
  if (add_statistics == TRUE) {
    stat.test <- qPCRstats(dataPool)
    stat.test <- stat.test %>% add_xy_position(x = "Sample_Gene")
    bxp + 
      stat_pvalue_manual(stat.test, hide.ns = TRUE, label = "p.adj.signif", tip.length = 0, step.increase = 0.1) +
      labs(
        title = expID,
        subtitle = get_test_label(res_aov, detailed = TRUE), 
        caption = get_pwc_label(res_posthoc),
        x = "Sample",
        y = paste0("relative Expression of ", expID, " [normalized /", RefG, "]")
        ) + 
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    bxp + 
      labs(
        title = expID,
        x = "Sample",
        y = paste0("relative Expression of ", expID, " [normalized /", RefG, "]")
      )
  }
}

qPCRstats <- function(data) {
  shapVal <- data %>% shapiro_test(normCt)
  levVal <- data %>% levene_test(normCt ~ Sample_Gene)
  if (shapVal$p > 0.05 & levVal$p > 0.05) {
    stat.test <- data %>% anova_test(normCt ~ Sample_Gene)
    posthoc_res <- data %>% pairwise_t_test(normCt ~ Sample_Gene, p.adjust.method = "bonferroni")
  } else {
    stat.test <- data %>% kruskal_test(normCt ~ Sample)
    posthoc_res <- data %>% pairwise_t_test(normCt ~ Sample_Gene, p.adjust.method = "bonferroni")
  }
  assign("res_aov", stat.test, envir = .GlobalEnv)
  assign("res_posthoc", posthoc_res, envir = .GlobalEnv)
}

saveData <- function(dataPool, extractedData, storePath, expID) {
  write.xlsx(dataPool, paste0(storePath, "/", expID, "_boxplot_Data.xlsx"))
  if (nrow(extractedData) > 0) {
    write.xlsx(extractedData, paste0(storePath, "/", expID, "_extractedData.xlsx"))
  }
}












































