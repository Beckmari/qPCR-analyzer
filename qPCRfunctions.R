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

plotData <- function(data, expID, add_statistics = FALSE) {
  #bxp <- ggplot(data = data, aes(x = reorder(desc(Sample_Gene), desc(Gene)), y = normCt)) + 
    #geom_boxplot()
    
  bxp <- ggplot(data = data, aes(x = Sample, y = normCt)) + geom_boxplot(fill = "lightblue") 
  #bxp <- ggplot(data = data, aes(x = Gene, y = relExp)) + geom_boxplot() 
  if (add_statistics == TRUE) {
    stat.test <- qPCRstats(data)
    stat.test <- stat.test %>% add_xy_position(x = "Sample_Gene")
    bxp + 
      stat_pvalue_manual(stat.test, hide.ns = TRUE, label = "p.adj.signif", tip.length = 0, step.increase = 0.1, y.position = 0.00105) +
      labs(
        title = expID,
        subtitle = get_test_label(res_aov, detailed = TRUE), 
        caption = get_pwc_label(res_posthoc),
        x = "Sample",
        y = paste0("mean normalized Expression \n of ", expID, " [normalized /", RefG, "]")
        ) + 
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 20),#angle = 45, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20)
            ) + 
      theme_minimal() +
      ylim(0.000375, 0.00105)
  } else {
    bxp + 
      labs(
        title = expID,
        x = "Sample",
        y = paste0("mean normalized Expression \n of ", expID, " [normalized /", RefG, "]")
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 20),#angle = 45, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title =  element_text(size = 20)
      ) +
      ylim(0, 0.16) +
      theme_minimal()
  }
  return(bxp)
}

qPCRstats <- function(data) {
  #data$normCt <- data$normCt * 100000 ###########################
  #data$relExp <- data$relExp * 100000 ###########################
  shapVal <- data %>% shapiro_test(normCt)
  #shapVal <- data %>% shapiro_test(relExp)
  #shapVal <- c()
  #shapVal$p <- 0.05
  levVal <- data %>% levene_test(normCt ~ Sample_Gene)
  #levVal <- data %>% levene_test(relExp ~ Sample_Gene)
  if (shapVal$p > 0.05 & levVal$p > 0.05) {
    stat.test <- data %>% anova_test(normCt ~ Sample_Gene)
    posthoc_res <- data %>% pairwise_t_test(normCt ~ Sample_Gene, p.adjust.method = "bonferroni")
  } else {
    stat.test <- data %>% kruskal_test(normCt ~ Sample)
    posthoc_res <- data %>% pairwise_t_test(normCt ~ Sample_Gene, p.adjust.method = "bonferroni")
  }
  #if (shapVal$p > 0.05 & levVal$p > 0.05) {
  #  stat.test <- data %>% anova_test(relExp ~ Sample_Gene)
  #  posthoc_res <- data %>% pairwise_t_test(relExp ~ Sample_Gene, p.adjust.method = "bonferroni")
  #} else {
  #  stat.test <- data %>% kruskal_test(relExp ~ Sample)
  #  posthoc_res <- data %>% pairwise_t_test(relExp ~ Sample_Gene, p.adjust.method = "bonferroni")
  #}
  assign("res_aov", stat.test, envir = .GlobalEnv)
  assign("res_posthoc", posthoc_res, envir = .GlobalEnv)
}

saveData <- function(dataPool, extractedData = NA, storePath, expID) {
  write.xlsx(dataPool, paste0(storePath, "/", expID, "_boxplot_Data.xlsx")) # "_", expGroup,
  
  if (!is.na(extractedData)) {
    write.xlsx(extractedData, paste0(storePath, "/", expID, "_extractedData.xlsx"))
  }
}



##### additional for relative expression

#load data from *.xlsx file with possible different sheets for biological replicates

load_Data <- function(loadPath) {
  
  DataTable <- c()
  biol_rep <- length(excel_sheets(loadPath))
  
  for (sheet in 1:biol_rep) {
    temp <- read_excel(loadPath, sheet = sheet)
    
    genes <- unique(temp$Gene)
    genes <- subset(genes, !grepl(RefG, genes))
    
    samples <- unique(temp$Sample)
    samples <- subset(samples, !grepl(refGroup, samples))
    
    temp$Cq <- gsub(",", ".", temp$Cq)
    temp$Cq <- gsub("-1", NA, temp$Cq)
    temp$Cq <- as.numeric(temp$Cq)
    #temp <- na.omit(temp)
    if (biol_rep >= 2) {
      temp <- aggregate(temp, list(temp$Sample, temp$Gene), mean)
      temp$Sample <- temp$Group.1
      temp$Gene <- temp$Group.2
    }
    DataTable <- rbind(DataTable, temp)
  }
  assign("genes", genes, envir = .GlobalEnv)
  assign("biol_rep", biol_rep, envir = .GlobalEnv)
  return(DataTable)
  
}


# dataPrep and Ct calc for rel Expression

rel_data_prep <- function(DataTable, RefG, genes, refGroup) {
  
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
  
  tempPool <- c()
  tempstatPoolaov <- c()
  tempstatPoolph <- c()
  for (gene in genes) {
    tempgene <- subset(dataPool, grepl(gene, dataPool$Gene))
    temp <- subset(tempgene, grepl(refGroup, tempgene$Sample))
    temp <- na.omit(temp)
    temp <- 1 / median(temp$normCt)
    tempgene$relExp <- tempgene$normCt * temp
    #tempstat <- tempgene %>% t_test(relExp ~ Sample)
    qPCRstatsrel(tempgene)
    #tempstatPoolaov <- rbind(tempstatPoolaov, as.data.frame(res_aov))
    tempstatPoolph <- rbind(tempstatPoolph,  as.data.frame(res_posthoc))
    tempPool <- rbind(tempPool, tempgene)
  }
  dataPool <- distinct(tempPool)
  assign("OA_res_posthoc", tempstatPoolph, envir = .GlobalEnv)
  return(dataPool)
}


# rel plotter

rel_plot <- function(dataPool, bars = "se") {
  
  dataPoolTest <- subset(dataPool, !grepl(refGroup, dataPool$Sample))
  dataPoolTest <- distinct(dataPoolTest)
  dataPoolTest$Gene <- factor(dataPoolTest$Gene, levels = level_list)
  
  KDeffsum <- subset(dataPool, grepl(level_list[[1]], dataPool$Gene))
  KDeffsum <- summarySE(data=KDeffsum, measurevar="relExp", groupvars="Sample", na.rm=FALSE, conf.interval=.95)
  KDeffsum <- subset(KDeffsum, !grepl(refGroup, KDeffsum$Sample))
  
  
  if(biol_rep > 1) {
    
    #plotData(dataPoolTest, expGroup)
    ggplot(dataPoolTest, aes(x=Gene, y = relExp*100)) + 
      geom_hline(aes(yintercept = 100), linetype = "dashed", linewidth = 1.2, colour = "blue")+
      geom_boxplot(fill = "lightblue") +
      theme_minimal() +
      labs(
        title = expID,
        x = "Gene",
        y = "relative Expression in %", 
        #caption = paste0("n = 3; efficancy: ", round((100 - (KDeffsum$relExp * 100)), 2), "% \u00B1 ", round(KDeffsum$sd[1]*100, 2), "%")
        caption = paste0("n = 3; efficancy: ", round(diff(KDeffsum$relExp) * 100, 2), "% \u00B1 ", round(KDeffsum$sd[1]*100, 2), "%")
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 15),# angle = 45),#, size = 20),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.caption = element_text(face = "italic", size = 15, hjust = 0.5)
      ) #+
      #ylim(0,100)
  } else if (biol_rep < 2) {
    
    dPbtw <- summarySE(data=dataPoolTest, measurevar="relExp", groupvars="Sample_Gene", na.rm=FALSE, conf.interval=.95)
    dPbtw$Sample_Gene <- gsub(paste0(expGroup, " "), "", dPbtw$Sample_Gene)
    #dPbtw$Sample_Gene <- genes
    dPbtw$Sample_Gene <- factor(dPbtw$Sample_Gene, levels = level_list)
    
    ggplot(dPbtw, aes(x=Sample_Gene, y = relExp*100)) + 
      geom_hline(aes(yintercept = 100), linetype = "dashed", linewidth = 1.2, colour = "blue")+
      geom_errorbar(aes(ymin=(relExp-get(bars))*100, ymax=(relExp+get(bars))*100),width=.15) +
      geom_point(size = 5) +
      theme_minimal() +
      labs(
        title = expID,
        x = "Gene",
        y = "relative Expression",
        caption = paste0("RNAi efficancy: ", round((100 - (KDeffsum$relExp * 100)), 2), "% \u00B1 ", round(KDeffsum$sd[1]*100, 2), "%")
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 15),# angle = 45),#, size = 20),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.caption = element_text(face = "italic", size = 10, hjust = 0.5)
      ) #+
      #ylim(0,100)
    
  
    }
  
  
  
  
  
  
}




qPCRstatsrel <- function(data) {
  #data$normCt <- data$normCt * 100000 ###########################
  #data$relExp <- data$relExp * 100000 ###########################
  #shapVal <- data %>% shapiro_test(normCt)
  shapVal <- data %>% shapiro_test(relExp)
  #shapVal <- c()
  #shapVal$p <- 0.05
  #levVal <- data %>% levene_test(normCt ~ Sample_Gene)
  levVal <- data %>% levene_test(relExp ~ Sample_Gene)
  #if (shapVal$p > 0.05 & levVal$p > 0.05) {
  #  stat.test <- data %>% anova_test(normCt ~ Sample_Gene)
  #  posthoc_res <- data %>% pairwise_t_test(normCt ~ Sample_Gene, p.adjust.method = "bonferroni")
  #} else {
  #  stat.test <- data %>% kruskal_test(normCt ~ Sample)
  #  posthoc_res <- data %>% pairwise_t_test(normCt ~ Sample_Gene, p.adjust.method = "bonferroni")
  #}
  if (shapVal$p > 0.05 & levVal$p > 0.05) {
    stat.test <- data %>% anova_test(relExp ~ Sample_Gene)
    posthoc_res <- data %>% pairwise_t_test(relExp ~ Sample_Gene, p.adjust.method = "bonferroni")
  } else {
    stat.test <- data %>% kruskal_test(relExp ~ Sample)
    posthoc_res <- data %>% pairwise_t_test(relExp ~ Sample_Gene, p.adjust.method = "bonferroni")
  }
  assign("res_aov", stat.test, envir = .GlobalEnv)
  assign("res_posthoc", posthoc_res, envir = .GlobalEnv)
}






























