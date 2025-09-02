#Expression data collection of different genes compared to CHS1 as Heatmap

#load dependency collection 
source("dependencies.R")

#load data collection
pDataCollect <- "C:/Users/beckm/Desktop/satging/ProfileDataElytra.xlsx"
DataCollect <- read_excel(pDataCollect)

#load new data and save in Collection
pNewData <- "C:/Users/beckm/Desktop/satging/FABP8X3 expression elytra_boxplot_Data.xlsx"
NewData <- read_excel(pNewData)
DataCollect <- rbind(DataCollect, NewData)
write.xlsx(DataCollect, pDataCollect)

#create heatmap
DataCollect$Gene <- factor(DataCollect$Gene, levels = c( "FABP8X3", "FABP8X2", "FABP8X1","CHS1"))

ggplot(data = DataCollect, aes(Sample, Gene, fill = normCt)) + 
  geom_tile() +
  scale_fill_gradientn(colours = col1(10), #limits = c(-1, 1),
                       guide = guide_colorbar(
                         direction = "horizontal",
                         title = "",
                         nbin = 1000,
                         ticks.colour = "black",
                         frame.colour = "black",
                         barwidth = 15,
                         barheight = 1.5)) +
  scale_colour_gradientn(colours = col1(10), #limits = c(-1, 1),
                         guide = guide_colorbar(
                           direction = "horizontal",
                           title = "",
                           nbin = 1000,
                           ticks.colour = "black",
                           frame.colour = "black",
                           barwidth = 15,
                           barheight = 1.5)) +
  
  #scale_fill_gradient(low="darkblue", high = "red") + 
  theme_minimal()+
  theme(legend.position = "bottom") 












test <- c()

test <- DataCollect %>%
  group_by(Sample, Gene) %>%
  summarise_at(vars(contains("Ct")), mean)



test2 <- subset(test, !grepl("P5", test$Sample))

datMat <- c()
datMat$CHS1 <- test2$normCt[test2$Gene == "CHS1"]
datMat$FABP8X1 <- test2$normCt[test2$Gene == "FABP8X1"]
datMat$FABP8X2 <- test2$normCt[test2$Gene == "FABP8X2"]
datMat$FABP8X3 <- test2$normCt[test2$Gene == "FABP8X3"]
#datMat$FABP15275 <- test2$normCt[test2$Gene == "FABP15275"]
#datMat$FABP15275_1 <- test2$normCt[test2$Gene == "FABP15275 1"]
#datMat$SERCA12671 <- test2$normCt[test2$Gene == "SERCA12671"]
#datMat$SERCA12671_1 <- test2$normCt[test2$Gene == "SERCA12671 1"]
datMat <- as.data.frame(datMat)
rownames(datMat) <- unique(test2$Sample)
corr <- cor(datMat)
p.mat2 <- cor_pmat(datMat)

#ggcorrplot(corr, hc.order = TRUE, p.mat = p.mat, type = "lower", method = "circle")

corr2 <- corr.test(datMat)
corr <- corr2$r
p.mat <- corr2$p
col1 <- colorRampPalette(c( "#00007F",  "blue", "#007FFF",
                            "lightblue", "#FF7F00", "red","#7F0000"))
p <- ggcorrplot.mixed(corr, upper = "ellipse", lower = "number", p.mat = p.mat2, insig = "label_sig", sig.lvl = c(0.05, 0.01, 0.001), pch = "*", pch.cex = 9) +
  scale_fill_gradientn(colours = col1(10), limits = c(-1, 1),
    guide = guide_colorbar(
    direction = "horizontal",
    title = "",
    nbin = 1000,
    ticks.colour = "black",
    frame.colour = "black",
    barwidth = 15,
    barheight = 1.5)) +
  scale_colour_gradientn(colours = col1(10), limits = c(-1, 1),
    guide = guide_colorbar(
    direction = "horizontal",
    title = "",
    nbin = 1000,
    ticks.colour = "black",
    frame.colour = "black",
    barwidth = 15,
    barheight = 1.5)) +
  theme(legend.position = "bottom")
p



r <- cor(x = test$normCt[test$Gene == "CHS1"], y = test$normCt[test$Gene == "FABP8X2"], method = "pearson")
R2 <- r^2
r
R2

ggplot() +
  #geom_point(aes(x = test2$normCt[test2$Gene == "CHS1"], y = test2$normCt[test2$Gene == "FABP8X1"]), color = "red") +
  geom_point(aes(x = test2$normCt[test2$Gene == "CHS1"], y = test2$normCt[test2$Gene == "FABP8X2"]), color = "blue")#+
  #geom_point(aes(x = test2$normCt[test2$Gene == "CHS1"], y = test2$normCt[test2$Gene == "FABP8X3"]), color = "yellow")+
  #geom_point(aes(x = test2$normCt[test2$Gene == "FABP8X1"], y = test2$normCt[test2$Gene == "FABP8X2"]), color = "green")+
  #geom_point(aes(x = test2$normCt[test2$Gene == "FABP8X1"], y = test2$normCt[test2$Gene == "FABP8X3"]), color = "violet")+
  #geom_point(aes(x = test2$normCt[test2$Gene == "FABP8X2"], y = test2$normCt[test2$Gene == "FABP8X3"]), color = "black")


r <- cor(x = test2$normCt[test2$Gene == "CHS1"], y = test2$normCt[test2$Gene == "FABP8X2"], method = "pearson")
R2 <- r^2
r
R2

ggplot() +
  #geom_point(aes(y = normCt, color = Gene)) +
  #geom_point(aes(x = CtOsi, y = CtFABP), color = "blue") 
  geom_point(aes(x = test2$normCt[test2$Gene == "CHS1"], y = test2$normCt[test2$Gene == "Ctl2_1323"]), color = "red")
