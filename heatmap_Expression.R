#Expression data collection of different genes compared to CHS1 as Heatmap

#load dependency collection 
source("dependencies.R")

#load data collection
pDataCollect <- "C:/Users/Marius Beck/OneDrive/Desktop/ExpressionProfile/ProfileDataElytra.xlsx"
DataCollect <- read_excel(pDataCollect)

#load new data and save in Collection
pNewData <- "C:/Users/Marius Beck/OneDrive/Desktop/ExpressionProfile/TsrProfileE_boxplot_Data.xlsx"
NewData <- read_excel(pNewData)
DataCollect <- rbind(DataCollect, NewData)
write.xlsx(DataCollect, pDataCollect)

#create heatmap

ggplot(data = DataCollect, aes(Sample, Gene, fill = deltaCt)) + 
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
  theme(legend.position = "bottom") 
  #scale_fill_gradient(low="darkblue", high = "yellow") + 
  theme_ipsum()

datmat2 <- as.matrix(t(datMat))
heatmap(datmat2, Colv = NA, scale = "row")











theme(legend.position = "bottom")theme(legend.position = "bottom")test <- c()

test <- DataCollect %>%
  group_by(Sample, Gene) %>%
  summarise(across(contains("Ct") ,mean))

test2 <- subset(test, !grepl("P5", test$Sample))

datMat <- c()
datMat$CHS1 <- test2$normCt[test2$Gene == "CHS1"]
datMat$Tsr <- test2$normCt[test2$Gene == "Tsr"]
datMat$Ctl2_1323 <- test2$normCt[test2$Gene == "Ctl2_1323"]
datMat$Osiris18 <- test2$normCt[test2$Gene == "Osiris18"]
datMat$FABP15275 <- test2$normCt[test2$Gene == "FABP15275"]
datMat$SERCA12671 <- test2$normCt[test2$Gene == "SERCA12671"]
datMat <- as.data.frame(datMat)
rownames(datMat) <- unique(test2$Sample)
corr <- cor(datMat)
p.mat2 <- cor_pmat(datMat)

#ggcorrplot(corr, hc.order = TRUE, p.mat = p.mat, type = "lower", method = "circle")

corr2 <- corr.test(datMat)
corr <- corr2$r
p.mat <- corr2$p
col1 <- colorRampPalette(c( "#00007F",  "blue", "#007FFF",
                            "cyan", "white", "yellow", "#FF7F00", "red","#7F0000"))
p <- ggcorrplot.mixed(corr, upper = "ellipse", lower = "number", p.mat = p.mat2, insig = "label_sig", sig.lvl = c(0.05, 0.01, 0.001), pch = "*", pch.cex = 4) +
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



r <- cor(x = test$normCt[test$Gene == "Tsr"], y = test$normCt[test$Gene == "Osiris18"], method = "pearson")
R2 <- r^2
r
R2

ggplot() +
  geom_point(aes(x = test2$normCt[test2$Gene == "Tsr"], y = test2$normCt[test2$Gene == "Ctl2_1323"]), color = "red") 
  geom_point(aes(x = test2$normCt[test2$Gene == "Tsr"], y = test2$normCt[test2$Gene == "Osiris18"]), color = "blue")+
  geom_point(aes(x = test2$normCt[test2$Gene == "Tsr"], y = test2$normCt[test2$Gene == "FABP15275"]), color = "yellow")+
  geom_point(aes(x = test2$normCt[test2$Gene == "Osiris18"], y = test2$normCt[test2$Gene == "Ctl2_1323"]), color = "green")+
  geom_point(aes(x = test2$normCt[test2$Gene == "Osiris18"], y = test2$normCt[test2$Gene == "FABP15275"]), color = "violet")+
  geom_point(aes(x = test2$normCt[test2$Gene == "FABP15275"], y = test2$normCt[test2$Gene == "Ctl2_1323"]), color = "black")


r <- cor(x = test2$normCt[test2$Gene == "CHS1"], y = test2$normCt[test2$Gene == "Ctl2_1323"], method = "pearson")
R2 <- r^2
r
R2

ggplot() +
  #geom_point(aes(y = normCt, color = Gene)) +
  #geom_point(aes(x = CtOsi, y = CtFABP), color = "blue") 
  geom_point(aes(x = test2$normCt[test2$Gene == "CHS1"], y = test2$normCt[test2$Gene == "Ctl2_1323"]), color = "red")
