library(ggplot2)

NormalizationStag1$Ct <- as.numeric(NormalizationStag1$Ct)
NormalizationStag1 <- na.omit(NormalizationStag1)
ans <- aggregate(NormalizationStag1$Ct, list(NormalizationStag1$Gene), FUN=mean)

ans <- ans[3,2]

anselm <- c()
ans3 <- c()

for (i in 1:nrow(NormalizationStag1)) {
  if (isTRUE(grepl("RpS6", NormalizationStag1$Gene[i]))) {
    ans2 <- ((ans - NormalizationStag1$Ct[i])*log(2))
    anselm <- c(anselm, ans2)
    
  }
  
}
anselm <- c(anselm, anselm, anselm)
anselma <- c()
for (i in 1:nrow(NormalizationStag1)) {
  ans3 <- NormalizationStag1$Ct[i] + anselm[i]
  anselma <- c(anselma, ans3)
}

NormalizationStag1 <- na.omit(NormalizationStag1)
CHS1 <- subset(NormalizationStag1, Gene == "CHS1")
Exp <- subset(NormalizationStag1, Gene == "Exp")
RpS6 <- subset(NormalizationStag1, Gene == "RpS6")

ans1 <- aggregate(CHS1$Ct, list(CHS1$StageRep), FUN=mean)
ans2 <- aggregate(Exp$Ct, list(Exp$StageRep), FUN=mean)
ans3 <- aggregate(RpS6$Ct, list(RpS6$StageRep), FUN=mean)
ansa <- rbind(ans1, ans2, ans3)
ida <- rep("CHS1", 9)
idb <- rep("Exp", 9)
idc <- rep("RpS6", 9)
ids <- matrix(c(ida, idb, idc))
meanData <- cbind(meanData, ids)
Stages <- c(rep("P0", 3), rep("P1", 3), rep("P2", 3))
Stages <- matrix(rep(Stages, 3))
meanData <- cbind(meanData, Stages)
meanData <- data.frame(
  Ct = ansa$x,
  Gene = ids,
  Stage = Stages
)


ggplot(data = meanData, aes(x= Stage, y= Ct, fill= Gene)) + geom_boxplot()  

ggplot(data = NormalizationStag1, aes(x = Stage, y = Ct, fill = Gene)) + geom_boxplot()