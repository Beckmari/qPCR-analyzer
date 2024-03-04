library(ggplot2)
library(dplyr)

effRPS6 <- 1.99
effCHS1 <- 2.00
effExp <- 1.95

ValueTable$Ct <- as.numeric(ValueTable$Ct)

CHS1 <- subset(ValueTable, Gene == "CHS1")
Exp <- subset(ValueTable, Gene == "Exp")
RpS6 <- subset(ValueTable, Gene == "RpS6")

effCt <- effCHS1^CHS1$Ct
CHS1 <- cbind(CHS1, effCt)
effCt <- effExp^Exp$Ct
Exp <- cbind(Exp, effCt)
effCt <- effRPS6^RpS6$Ct
RpS6 <- cbind(RpS6, effCt)

nCHS1 <- aggregate(effCt ~ Samplename + Stage + Gene, data = CHS1, FUN = mean)
nExp <- aggregate(effCt ~ Samplename + Stage + Gene, data = Exp, FUN = mean)
nRpS6 <- aggregate(effCt ~ Samplename + Stage + Gene, data = RpS6, FUN = mean)
ans <- rbind(nCHS1$Samplename, nExp$Samplename, nRpS6$Samplename)

ans <- cbind(CHS1, Exp, RpS6)
ans <- na.omit(ans)
CHS1 <- ans[,1:5]
Exp <- ans[,6:10]
RpS6 <- ans[,11:15]
normCt <- c()
for (i in 1:52) {
  normCtans <- c()
  ans1 <- which(nRpS6 == nCHS1$Samplename[i])
  n <- length(ans1)
  for (x in 1:n) {
    if (grepl(nCHS1$Stage[i], nRpS6$Stage[ans1[x]]) == TRUE) {
      ans <- x
      normCtans <- nRpS6$effCt[i]/nCHS1$effCt[ans]
    }
  }
  if (is.null(normCtans))
    normCtans <- NA
  
  normCt <- c(normCt, normCtans)
}
nCHS1 <- cbind(nCHS1, normCt)

normCt <- c()
for (i in 1:54) {
  normCtans <- c()
  ans1 <- which(nRpS6 == nExp$Samplename[i])
  n <- length(ans1)
  for (x in 1:n) {
    if (grepl(nExp$Stage[i], nRpS6$Stage[ans1[x]]) == TRUE) {
      ans <- x
      normCtans <- nRpS6$effCt[i]/nExp$effCt[ans]
    }
  }
  if (is.null(normCtans))
    normCtans <- NA
  
  normCt <- c(normCt, normCtans)
}
nExp <- cbind(nExp, normCt)

#normCt <- RpS6$effCt/CHS1$effCt
#CHS1 <- cbind(CHS1, normCt)
#normCt <- RpS6$effCt/Exp$effCt
#Exp <- cbind(Exp, normCt)
finalValues <- rbind(nCHS1, nExp)
#ans <- finalValues[!finalValues$Stage == 'P4.16',]
ans <- finalValues[grep("P0", finalValues$Samplename),]
ggplot(data = ans, aes(x = Stage, y = effCt, fill = Gene)) + geom_boxplot()


ans <- aggregate(normCt ~ Stage + Samplename, data = CHS1, FUN = mean)
