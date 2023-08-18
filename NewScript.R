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

ans <- cbind(CHS1, Exp, RpS6)
ans <- na.omit(ans)
CHS1 <- ans[,1:5]
Exp <- ans[,6:10]
RpS6 <- ans[,11:15]

normCt <- RpS6$effCt/CHS1$effCt
CHS1 <- cbind(CHS1, normCt)
normCt <- RpS6$effCt/Exp$effCt
Exp <- cbind(Exp, normCt)
finalValues <- rbind(CHS1, Exp)
#ans <- finalValues[!finalValues$Stage == 'P116',]
ggplot(data = finalValues, aes(x = Stage, y = normCt, fill = Gene)) + geom_boxplot()


ans <- aggregate(normCt ~ Stage + Samplename, data = CHS1, FUN = mean)
