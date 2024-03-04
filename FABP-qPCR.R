library(ggplot2)
library(dplyr)

#effRPS6 <- 1.99 #empirical
#effCHS1 <- 2.00 #empirical
#effExp <- 1.95 #empirical
effRpS18 <- 2 #theoretical
effFABP15275 <- 2 #theoretical


#ValueTable$Cq <- as.numeric(ValueTable$Cq)

FABP <- subset(ValueTable, grepl("FABP", ValueTable$Sample))
FABP <- subset(FABP, !grepl("H", FABP$Well))
FABP$Cq <- as.numeric(FABP$Cq)
FABP$Result <- NULL
#FABP <- na.omit(FABP)

RpS18 <- subset(ValueTable, grepl("RpS18", ValueTable$Sample))
RpS18 <- subset(RpS18, !grepl("H", RpS18$Well))
RpS18$Cq <- as.numeric(RpS18$Cq)
RpS18$Result <-NULL
#RpS18 <- na.omit(RpS18)

effCq <- effFABP15275^FABP$Cq
FABP <- cbind(FABP, effCq)
effCq <- effRpS18^RpS18$Cq
RpS18 <- cbind(RpS18, effCq)

nFABP <- aggregate(effCq ~ Sample, data = FABP, FUN = mean)
sdv <- aggregate(effCq ~ Sample, data = FABP, FUN = sd)
nFABP <- cbind(nFABP, sdv)
nRpS18 <- aggregate(effCq ~ Sample, data = RpS18, FUN = mean)
sdv <- aggregate(effCq ~ Sample, data = RpS18, FUN = sd)
nRpS18 <- cbind(nRpS18, sdv)

normCq <- RpS18$effCq/FABP$effCq
deltaCq <- RpS18$Cq - FABP$Cq
FABP <- cbind(FABP, normCq, deltaCq)

ggplot(data = FABP, aes(x = Sample, y = deltaCq)) + geom_boxplot()

