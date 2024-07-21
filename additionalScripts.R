###extra script for different tasks


##Reference gene stability in tested samples

RefTest <- subset(DataTable, grepl("RpS18", DataTable$Gene))

RefTest$Sample <- factor(RefTest$Sample, levels = c("L1", "L2", "L3", "L4", "L5", "PP", "P0", "P1", "P2", "P3", "P4", "P5", "A"))


#1st try -> with mean
RefTest$Cq <- gsub(",", ".", RefTest$Cq)
RefTest$Cq <- gsub("-1", NA, RefTest$Cq)
RefTest <- na.omit(RefTest)
RefTest$Cq <- as.numeric(RefTest$Cq)
RefTest <- aggregate(RefTest, list(RefTest$Sample), mean)
#RefTest$SD <- sd(RefTest)
RefTest$Sample <- 1:nrow(RefTest) #enter values for 
fm <- lm(RefTest$Cq ~ RefTest$Sample)
cfm <- coef(fm)
sfm <- summary(fm)
r_square <- sfm$adj.r.squared
ggplot(data = RefTest, aes(x = Group.1, y = Cq)) +
  geom_point() + ylim(0, 35) #+
  geom_abline(intercept = cfm[1], slope = cfm[2], color = "red") + ylim(0, 30)

#2nd try -> without mean
RefTest$Number <-  NA
RefTest$Number[which(RefTest$Sample == "Elytra P0")] <-  1
RefTest$Number[which(RefTest$Sample == "Elytra P1")] <-  2
RefTest$Number[which(RefTest$Sample == "Elytra P2")] <-  3
RefTest$Number[which(RefTest$Sample == "Elytra P3")] <-  4
RefTest$Number[which(RefTest$Sample == "Elytra P4")] <-  5
RefTest$Number[which(RefTest$Sample == "Elytra P5")] <-  6
RefTest$Cq <- gsub(",", ".", RefTest$Cq)
RefTest$Cq <- as.integer(RefTest$Cq)
fm <- lm(RefTest$Cq ~ RefTest$Number)
cfm <- coef(fm)
sfm <- summary(fm)
r_square <- sfm$adj.r.squared
ggplot(data = RefTest, aes(x = Sample, y = Cq)) +
  geom_point() +
  geom_abline(intercept = cfm[1], slope = cfm[2], color = "red") + ylim(0, 35)
