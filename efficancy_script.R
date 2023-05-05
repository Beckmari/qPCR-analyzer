effend <- which(ctvalues$Well == "G12")
efficancy <- ctvalues[1:effend,]
efficancy <- cbind(efficancy$Well, efficancy$Ct)
efficancy <- transform(efficancy, X2 = as.numeric(X2))
efficancystat <- c()
eff_rows <- nrow(efficancy)
for (n in seq(from = 1, to = effend , by = 3)) {
  a <- template$template[which(template$plate == efficancy$X1[n])]
  nend <- n+2
  if (!is.na(efficancy[n, 2])){
    i <- n+1
    if (!is.na(efficancy[i,2])){
      if (!is.na(efficancy[nend,2])){
        b <- efficancy[n:nend,2]
        k <- mean(efficancy[n:nend,2])
        k_sd <- sd(efficancy[n:nend,2])
        k_var <- var(efficancy[n:nend,2])
        result <- c(a, b, k, k_sd, k_var)
        efficancystat <- cbind(efficancystat, result)
        rm(result)
      }
    }
  }
  
}

efficancystat <- efficancystat[, -c(1,3)]
efficancystat <- t(efficancystat)
efficancystat <- data.frame(efficancystat)
effCTVal <- efficancystat[, 2:5]
effCTVal <- apply(effCTVal, 2, function(x) as.numeric(as.character(x)))
effConc <- efficancystat$X1
effConc <- as.numeric(effConc)
effConc <- effConc * 50
effConc <- log10(effConc)
tryout <- as.vector(t(effCTVal))
effConc2 <- c(effConc, effConc, effConc, effConc)
effConc2 <- sort(effConc2, decreasing = TRUE)
df <- cbind(effConc2, tryout)
df <- data.frame(df)
fm <- lm(df$effConc2 ~ df$tryout)
plot(df$effConc2, df$tryout)
abline(fm, col = "red")



relefficancy <- 10^(-1/fm$coefficients[2])




#groups in list are needed to be separated in each loop to determine whether they belong together

