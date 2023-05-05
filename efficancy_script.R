effend <- which(ctvalues$Well == "G12")
efficancy <- ctvalues[1:effend,]
my_list <- list()
efficancy <- cbind(efficancy$Well, efficancy$Ct)
efficancy <- transform(efficancy, X2 = as.numeric(X2))
for (n in 1:nrow(efficancy)) {
  nend <- n+2
  k <- mean(efficancy[n:nend,2])
  k_sd <- sd(efficancy[n:nend,2])
  k_var <- var(efficancy[n:nend,2])
  list_var <- c(efficancy[n:nend, 2], k, k_sd, k_var)
  my_list <- c(my_list, list_var) #try append
  n <- n+2
}
#groups in list are needed to be seperated in each loop to determine wether they belong together
