#Script function library

#delta Ct calculation
dCt_calc <- function(a, b) {
  dCt <- a - b #experimental Gene Ct - Referencegene Ct
  return(dCt)
}


#delta delta Ct
ddCt_calc <- function(a, b) {
  ddCt <- a - b #experimental gene dCt - Positivecontrol dCt
  rel_ddCt <- 2^-ddCt #expression in %
  ddCt_val <- list("ddCt" = ddCt, "rel_ddCt" = rel_ddCt)
  return(ddCt_val)
}