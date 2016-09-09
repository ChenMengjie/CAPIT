tapering_estimator <- function(x, k){
  p <- ncol(x)
  y <- matrix(0, p, p)
  for(i in 1:p){
    aa <- abs(i - 1:p) 
    omega <- ifelse(aa <= k/2, 1, ifelse(aa > k/2 & aa < k, 2- aa*2/k, 0))
    y[i, ] <- x[i, ]*omega
  }
  return(y)
}