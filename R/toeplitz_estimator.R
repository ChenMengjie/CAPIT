toeplitz_estimator <- function(x, k){
  p <- ncol(x)
  allelements <- c(x)
  for(i in 0:(p-1)){
    omega <- ifelse(i <= k/2, 1, ifelse(i <= k & i > k/2, 2-2*i/k, 0))
    aa <- allelements[seq(1+i, length(allelements) - p*(i-1), by = p + 1)] 
    allelements[seq(1+i, length(allelements) - p*(i-1), by = p + 1)] <- mean(aa)*omega
    allelements[seq(1+p*i, length(allelements), by = p + 1)] <- mean(aa)*omega
  }
  y <- matrix(allelements, p, p)
  return(y)
}