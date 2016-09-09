covariance_toeplitz <- function(X, tuneX, select.covariance = "loglikelihood"){
  
  loss_likelihood <- function(Sigma, Omega){
    tmp <- sum(diag(Sigma%*%Omega)) - log(det(Omega))
    if(is.finite(tmp) & tmp > 0) return(tmp)
    else tmp <- 10000
    return(tmp)
  }
  
  floss <- function(x, y){
    return(sum((c(x-y))^2))
  }
  covX <- cov(X)
  covtuneX <- cov(tuneX)
  
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
  
  p <- ncol(X)
  largestk <- p 
  klist <- seq(p, 2, by=-2)
  floss_list <- rep(0, length(klist))
  for(i in 1:length(klist)){
    estCov <- toeplitz_estimator(covX, klist[i])  
    if(select.covariance == "loglikelihood"){
      floss_list[i] <- loss_likelihood(covtuneX, solve(estCov))
    }
    else{
      floss_list[i] <- floss(covtuneX, estCov)
    }
  }
  selected <- klist[which.min(floss_list)]
  estCov <- toeplitz_estimator(covX, selected)
  return(estCov)
}