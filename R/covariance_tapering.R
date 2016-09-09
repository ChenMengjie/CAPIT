covariance_tapering <- function(X, tuneX, alpha, select.covariance = "loglikelihood"){
  
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
  p <- ncol(X)
  n <- nrow(X)
  
  tapering_estimator <- function(x, k){
    p <- ncol(x)
    y <- matrix(0, p, p)
    for(i in 1:p){
      for(j in 1:p){
        aa <- abs(i - j) 
        omega <- ifelse(aa <= k/2, 1, ifelse(aa > k/2 & aa < k, 2- aa*2/k, 0))
        y[i, j] <- x[i, j]*omega
      }
    }
    return(y)
  }
  
  optimalk <- floor(n^(1/(2*alpha+1)))  
  klist <- seq(optimalk-20, optimalk+40, by=1)
  floss_list <- rep(0, length(klist))
  for(i in 1:length(klist)){
    estCov <- tapering_estimator(covX, klist[i])
    if(select.covariance == "loglikelihood"){
      floss_list[i] <- loss_likelihood(covtuneX, solve(estCov))
    }
    else{
      floss_list[i] <- floss(covtuneX, estCov)
    }
  }
  selected <- klist[which.min(floss_list)]
  estCov <- tapering_estimator(covX, selected)
  return(estCov)
}