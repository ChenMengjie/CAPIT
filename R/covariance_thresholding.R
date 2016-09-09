covariance_thresholding <- function(X, tuneX, lambda, selection, select.covariance){
   
   loss_likelihood <- function(Sigma, Omega){
      tmp <- sum(diag(Sigma%*%Omega)) - log(det(Omega))
      if(is.finite(tmp) & tmp > 0) return(tmp)
      else tmp <- 10000
      return(tmp)
   }
   covX <- cov(X)
   floss <- function(x, y){
     return(sum((c(x-y))^2))
   }
   covtuneX <- cov(tuneX)
   
   floss_list <- rep(0, length(lambda))
   for(i in 1:length(lambda)){
      estCov <- apply(covX, 2, function(x){
          return(Thresholding(x, lambda[i], selection = selection))
      })
      if(select.covariance == "loglikelihood"){
        floss_list[i] <- loss_likelihood(covtuneX, solve(estCov))
      }
      else{
        floss_list[i] <- floss(covtuneX, estCov)
      }
   }
   selected <- lambda[which.min(floss_list)]
   estCov <- apply(covX, 2, function(x){
     return(Thresholding(x, selected, selection = selection))
   })
   return(estCov)
}