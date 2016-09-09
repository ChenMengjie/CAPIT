sCCA_zero_padding <- function(A_tuta, lambda1, lambda2, eta, previous = TRUE){

  normalize <- function(vec){
    vec_l2norm <- apply(vec, 2, function(z){
      return(sqrt(sum(z^2)))
    })
    stan_vec <- matrix(0, nrow = nrow(vec), ncol = ncol(vec))
    for(i in 1:ncol(vec)) stan_vec[, i] <- vec[, i]/vec_l2norm[i]
    return(stan_vec)
  }

  #### zero-padding ####

  maxI <- apply(A_tuta, 1, max)
  maxJ <- apply(A_tuta, 2, max)

  Is <- which(abs(maxI) > lambda1*eta)
  Js <- which(abs(maxJ) > lambda1*eta)
  if(length(Is)==0){ Is <- order(maxI, decreasing=T)[1:10] }
  if(length(Js)==0){ Js <- order(maxJ, decreasing=T)[1:10] }
  A_tuta_IJ <- A_tuta[Is, Js]

  A_tuta_IJ_svd <- svd(A_tuta_IJ)
  alpha_0 <- rep(0, length(maxI))
  alpha_0[Is] <- A_tuta_IJ_svd$u[, 1]
  beta_0 <- rep(0, length(maxJ))
  beta_0[Js] <- A_tuta_IJ_svd$v[, 1]

  P_alpha_0 <- alpha_0%*%t(alpha_0)
  P_beta_0 <- beta_0%*%t(beta_0)
  A_tuta_prime <- t(A_tuta)

  if(previous == TRUE){
    for(i in 1:5){
      #### left thresholding
      alpha_1 <- A_tuta%*%beta_0
      alpha_1 <- ifelse(abs(alpha_1) > lambda2*eta, alpha_1, 0)
      if(any(alpha_1 != 0)){
        alpha_1 <- normalize(alpha_1)
      } else {
        alpha_1 <- rep(0, length(alpha_1))
      }

      beta_1 <- A_tuta_prime%*%alpha_0
      beta_1 <- ifelse(abs(beta_1) > lambda2*eta, beta_1, 0)
      if(any(beta_1 != 0)){
        beta_1 <- normalize(beta_1)
      } else {
        beta_1 <- rep(0, length(beta_1))
      }

      alpha_0 <- alpha_1
      beta_0 <- beta_1
    }
  } else {

    tol <- 1

    while(tol > 10^-6){
      #### left thresholding
      alpha_1 <- A_tuta%*%beta_0
      alpha_1 <- ifelse(abs(alpha_1) > lambda2*eta, alpha_1, 0)
      if(any(alpha_1 != 0)){
        alpha_1 <- normalize(alpha_1)
      } else {
        alpha_1 <- rep(0, length(alpha_1))
      }

      beta_1 <- A_tuta_prime%*%alpha_1
      beta_1 <- ifelse(abs(beta_1) > lambda2*eta, beta_1, 0)
      if(any(beta_1 != 0)){
        beta_1 <- normalize(beta_1)
      } else {
        beta_1 <- rep(0, length(beta_1))
      }
      if(all(alpha_1 == 0) | all(beta_1 == 0)) {
        tol <- 10^-7
      } else {
        P_alpha_1 <- alpha_1%*%t(alpha_1)
        P_beta_1 <- beta_1%*%t(beta_1)
        tol <- max(sum(c(P_alpha_1-P_alpha_0)^2), sum(c(P_beta_1-P_beta_0)^2))
        alpha_0 <- alpha_1
        beta_0 <- beta_1
        P_alpha_0 <- P_alpha_1
        P_beta_0 <- P_beta_1
      }
    }

  }

  D <- t(alpha_1)%*%A_tuta%*%beta_1
  if(D < 0){
    beta_1 <- -beta_1
  }
  res <- list(alpha_1, beta_1)
  names(res) <- c("alpha", "beta")
  return(res)
}

