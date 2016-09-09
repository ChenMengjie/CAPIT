sCCA <- function(X1, Y1, X2, Y2, lambda1, lambda2, method,
                            lambda, selection, alpha, select.covariance){
  normalize <- function(vec){
    vec_l2norm <- apply(vec, 2, function(z){
      return(sqrt(sum(z^2)))
    })
    stan_vec <- matrix(0, nrow = nrow(vec), ncol = ncol(vec))
    for(i in 1:ncol(vec)) stan_vec[, i] <- vec[, i]/vec_l2norm[i]
    return(stan_vec)
  }

  p1 <- ncol(X1)
  p2 <- ncol(Y1)
  pp <- p1 + p2
  n <- nrow(X1)
  training_size <- floor(2*n/3)
  training <- sample(1:n, training_size)
  eta <- sqrt(log(pp)/n)

  trainingX <- X1[training, ]
  trainingY <- Y1[training, ]
  tuneX <- X1[-training, ]
  tuneY <- Y1[-training, ]

  if(method == "thresholding"){
    covX <- covariance_thresholding(trainingX, tuneX, lambda = lambda, selection = selection, select.covariance = select.covariance)
    covY <- covariance_thresholding(trainingY, tuneY, lambda = lambda, selection = selection, select.covariance = select.covariance)
    preX <- solve(covX)
    preY <- solve(covY)
  }
  if(method == "toeplitz"){
    covX <- covariance_toeplitz(trainingX, tuneX, select.covariance = select.covariance)
    covY <- covariance_toeplitz(trainingY, tuneY, select.covariance = select.covariance)
    preX <- solve(covX)
    preY <- solve(covY)
  }
  if(method == "tapering"){
    covX <- covariance_tapering(trainingX, tuneX, alpha = alpha, select.covariance = select.covariance)
    covY <- covariance_tapering(trainingY, tuneY, alpha = alpha, select.covariance = select.covariance)
    preX <- solve(covX)
    preY <- solve(covY)
  }

  trans.X2 <- apply(X2, 1, function(z){ preX%*%z })
  trans.Y2 <- apply(Y2, 1, function(z){ preY%*%z })

  Z <- cbind(t(trans.X2), t(trans.Y2))
  cov.Z <- cov(Z)
  A_tuta <- cov.Z[1:p1, (p1+1):(pp)]

  #### zero-padding ####

  maxI <- apply(A_tuta, 1, max)
  maxJ <- apply(A_tuta, 2, max)

  Is <- which(abs(maxI) > lambda1*eta)
  Js <- which(abs(maxJ) > lambda1*eta)
  if(length(Is)==0){ Is <- order(maxI, decreasing=T)[1:10] }
  if(length(Js)==0){ Js <- order(maxJ, decreasing=T)[1:10] }
  A_tuta_IJ <- A_tuta[Is, Js]

  A_tuta_IJ_svd <- svd(A_tuta_IJ)
  alpha_0 <- rep(0, p1)
  alpha_0[Is] <- A_tuta_IJ_svd$u[, 1]
  beta_0 <- rep(0, p2)
  beta_0[Js] <- A_tuta_IJ_svd$v[, 1]

  P_alpha_0 <- alpha_0%*%t(alpha_0)
  P_beta_0 <- beta_0%*%t(beta_0)
  A_tuta_prime <- t(A_tuta)
  tol <- 1
  while(tol > 10^-6){
    #### left thresholding
    alpha_1 <- A_tuta%*%beta_0
    alpha_1 <- ifelse(abs(alpha_1) > lambda2*eta, alpha_1, 0)
    alpha_1 <- normalize(alpha_1)

    beta_1 <- A_tuta_prime%*%alpha_1
    beta_1 <- ifelse(abs(beta_1) > lambda2*eta, beta_1, 0)
    beta_1 <- normalize(beta_1)

    P_alpha_1 <- alpha_1%*%t(alpha_1)
    P_beta_1 <- beta_1%*%t(beta_1)
    tol <- max(sum(c(P_alpha_1-P_alpha_0)^2), sum(c(P_beta_1-P_beta_0)^2))
    alpha_0 <- alpha_1
    beta_0 <- beta_1
    P_alpha_0 <- P_alpha_1
    P_beta_0 <- P_beta_1
  }
  D <- t(alpha_1)%*%A_tuta%*%beta_1
  if(D < 0){
    beta_1 <- -beta_1
  }
  res <- list(alpha_1, beta_1)
  names(res) <- c("alpha", "beta")
  return(res)
}
