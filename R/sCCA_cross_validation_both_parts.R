sCCA_cross_validation_both_parts <- function(X1, Y1, X2, Y2, method, lambda, selection,
                                             alpha, select.covariance, search.grid, same.threshold, previous){

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

  p1_prime <- ncol(X2)
  p2_prime <- ncol(Y2)
  pp_prime <- p1_prime + p2_prime
  n_prime <- nrow(X2)
  training_size_prime <- floor(2*n_prime/3)
  training_prime <- sample(1:n_prime, training_size_prime)
  eta_prime <- sqrt(log(pp_prime)/n_prime)

  trainingX_prime <- X2[training_prime, ]
  trainingY_prime <- Y2[training_prime, ]
  tuneX_prime <- X2[-training_prime, ]
  tuneY_prime <- Y2[-training_prime, ]

  res <- sCCA_zero_padding_cross_validation(trainingX_prime, tuneX_prime, trainingY_prime, tuneY_prime, preX, preY,
                                            p1_prime, pp_prime, search.grid, eta_prime, same.threshold, previous)
  return(res)
}
