sCCA_zero_padding_cross_validation <- function(trainingX, tuneX, trainingY, tuneY, preX, preY,
                                               p1, pp, search.grid, eta, same.threshold = FALSE, previous = TRUE){

  distance_cal <- function(alpha, beta, covX, covY, covXY){
    aa <- t(alpha)%*%covXY%*%beta
    bb <- t(alpha)%*%covX%*%alpha%*%t(beta)%*%covY%*%beta
    return(aa/sqrt(bb))
  }

  trans.X2 <- apply(trainingX, 1, function(z){ preX%*%z })
  trans.Y2 <- apply(trainingY, 1, function(z){ preY%*%z })

  Z <- cbind(t(trans.X2), t(trans.Y2))
  cov.Z <- cov(Z)
  A_tuta <- cov.Z[1:p1, (p1+1):(pp)]

  cov.tuneXY <- cov(tuneX, tuneY)
  cov.tuneX <- cov(tuneX)
  cov.tuneY <- cov(tuneY)

  dis_list <- rep(0, length(search.grid))
  if(same.threshold == TRUE){
    for(i in 1:length(search.grid)){
      lambda1 <- lambda2 <- search.grid[i]
      alpha_beta_pair <- sCCA_zero_padding(A_tuta, lambda1, lambda2, eta, previous)
      if(any(alpha_beta_pair[[1]] != 0) & any(alpha_beta_pair[[2]] != 0)){
        dis_list[i] <- distance_cal(alpha_beta_pair[[1]], alpha_beta_pair[[2]], cov.tuneX, cov.tuneY, cov.tuneXY)
      }
    }
    selected <- search.grid[which.max(abs(dis_list))]
    final_alpha_beta_pair <- sCCA_zero_padding(A_tuta, selected, selected, eta)
  } else {

    max.i <- 1
    max.j <- 1
    current.dis <- 0
    max.dis <- 0
    for(i in 1:length(search.grid)){
      lambda1 <- search.grid[i]
      for(j in 1:length(search.grid)){
        lambda2 <- search.grid[i]
        alpha_beta_pair <- sCCA_zero_padding(A_tuta, lambda1, lambda2, eta)
        if(any(alpha_beta_pair[[1]] != 0) & any(alpha_beta_pair[[2]] != 0)){
          current.dis <- distance_cal(alpha_beta_pair[[1]], alpha_beta_pair[[2]], cov.tuneX, cov.tuneY, cov.tuneXY)
          if(current.dis > max.dis){
            max.i <- i
            max.j <- j
          }
        }
      }
    }
    final_alpha_beta_pair <- sCCA_zero_padding(A_tuta, search.grid[max.i], search.grid[max.j], eta)
  }

  return(final_alpha_beta_pair)
}
