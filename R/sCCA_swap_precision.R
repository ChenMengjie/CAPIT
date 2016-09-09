sCCA_swap_precision <- function(X, Y, lambda1, lambda2, preX1, preY1, preX2, preY2){

  svd.regression <- function(u, s, X, Y){
    uID <- which(u != 0)
    aa <- u[uID]%*%t(X[, uID])
    sID <- which(s != 0)
    kk <- solve(t(Y[, sID])%*%Y[, sID])%*%t(Y[, sID])%*%t(aa)
    bb <- t(kk)%*%t(Y[, sID])
    mm <- solve(t(X[, uID])%*%X[, uID])%*%t(X[, uID])%*%t(bb)
    u[uID] <- mm
    s[sID] <- kk
    return(c(u, s))
  }

  ####### divide data into two parts ########
  if(nrow(X) != nrow(Y))
    stop("X and Y should have same number of rows!")
  n <- nrow(X)
  X1 <- X[1:floor(n/2), ]
  X2 <- X[floor(n/2+1):n, ]
  Y1 <- Y[1:floor(n/2), ]
  Y2 <- Y[floor(n/2+1):n, ]

  XY1 <- sCCA_precision(X1, Y1, X2, Y2, lambda1 = lambda1, lambda2 = lambda2, preX1, preY1)
  XY2 <- sCCA_precision(X2, Y2, X1, Y1, lambda1 = lambda1, lambda2 = lambda2, preX2, preY2)

  res <- cbind(XY1$alpha, XY2$alpha)
  ll <- apply(res, 2, function(x){
    if(length(x[x>0]) < length(x[x<0]))
      return(-x)
    else return(x)
  })
  alpha <- apply(ll, 1, mean)

  res <- cbind(XY1$beta, XY2$beta)
  ll <- apply(res, 2, function(x){
    if(length(x[x>0]) < length(x[x<0]))
      return(-x)
    else return(x)
  })
  beta <- apply(ll, 1, mean)
  res <- list(Theta = alpha, Eta = beta)

  after.regression <- svd.regression(alpha, beta, X, Y)
  p1 <- ncol(X)
  p2 <- ncol(Y)
  resOLS <- list(ThetaOLS = after.regression[1:p1], EtaOLS = after.regression[(p1+1):(p1+p2)])
  result <- list(res = res, resOLS = resOLS)

  return(result)
}
