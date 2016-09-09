FLoss <- function(u, v, p1, root=TRUE){
  p <- length(u)
  u1 <- u[1:p1]
  u2 <- u[(p1+1):p]
  v1 <- v[1:p1]
  v2 <- v[(p1+1):p]
  U1 <- sum(u1^2)
  U2 <- sum(u2^2)
  V1 <- sum(v1^2)
  V2 <- sum(v2^2)
  mat1 <- u1%*%t(u1)/U1 - v1%*%t(v1)/V1
  mat2 <- u2%*%t(u2)/U2 - v2%*%t(v2)/V2
  loss <- max(norm(mat1, "F")^2, norm(mat2, "F")^2)
  return(sqrt(loss))
}
