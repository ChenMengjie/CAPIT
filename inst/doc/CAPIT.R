### R code from vignette source 'CAPIT.Rnw'

###################################################
### code chunk number 1: CAPIT.Rnw:56-57
###################################################
options(width=60)


###################################################
### code chunk number 2: CAPIT.Rnw:99-135
###################################################
n <- 750*2
p1 <- 200
p2 <- 200

s1 <- 5
s2 <- 5
rho <- 0.2
Sigma1 <- matrix(0, ncol = p1, nrow = p1)
for(i in 1:p1){
	for(j in 1:p1){
		Sigma1[i, j] <- rho^(abs(i - j))
	}
}

Sigma2 <- matrix(0, ncol = p2, nrow = p2)
for(i in 1:p2){
	for(j in 1:p2){
		Sigma2[i, j] <- rho^(abs(i - j))
	}
}

theta <- as.matrix(c(rep(c(1, 0, 0, 0, 0), s1), rep(0, p1-5*s1)))
theta <- theta/as.numeric(sqrt(t(theta)%*%Sigma1%*%theta))
eta <- as.matrix(c(rep(c(1, 0, 0, 0, 0), s2), rep(0, p2-5*s2)))
eta <- eta/as.numeric(sqrt(t(eta)%*%Sigma2%*%eta))

lambda <- 0.9
sigma_cov <- rbind(cbind(Sigma1, lambda*Sigma1%*%theta%*%t(eta)%*%Sigma2),
                   cbind(lambda*Sigma2%*%eta%*%t(theta)%*%Sigma1, Sigma2))
require(MASS)

set.seed(100)
Z <- mvrnorm(n, rep(0, p1+p2), sigma_cov)
X <- Z[, 1:p1]
Y <- Z[, (p1+1):(p1+p2)]
v <- c(theta, eta)


###################################################
### code chunk number 3: CAPIT.Rnw:143-145
###################################################
library(CAPIT)
u1 <- CAPIT(X, Y)


###################################################
### code chunk number 4: CAPIT.Rnw:151-153
###################################################
FLoss(unlist(u1[[1]]), v, p1) #CAPIT
FLoss(unlist(u1[[2]]), v, p1) #CAPIT + refinement by OLS


###################################################
### code chunk number 5: CAPIT.Rnw:158-162
###################################################
u2 <- CAPIT(X, Y, method = "toeplitz")
FLoss(unlist(u2[[1]]), v, p1) #CAPIT
FLoss(unlist(u2[[2]]), v, p1) #CAPIT + refinement by OLS



###################################################
### code chunk number 6: CAPIT.Rnw:165-168
###################################################
u3 <- CAPIT(X, Y, method = "tapering")
FLoss(unlist(u3[[1]]), v, p1) #CAPIT
FLoss(unlist(u3[[2]]), v, p1) #CAPIT + refinement by OLS


