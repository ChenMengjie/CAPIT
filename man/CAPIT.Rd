% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CAPIT.R
\name{CAPIT}
\alias{CAPIT}
\title{A function implements CAPIT algorithm for sparse CCA.}
\usage{
CAPIT(X, Y, cross.validation = TRUE, precision = FALSE, preX1 = NULL,
  preY1 = NULL, preX2 = NULL, preY2 = NULL, method = "thresholding",
  lambda = seq(0.001, 0.3, length.out = 20), selection = "hard",
  alpha = 0.3, select.covariance = "loglikelihood", search.grid = seq(0.5,
  3, by = 0.5), lambda1 = NULL, lambda2 = NULL, same.threshold = TRUE)
}
\arguments{
\item{X}{A $n \times p1$ data matrix.}

\item{Y}{A $n \times p2$ data matrix.}

\item{cross.validation}{Whether the threholding paramters for zero padding are selected by cross validation. The default is TRUE.}

\item{precision}{Whether to input estimated precision matrix. The default is FALSE.}

\item{preX1}{When ``precision=TURE", users can input precision matrix estimated from other methods. \code{preX1, preX2, preY1, preY2}  are used to input the pre-calculated precison matrix estimate for X1, X2, Y1, Y2, respectively,
where X1 and X2 are two halves of X, similarly Y1 and Y2 are two halves of Y.}

\item{method}{The method used to estimate sparse covariance matrix. The default is ``thresholding" algorithm, which further contains hard thresholding and soft thresholding.
Other choices include ``toeplitz" algorithm: implements the method proposed in Cai et al. (2013), assuming the Toeplitz structure is known.
``tapering"  algorithm: implements the tapering procedure proposed in Cai et al. (2010), assuming covariance decay as they move away from the diagonal.}

\item{lambda}{A sequence of thresholding levels will be tested in cross-validation in ``thresholding" algorithm.}

\item{selection}{The thresholding method used to estimate sparse covariance matrix. Choices include ``hard" and ``soft".}

\item{alpha}{Specified decay rate in ``tapering"  algorithm. The default is 0.3.}

\item{search.grid}{A sequence of thresholding levels will be tested in cross-validation for zero padding.}

\item{lambda1}{when ``cross.validation = FALSE'', the thresholding level used in the initialization is required as input.}

\item{lambda2}{when ``cross.validation = FALSE'', the thresholding level for iterative thresholding is required as input.}

\item{same.threshold}{When searching for the best performed thresholding levels through cross validation, whether to assume the thresholding levels used in the initialization and iterative thresholding are the same. The default is TRUE.}

\item{select.covariance:}{Criterion used in cross validation of sparse covariance matrix estimation. The default is ``loglikelihood".
When specified to "floss", frobenius loss between training and testing is minimized.}
}
\value{
The output of \code{CAPIT} is a list, which further contains two list \code{res} and \code{resOLS}.
\code{res} contains \code{alpha} and \code{beta}, which returns canonical vectors of X and Y estimated from vanilla CAPIT algorithm, respectively.
\code{res} contains \code{alpha} and \code{beta}, which returns canonical vectors of X and Y estimated from CAPIT algorithm with refinement using ordinary least squares, respectively.
}
\description{
Sparse Canonical Correlation Analysis (CCA), as an approach to study the maximal correlation between two sets of random variables,
has received considerable attention in high-dimensional data analysis. Despite of its popularity in application, there has been remarkably
few theoretical studies on sparse CCA. We propose Sparse CCA via Precision Adjusted Iterative Thresholding (CAPIT), a procedure that is rate-optimal
under various assumptions on nuisance parameters with theoretical guarantee.
Our algorithm has two steps, first step is to estimate covariance/precison matrix and the second step is iterative thresholding.
For the first step, we implement three approaches for sparse covariance matrix estimation, hard/soft thresholding, tapering and toeplitz.
We also allow direct input of precision matrix from other methods.
For the second step, we implements iterative thresholding algorithm, which selects thresholding levels through cross validation. Pre-specified thresholding levels are also allowed.
}
\examples{
### Generate simulate data
n <- 750*2
p1 <- 200
p2 <- 200
s1 <- 5
s2 <- 5
rho <- 0.3
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
theta <- theta/as.numeric(sqrt(t(theta)\%*\%Sigma1\%*\%theta))
eta <- as.matrix(c(rep(c(1, 0, 0, 0, 0), s2), rep(0, p2-5*s2)))
eta <- eta/as.numeric(sqrt(t(eta)\%*\%Sigma2\%*\%eta))
lambda <- 0.9
a1 <- cbind(Sigma1, lambda*Sigma1\%*\%theta\%*\%t(eta)\%*\%Sigma2)
a2 <- cbind(lambda*Sigma2\%*\%eta\%*\%t(theta)\%*\%Sigma1, Sigma2)
sigma_cov <- rbind(a1, a2)
library(MASS)
v <- c(theta, eta)
set.seed(20453)
Z <- mvrnorm(n, rep(0, p1+p2), sigma_cov)
X <- Z[, 1:p1]
Y <- Z[, (p1+1):(p1+p2)]
## u1.1 <- CAPIT(X, Y, method = "toeplitz")
## u1.2 <- CAPIT(X, Y, method = "tapering")
u1.3 <- CAPIT(X, Y)
}

