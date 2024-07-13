#' Multivariate normal random number generator
#'
#' @description
#' This function generates random numbers from a multivariate normal distribution, required for parametric bootstrapping (prediction and simulation of state space models).
#'
#'
#' @param n the number of random samples of a multivariate normal vector
#' @param mu.vec mean vector of the multivariate normal distribution to sample from
#' @param cov.mat Variance-covariance matrix of the multivariate normal distribution to sample from
#'
#' @return A matrix from a Cholesky decomposition (factorization of a real symmetric positive-definite sqr matriz)
#' @export
#'
#' @examples
#' See `ouss_sim()`
#'
randmvn <- function(n, mu.vec, cov.mat){

  # Save the length of the mean vector of the multivariate normal distribution to sample
  p         <- length(mu.vec);
  # The Cholesky decomposition (factorization of a real symmetric positive-definite sqr matriz)
  Tau       <- chol(cov.mat, pivot=TRUE);
  # generate normal deviates outside loop
  Zmat      <- matrix(rnorm(n=p*n,mean=0,sd=1),nrow=p,ncol=n);

  # empty matrix
  out       <- matrix(0,nrow=p,ncol=n);
  # iterate
  for(i in 1:n){
    Z       <- Zmat[,i];
    out[,i] <- t(Tau)%*%Z + mu.vec
  }

  return(out)
}
