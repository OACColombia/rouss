#' Numerical optimization for MLE of the OUSS
#'
#'@description
#' Function needed to compute the Maximum Likelihood Estimation of parameters within the multivariate log-likelihood for the stationary Ornstein-Uhlenbeck State-Space (OUSS) population model, see Dennis & Ponciano (2014).
#'
#' @param yt A vector of log-abundance observations in the time-series.
#' @param tt A vector of observations times.
#' @param fguess A vector of first guess of the four parameters, from `guess_ouss()`
#'
#' @return the value of negative log-Likelihood for the OUSS MLE
#' @export
#'
#' @examples
#' yt1 = log(c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6))
#' tt1 = c(1966:1995)
#' negloglike_ouss_mle(yt = yt1, tt = tt1, fguess = guess_ouss(yt = yt1, tt = tt1))

negloglike_ouss_mle <- function(yt,tt,fguess){

  # log-abundance stationary distribution mean (Eq10 in Dennis & Ponciano 2014)
  mu        <- fguess[1];
  # Constrains parameters theta, beta^2, and tau^2 > 0
  guess     <- exp(fguess[2:4]);
  # speed of equilibration (Eq1 in Dennis & Ponciano 2014)
  theta     <- guess[1];
  # variability of process noise
  betasq    <- guess[2];
  # variability of sampling
  tausq     <- guess[3];
  # number of time-series transitions
  q         <- length(yt) - 1;
  # length of time-series
  qp1       <- q+1;
  # Variance (Eq11 in Dennis & Ponciano 2014)
  Var.inf   <- betasq/(2*theta);
  # time intervals (not used here?)
  t.s       <- tt[2:qp1] - tt[1:q];
  # part of Eq18 in Dennis & Ponciano 2014
  t.cols    <- matrix(rep(tt,each=qp1),
                      nrow=qp1,
                      ncol=qp1,
                      byrow=FALSE);
  # (part of Eq18 in Dennis & Ponciano 2014)
  t.rows    <- t(t.cols);
  # (part of Eq18 in Dennis & Ponciano 2014)
  abs.diffs <- abs(t.rows-t.cols);
  # Covariance (Eq18 in Dennis & Ponciano 2014)
  V         <- Var.inf*exp(-theta*abs.diffs);
  diag(V)   <- diag(V) + rep(tausq,qp1);
  # column vector **m** (from Eq16 in Dennis & Ponciano 2014)
  mu.vec    <- rep(mu,qp1);

  #note the signs change because we want here the negative log-likelihood (Eq19 in Dennis & Ponciano 2014)
  neglogl   <- (qp1/2)*log(2*pi) + (1/2)*log(det(V)) + (1/2)*(yt-mu.vec)%*%MASS::ginv(V)%*%(yt-mu.vec);

  return(neglogl)
}
