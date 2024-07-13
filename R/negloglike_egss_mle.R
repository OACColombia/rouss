#' Numerical optimization for MLE of the EGSS
#'
#'@description
#' Function needed to compute the Maximum Likelihood Estimation of parameters within the multivariate log-likelihood for the Exponential Growth State-Space (EGSS) population model, see Humbert _et al._ (2009).
#'
#' @param yt A vector of log-abundance observations in the time-series.
#' @param tt A vector of observations times.
#' @param fguess_egss A vector of first guess of the four parameters, from `guess_egss()`
#'
#' @return the value of negative log-Likelihood for the OUSS MLE
#' @export
#'
#' @examples
#' yt1 = log(c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6))
#' tt1 = c(1966:1995)
#' negloglike_egss_mle(yt = yt1, tt = tt1, fguess_egss = guess_egss(yt = yt1, tt = tt1))

negloglike_egss_mle <- function(yt,tt,fguess_egss){

  theta         <- fguess_egss[1];
  sigmasq       <- exp(fguess_egss[2]);
  tausq         <- exp(fguess_egss[3]);
  xo            <- fguess_egss[4];
  q             <- length(yt) - 1;
  qp1           <- q+1;
  yt            <- matrix(yt,nrow=qp1,ncol=1); # makes data a matrix object
  vx            <- matrix(0,qp1,qp1);
  for(i in 1:q){
    vx[((i+1):qp1),((i+1):qp1)] <- matrix(1,(qp1-i),(qp1-i))*tt[(i+1)];
  }
  Sigma.mat     <- sigmasq*vx;
  Itausq        <- matrix(rep(0,(qp1*qp1)), nrow=qp1, ncol=qp1);
  diag(Itausq)  <- rep(tausq,qp1);
  V             <- Sigma.mat + Itausq;
  theta.vec     <- matrix((xo+theta*tt), nrow=qp1,ncol=1);

  return((qp1/2)*log(2*pi) + 0.5*log(det(V)) + 0.5*(t(yt-theta.vec)%*%ginv(V)%*%(yt-theta.vec)))
}
