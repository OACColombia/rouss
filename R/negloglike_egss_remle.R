#' Numerical optimization for REMLE of the EGSS
#'
#'@description
#' Function needed to compute the Restricted Maximum Likelihood Estimation of parameters within the multivariate log-likelihood for the Exponential Growth State-Space (EGSS) population model, see Humbert _et al._ (2009).
#'
#' @param fguess_egss A vector of first guess of the four parameters, from `guess_egss()`
#' @param yt A vector of log-abundance observations in the time-series.
#' @param tt A vector of observations times.
#'
#' @return the value of negative log-Likelihood for the EGSS REMLE
#' @export
#'
#' @examples
#' yt1 = log(c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6))
#' tt1 = c(1966:1995)
#' negloglike_egss_remle(yt = yt1, tt = tt1, fguess_egss = guess_egss(yt = yt1, tt = tt1))

negloglike_egss_remle <- function(fguess_egss,yt,tt){

  sigmasq     <- exp(fguess_egss[1]); #in egss_remle, this have only the two parameters
  tausq       <- exp(fguess_egss[2]);
  q           <- length(yt) - 1;
  qp1         <- q+1;

  ss          <- tt[2:qp1]-tt[1:q];
  wt          <- (yt[2:qp1]-yt[1:q])/ss;
  ut          <- wt[2:q]-wt[1:q-1];
  vx          <- matrix(0,qp1,qp1);
  for(i in 1:q){
    vx[(i+1):qp1,(i+1):qp1] <- matrix(1,(qp1-i),(qp1-i))*tt[i+1];
  }
  Sigma.mat   <- sigmasq*vx;
  Itausq      <- matrix(rep(0,(qp1*qp1)), nrow=qp1, ncol=qp1);
  diag(Itausq)<- rep(tausq,qp1);
  V           <- Sigma.mat + Itausq;

  D1mat       <- cbind(-diag(1/ss),matrix(0,q,1))+cbind(matrix(0,q,1),diag(1/ss));
  D2mat       <- cbind(-diag(1,(q-1)),matrix(0,(q-1),1)) + cbind(matrix(0,(q-1),1),diag(1,(q-1)));
  V2          <- D2mat%*%D1mat%*%V%*%t(D1mat)%*%t(D2mat);

  ofn=((q-1)/2)*log(2*pi)+(0.5*log(det(V2))) + (0.5*(ut%*%ginv(V2)%*%ut));

  return(ofn)
}
