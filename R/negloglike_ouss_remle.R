#' Numerical optimization for REMLE of the OUSS
#'
#'@description
#' Function needed to compute the Restricted Maximum Likelihood Estimation of parameters within the multivariate log-likelihood for the stationary Ornstein-Uhlenbeck State-Space (OUSS) population model, see Dennis & Ponciano (2014).
#'
#' @param yt A vector of log-abundance observations in the time-series.
#' @param tt A vector of observations times.
#' @param fguess A vector of first guess of the four parameters, from `guess_ouss()` (this function will use only 2:4 parameters)
#'
#' @return the value of negative log-Likelihood for the OUSS REMLE
#' @export
#'
#' @examples
#' yt1 = log(c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6))
#' tt1 = c(1966:1995)
#' negloglike_ouss_remle(yt = yt1, tt = tt1, fguess = guess_ouss(yt = yt1, tt = tt1))

negloglike_ouss_remle=function(yt,tt,fguess){
  # Constrains parameters theta, beta^2, and tau^2 > 0

  # speed of equilibration (Eq1 in DP_E)
  theta         <- exp(fguess[2]);
  # variability of process noise
  betasq        <- exp(fguess[3]);
  # variability of sampling
  tausq         <- exp(fguess[4]);
  # number of time-series transitions
  q             <- length(yt) - 1;
  # length of time-series
  qp1           <- q+1;
  # Variance (Eq11 in DP_E)
  Var.inf       <- betasq/(2*theta);
  # time intervals (not used here?)
  t.s           <- tt[2:qp1] - tt[1:q];
  # part of Eq18 in DP_E
  t.cols        <- matrix(rep(tt,each=qp1),
                          nrow=qp1,
                          ncol=qp1,
                          byrow=FALSE);
  # (part of Eq18 in DP_E)
  t.rows        <- t(t.cols);
  # (part of Eq18 in DP_E)
  abs.diffs     <- abs(t.rows-t.cols);

  # Covariance of the process (Eq18 in DP_E)
  Sigma.mat     <- Var.inf*exp(-theta*abs.diffs);
  # Create a matrix full of 0s of the length of time series
  Itausq        <- matrix(0,qp1,qp1);
  # Repeat the observation error variance guess in the diagonal of the matrix
  diag(Itausq)  <- rep(tausq,qp1);
  # add Covariance with the matrix
  V             <- Sigma.mat+Itausq;
  # Create the differencing matrix **D**
  Dmat          <- cbind(-diag(1,q),matrix(0,q,1)) + cbind(matrix(0,q,1),diag(1,q));
  # Variance-covariance matrix **Phi** (Eq20 DP_E)
  Phi.mat       <- Dmat%*%V%*%t(Dmat);
  # simple differencing of the observations (W_i? )
  wt            <- yt[2:qp1]-yt[1:q];

  # note the signs change because we want here the negative log-likelihood (Eq22*-1)
  neglogl         <- (q/2)*log(2*pi) + (1/2)*log(det(Phi.mat)) + (1/2)*wt%*%ginv(Phi.mat)%*%wt;

  # What to do if the `neglogl` is not finite? assign a big number of 50000
  if(is.infinite(neglogl)==TRUE){
    return(50000)}else{
      return(neglogl)}
}
