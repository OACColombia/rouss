#' Compute the REMLE of the OUSS model
#'
#'@description
#' This function computes the Restricted Maximum Likelihood Estimates for the stationary Ornstein-Uhlenbeck State-Space (OUSS) population model.
#'
#' @param yt A vector of log-abundance observations in the time-series.
#' @param tt A vector of observations times.
#' @param fguess A vector of first guess of the four parameters, from `guess_ouss()` (this function will use only 2:4 parameters)
#'
#' @return A list of MLEs of four parameters, estimated log-likelihood, and Akaike Information Criteria.
#' @export
#'
#' @examples
#' yt1 = log(c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6))
#' tt1 = c(1966:1995)
#' ouss_remle(yt = yt1, tt = tt1, fguess = guess_ouss(yt = yt1, tt = tt1))

ouss_remle <- function(yt, tt, fguess){

  # Time-vector starting in 0.
  t.i           <- tt-tt[1];
  # Number of time-series transitions
  # length of time-series
  q             <- length(yt)-1;
  qp1           <- q+1;
  # time intervals
  t.s           <- t.i[2:qp1]-t.i[1:q];
  # initial guesses (all, but negloglike.OU.remle will use only fguess[2:4])
  guess.optim   <- c(fguess[1],
                     log(fguess[2:4]));
  # numerical optimization
  optim.out     <- optim(par = guess.optim,
                         fn=negloglike_ouss_remle,
                         method="Nelder-Mead",
                         yt=yt,
                         tt=t.i);
  # Restricted maximum likelihood estimates (REMLE) and lnL.hat
  remles        <- exp(optim.out$par);
  theta.remle   <- remles[2];
  betasq.remle  <- remles[3];
  tausq.remle   <- remles[4];

  lnL.hat       <- -optim.out$value[1];

  # Variance (Eq11 in DP_E)
  Var.inf       <- betasq.remle/(2*theta.remle)
  # creates an matrix full of 1 dim qp1 x qp1
  vx            <- matrix(1,qp1,qp1);
  # iterate to fill the matrix (couldn't find vx in DP_E!)
  for (t.i in 1:q){
    vx[(t.i+1):qp1,t.i]=exp(-theta.remle*cumsum(t.s[t.i:q]));
    vx[t.i,(t.i+1):qp1]=vx[(t.i+1):qp1,t.i];
  }
  # ?
  Sigma.mat     <- vx*Var.inf;
  # Create a matrix full of 0s of the length of time series
  Itausq        <- matrix(0,qp1,qp1);
  # Repeat the observation error variance remle in the diagonal of the matrix
  diag(Itausq)  <- rep(tausq.remle,qp1);
  # Variance-covariance matrix (V.hat) evaluated with remles to estimate mu.hat
  V.remle       <- Sigma.mat+Itausq;
  # column vector matrix of ones
  j             <- matrix(1,qp1,1);
  # Inverse matrix (part of Eq23 in DP_E)
  Vinv          <- ginv(V.remle);
  # REMLE of mu (mu.hat) with Eq23 in DP_E
  mu.remle      <- (t(j)%*%Vinv%*%yt)/(t(j)%*%Vinv%*%j);
  #AIC
  AIC           <- -2*lnL.hat + 2*4 #where 4 = length(mles)...

  #Results
  out           <- list(remles = c(mu.remle,
                                   theta.remle,
                                   betasq.remle,
                                   tausq.remle),
                        lnLhat = lnL.hat,
                        AIC = AIC)
  return(out)
}
