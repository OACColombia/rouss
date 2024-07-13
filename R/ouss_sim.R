#' Simulation of OUSS model
#'
#' @description
#' Function to simulate stationary stochastic population dynamics from parameters fitted under the Ornstein-Uhlenbeck State-Space (OUSS) population model.
#'
#' @param nsims The number of bootstrap replicates to simulate (≥2000)
#' @param tt The ORIGINAL vector of observation times (𝑡_0, 𝑡_1, 𝑡_2, …, 𝑡_𝑞)
#' @param parms A vector of parameters values estimated from `ouss_remle()` (it has better statistical properties than `ouss_mle()`)
#'
#' @return A matrix of size `nsims` with simulations of log abundance in each time step of the timeseries (using `randmvn()`)
#' @export
#'
#' @examples
#' yt1 = log(c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6))
#' tt1 = c(1966:1995)
#' B = 10 #it was recommended 2000 for estimate CI
#'
#' OUSS.REMLE.model = ouss_mle(yt = yt1, tt = tt1, fguess = guess_ouss(yt = yt1, tt = tt1))
#'
#' ouss_sim(nsims = B, tt = tt1, parms = OUSS.REMLE.model$remles)

ouss_sim <- function(nsims,tt,parms){

  # Time-vector starting in 0.
  t.i       <- tt-tt[1];
  # Number of time-series transitions
  q         <- length(t.i)-1;
  # length of time-series
  qp1       <- q+1;

  # parameters
  mu        <- parms[1];
  theta     <- parms[2];
  betasq    <- parms[3];
  tausq     <- parms[4];

  Var.inf   <- betasq/(2*theta);
  t.s       <- t.i[2:qp1] - t.i[1:q];
  t.cols    <- matrix(rep(t.i,each=qp1),
                      nrow=qp1,
                      ncol=qp1,
                      byrow=FALSE);
  t.rows    <- t(t.cols);
  abs.diffs <- abs(t.rows-t.cols);
  V         <- Var.inf*exp(-theta*abs.diffs);
  diag(V)   <- diag(V) + rep(tausq,qp1);
  m.vec     <- rep(mu,qp1);
  out       <- randmvn(n=nsims,
                       mu.vec=m.vec,
                       cov.mat = V)
  return(out)
}
