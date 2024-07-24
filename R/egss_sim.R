#' Simulation of EGSS model
#'
#' @description
#' Function to simulate exponential stochastic population dynamics from parameters fitted under the Exponential Growth State-Space models.
#'
#' @param nsims The number of bootstrap replicates to simulate (≥2000)
#' @param tt The ORIGINAL vector of observation times (𝑡_0, 𝑡_1, 𝑡_2, …, 𝑡_𝑞)
#' @param parms A vector of parameters values estimated from `egss_remle()` (recommended, althoug also works with `egss_mle()`)
#'
#' @return A matrix of size `nsims` with simulations of log abundance in each time step of the time-series (using `randmvn()`)
#' @export
#'
#' @examples
#' yt1 = log(c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6))
#' tt1 = c(1966:1995)
#' B = 10 #it was recommended 2000 for estimate CI
#'
#' EGSS.REMLE.model = egss_remle(yt = yt1, tt = tt1, fguess_egss = guess_egss(yt = yt1, tt = tt1))
#'
#' egss_sim(nsims = B, tt = tt1, parms = EGSS.REMLE.model$remles)

egss_sim <- function(nsims,tt,parms){

  # time and temporal scale
  t.i           <- tt-tt[1];
  q             <- length(t.i)-1;
  qp1           <- q+1;

  # parameters
  theta         <- parms[1];
  sigmasq       <- parms[2];
  tausq         <- parms[3];
  x0            <- parms[4];

  vx      <- matrix(0,qp1,qp1);
  for(i in 1:q){
    vx[((i+1):qp1),((i+1):qp1)] <- matrix(1,(qp1-i),(qp1-i))*t.i[(i+1)];
  }

  Sigma.mat     <- sigmasq*vx;
  Itausq        <- matrix(rep(0,(qp1*qp1)),
                          nrow=qp1,
                          ncol=qp1);
  diag(Itausq)  <- rep(tausq,qp1);
  V             <- Sigma.mat + Itausq;
  theta.vec     <- matrix((x0+theta*t.i),
                          nrow=qp1,
                          ncol=1);
  out           <- randmvn(n=nsims,
                           mu.vec=theta.vec,
                           cov.mat=V);

  return(out)
}
