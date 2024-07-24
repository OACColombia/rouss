#' Simulation of EGSS model
#'
#' @description
#' Function to simulate exponential stochastic population dynamics from parameters fitted under the Exponential Growth State-Space models.
#'
#' @param nsims The number of bootstrap replicates to simulate (≥2000)
#' @param tt The ORIGINAL vector of observation times (𝑡_0, 𝑡_1, 𝑡_2, …, 𝑡_𝑞)
#' @param parms A vector of parameters values estimated from `egss_remle()`
#'
#' @return A matrix of size `nsims` with simulations of log abundance in each time step of the timeseries (using `randmvn()`)
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

  t.i           <- tt-tt[1];
  q             <- length(t.i)-1;
  qp1           <- q+1;
  t.s     <- t.i[2:qp1] - t.i[1:q];

  # parameters
  sigmasq       <- parms[1];
  tausq         <- parms[2];

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

  D1mat=cbind(-diag(1/t.s),
              matrix(0,q,1))+cbind(matrix(0,q,1),
                                   diag(1/t.s));

  V1mat=D1mat%*%V%*%t(D1mat);

  W.t=(yt[2:qp1]-yt[1:q])/t.s;
  j1=matrix(1,q,1);
  V1inv=ginv(V1mat);

  theta.remle=(t(j1)%*%V1inv%*%W.t)/(t(j1)%*%V1inv%*%j1);
  j=matrix(1,qp1,1);
  Vinv=ginv(V);
  x0.remle=(t(j)%*%Vinv%*%(yt-theta.remle*t.i))/(t(j)%*%Vinv%*%j);

  theta.vec     <- matrix((c(x0.remle)+c(theta.remle)*t.i),
                          nrow=qp1,
                          ncol=1);
  out           <- randmvn(n=nsims,
                           mu.vec=theta.vec,
                           cov.mat=V);

  return(out)
}
