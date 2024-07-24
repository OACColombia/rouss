#' Compute the REMLE of the EGSS model
#'
#'@description
#' This function computes the Restricted Maximum Likelihood Estimates for the Exponential Growth State-Space (EGSS) population model.
#'
#' @param yt A vector of log-abundance observations in the time-series.
#' @param tt A vector of observations times.
#' @param fguess_egss A vector of first guess of the four parameters, from `guess_egss()`
#'
#' @return A list of REMLEs of two parameters (sigmasqr and tausqr), estimated log-likelihood, and Akaike Information Criteria.
#' @export
#'
#' @examples
#' yt1 = log(c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6))
#' tt1 = c(1966:1995)
#' egss_remle(yt = yt1, tt = tt1, fguess_egss = guess_egss(yt = yt1, tt = tt1))

egss_remle <- function(yt,tt,fguess_egss){

  #Temporal vectors
  t.i         <- tt-tt[1];
  q           <- length(t.i)-1;
  qp1         <- q+1;
  t.s     <- t.i[2:qp1] - t.i[1:q];

  # initial guesses (sigmasq and tausq at log scale)
  guess.optim <- c(log(fguess_egss[2:3]))
  # numerical optimization
  optim.out   <- optim(par=guess.optim,
                       fn=negloglike_egss_remle,
                       method="Nelder-Mead",
                       yt=yt,
                       tt=t.i)

  #extract parameters estimated by REML
  sigmasq <- exp(optim.out$par)[1]
  tausq <- exp(optim.out$par)[1]

  #to estimate trend parameter (theta) and initial population (x0)
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

  #Trend parameter
  theta.remle=(t(j1)%*%V1inv%*%W.t)/(t(j1)%*%V1inv%*%j1);

  j=matrix(1,qp1,1);
  Vinv=ginv(V);

  #initial population
  x0.remle=(t(j)%*%Vinv%*%(yt-as.numeric(theta.remle)*t.i))/(t(j)%*%Vinv%*%j);

  #Extract REMLEs and AIC
  remles        <- c(theta.remle,exp(optim.out$par[1:2]),x0.remle)
  lnL.hat     <- - optim.out$value[1]
  AIC         <- -2*lnL.hat + 2*2 #where 2 = length(REMLEs)...

  out         <- list(remles=remles,
                      lnL.hat = lnL.hat,
                      AIC=AIC)
  return(out)
}
