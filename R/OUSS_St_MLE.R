#' Compute the MLE of the OUSS model
#'
#'@description
#' This function comput the Maximum Likelihood Estimates for the stationary Ornstein-Uhlenbeck State-Space (OUSS) population model.
#'
#' @param yt A vector of log-abundance observations in the time-series.
#' @param tt A vector of observations times.
#' @param fguess A vector of first guess of the four parameters, from `guess_ouss()`
#'
#' @return list of MLEs of four parameters, estimated log-likelihood, and Akaike Information Criteria.
#' @export
#'
#' @examples
#' yt1 = log(c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6))
#' tt1 = c(1966:1995)
#' OUSS_St_MLE(yt = yt1, tt = tt1, fguess = guess_ouss(yt = yt1, tt = tt1))

OUSS_St_MLE <- function(yt,tt,fguess){

  # Time-vector starting in 0.
  t.i         <- tt-tt[1];
  # Number of time-series transitions
  q           <- length(yt)-1;
  # length of time-series
  qp1         <- q+1;
  # initial guesses (mu as estimated, theta, betasq, and tausq at log scale (?))
  guess.optim <- c(fguess[1], log(fguess[2:4]))
  # numerical optimization
  optim.out   <- MASS::optim(par=guess.optim,
                       fn=negloglike_OU_mle,
                       method="Nelder-Mead",
                       yt=yt,
                       tt=t.i)
  #Extract MLEs and AIC
  mles        <- c(optim.out$par[1],
                   exp(optim.out$par[2:4]))
  lnL.hat     <- - optim.out$value[1]
  AIC         <- -2*lnL.hat + 2*4 #where 4 = length(MLEs)...

  out           <- list(mles=mles,
                        lnL.hat = lnL.hat,
                        AIC=AIC)
  return(out)
}
