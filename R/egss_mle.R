#' Compute the MLE of the EGSS model
#'
#'@description
#' This function computes the Maximum Likelihood Estimates for the Exponential Growth State-Space (EGSS) population model.
#'
#' @param yt A vector of log-abundance observations in the time-series.
#' @param tt A vector of observations times.
#' @param fguess_egss A vector of first guess of the four parameters, from `guess_egss()`
#'
#' @return A list of MLEs of four parameters, estimated log-likelihood, and Akaike Information Criteria.
#' @export
#'
#' @examples
#' yt1 = log(c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6))
#' tt1 = c(1966:1995)
#' egss_mle(yt = yt1, tt = tt1, fguess_egss = guess_egss(yt = yt1, tt = tt1))

egss_mle <- function(yt,tt,fguess_egss){

  # Time-vector starting in 0.
  t.i         <- tt-tt[1];
  # Number of time-series transitions
  q           <- length(t.i)-1;
  # length of time-series
  qp1         <- q+1;
  # initial guesses (r as estimated, sigmasq and tausq at log scale (?), x0 as estimated)
  guess.optim <- c(fguess_egss[1],
                   log(fguess_egss[2:3]),
                   fguess_egss[4])
  # numerical optimization
  optim.out   <- optim(par=guess.optim,
                       fn=negloglike_egss_mle,
                       method="Nelder-Mead",
                       yt=yt,
                       tt=t.i)
  #Extract MLEs and AIC
  mles        <- c(optim.out$par[1],
                   exp(optim.out$par[2:3]),
                   optim.out$par[4])
  lnL.hat     <- - optim.out$value[1]
  AIC         <- -2*lnL.hat + 2*4 #where 4 = length(MLEs)...

  out         <- list(mles=mles,
                      lnL.hat = lnL.hat,
                      AIC=AIC)
  return(out)
}
