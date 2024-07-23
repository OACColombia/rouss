egss_remle <- function(yt,tt,fguess_egss){

  # Time-vector starting in 0.
  t.i         <- tt-tt[1];
  # Number of time-series transitions
  q           <- length(t.i)-1;
  # length of time-series
  qp1         <- q+1;
  # initial guesses (r as estimated, sigmasq and tausq at log scale (?), x0 as estimated)
  guess.optim <- c(log(fguess_egss[2:3]))
  # numerical optimization
  optim.out   <- optim(par=guess.optim,
                       fn=negloglike_egss_remle,
                       method="Nelder-Mead",
                       yt=yt,
                       tt=t.i)
  #Extract MLEs and AIC
  remles        <- c(exp(optim.out$par[2:3]))
  lnL.hat     <- - optim.out$value[1]
  AIC         <- -2*lnL.hat + 2*2 #where 2 = length(REMLEs)...

  out         <- list(remles=remles,
                      lnL.hat = lnL.hat,
                      AIC=AIC)
  return(out)
}
