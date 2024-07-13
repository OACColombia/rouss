#' Generate a vector of parameter guess for EGSS
#'
#'@description
#' The numerical optimization (Maximum Likelihood Estimate - MLE) for the Exponential Growth State-Space (EGSS) population model requires a first guess of the four parameters:
#' * 𝜃 (𝑟−1/2𝛽^2; or population growth rate, as defined by 𝜇91 in Dennis _et al._ (1991) and 𝜇 in Humbert _et al._ (2009)),
#' * 𝜎^2 (variability of process noise),
#' * 𝜏^2 (variability of sampling).
#' * 𝑥_0 (initial population in log-scale)
#'
#' This function provides the first guess of these parameters, roughly computed by provide the vector of log-abundance observations `yt` (𝐲) and the vector of observation times `tt`
#'
#' @param yt A vector of log-abundance observations in the time-series.
#' @param tt A vector of observations times.
#'
#' @return A vector of four values, the first guess of the parameters to conduct numerical optimization.
#' @export
#'
#' @examples
#' yt1 = log(c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6))
#' tt1 = c(1966:1995)
#' guess_egss(yt = yt1, tt = tt1)
#'
guess_egss <- function(yt,tt){

  # Time-vector starting in 0.
  t.i       <- tt-tt[1];
  # Number of time-series transitions
  q         <- length(yt)-1;
  # length of time-series
  qp1       <- q+1;
  # time intervals (named as S.t in H_O and sometimes in DP_E)
  t.s       <- t.i[2:qp1]-t.i[1:q];

  #The Exponential Growth Observation Error (EGOE in H_O) initial values

  # mean of the observations as assumed to arise from stationary distribution
  Ybar      <- mean(yt);
  # mean of the time series
  Tbar      <- mean(t.i)
  # trend parameter for EGOE (theta = ln(lambda) in H_O)
  theta.egoe    <- sum((t.i-Tbar)*(yt-Ybar))/sum((t.i-Tbar)*(t.i-Tbar));
  # Initial population  of EGOE
  x0.egoe   <- Ybar-theta.egoe*Tbar
  # sigma square for EGOE is 0 (assume no ecological process variation)
  ssq.egoe  <- 0
  # estimate of initial population observed under EGOE
  Yhat.egoe <- x0.egoe+theta.egoe*t.i;
  # initial value for tau^2
  tsq.egoe  <- sum((yt-Yhat.egoe)*(yt-Yhat.egoe))/(q-1);

  #The Exponential Growth Process Noise (EGPN in H_O) initial values
  # Square root of time intervals (time trend?)
  Ttr       <- sqrt(t.s);
  # Observed trend?
  Ytr       <- (yt[2:qp1]- yt[1:q])/Ttr;
  # trend parameter for EGPN (mu = ln(lambda) in H_O)
  theta.egpn    <- sum(Ttr*Ytr)/sum(Ttr*Ttr);
  # Trend of observed estimated
  Ytrhat    <- theta.egpn*Ttr;
  # initial value for sigma^2
  ssq.egpn  <- sum((Ytr-Ytrhat)*(Ytr-Ytrhat))/(q-1);
  # tau square for EGPN is 0 (assume no observation variation)
  tsq.egpn  <- 0;
  # Initial population  of EGPN is the first observation
  x0.egpn   <- yt[1];

  #four parameters needed in EGSS and OUSS.NoSt
  theta0    <- (theta.egoe+theta.egpn)/2;
  ssq0      <- ssq.egpn/2;
  tsq0      <- tsq.egoe/2;
  x0.out    <- (x0.egoe+x0.egpn)/2;

  return(c(theta0, ssq0, tsq0, x0.out))
}
