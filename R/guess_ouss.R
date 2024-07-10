#' Generate a vector of parameter guess for OUSS
#'
#'@description
#' The numerical optimization (Maximum Likelihood Estimate - MLE and Restricted Maximum Likelihood Estimates - REMLE) for the Ornstein-Uhlenbeck State-Space (OUSS) population model requires a first guess of the four parameters:
#' * 𝜇 (in stationary distribution is the expected value of stationary distribution of log-abundance mean; it is different than 𝜇91 in Dennis _et al._ (1991) and 𝜇 in Humbert _et al._ (2009), named here 𝜃 in EGSS),
#' * 𝜃 (rate to approach to stationarity),
#' * 𝛽^2 (variability of process noise),
#' * 𝜏^2 (variability of sampling).
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
#' guess_ouss(yt = yt1, tt = tt1)
guess_ouss <- function(yt,tt){

  # Time-vector starting in 0.
  t.i     <- tt-tt[1];
  # Number of time-series transitions
  q       <- length(yt)-1;
  # length of time-series
  qp1     <- q+1;
  # time intervals
  t.s     <- t.i[2:qp1]-t.i[1:q];
  # mean of the observations as assumed to arise from stationary distribution
  Ybar    <- mean(yt);
  # Variance of the observations
  Yvar    <- sum((yt-Ybar)*(yt-Ybar))/q;
  # Initial mu estimate (at stationary distribution)
  mu1     <- Ybar;

  # Kludge an initial value for theta based on mean of Y(t+s) given Y(t).
  th1     <- -mean(log(abs((yt[2:qp1]-mu1)/(yt[1:q]-mu1)))/t.s);
  # Moment estimate using stationary distribution
  bsq1    <- 2*th1*Yvar/(1+2*th1);
  # Observation error variance, assumed as first guess as betasq=tausq.
  tsq1    <- bsq1;

  # What to do if initial guesses is three 0's (or NAs)? Assume arbitrary values
  three0s <- sum(c(th1,bsq1,tsq1))

  if(three0s==0|is.na(three0s)){
    th1   <- 0.5;
    bsq1  <- 0.09;
    tsq1  <- 0.23;}

  out1    <- c(th1,bsq1,tsq1);

  # What to do if initial guesses are too little? Assume arbitrary values
  if(sum(out1<1e-7)>=1){
    out1  <- c(0.5,0.09,0.23)}

  out     <- c(mu1,out1);

  return(abs(out))
}
