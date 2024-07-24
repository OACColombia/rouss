#' Predict trajectory for EGSS
#'
#' @description
#' UNFINISHED function that provides prediction trajectory of a population under EGSS model. It is still in development.
#'
#' @param yt A vector of log-abundance observations in the time-series.
#' @param tt A vector of observations times.
#' @param parms Parameters estimates of REMLEs from `egss_remle()
#'
#' @return A dataframe of the original tt vector, Prediction (REMLE), and observed vector (real)
#' @export
#'
#' @examples
#' yt1 = log(c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6))
#' tt1 = c(1966:1995)
#' parms1 = egss_remle(yt = yt1, tt = tt1, fguess_egss = guess_egss(yt = yt1, tt = tt1))
#'
#' egss_predict(yt = yt1, tt = tt1, parms = parms1$remles)
#'
egss_predict <- function(yt,tt,parms, plot.it="TRUE"){

  # Time-vector starting in 0.
  t.i     <- tt-tt[1];
  q       <- length(t.i)-1;
  qp1     <- q+1;
  t.s     <- t.i[2:qp1] - t.i[1:q];

  # parameters ()
  theta.remle <- parms[1];
  sigmasq     <- parms[2];
  tausq       <- parms[3];
  x0.remle    <- parms[4];

  #Calculate estimated population size for EGSS model

  m=rep(1,qp1); # Will contain Kalman means for Kalman calculations.
  v=rep(1,qp1); # Will contain variances for Kalman calculations.

  m[1]=x0.remle; # Initial mean of Y(t).
  v[1]=tausq; # Initial variance of Y(t).

  for (ti in 1:q) # Loop to generate estimated population abundances
  { # using Kalman filter (see equations 6 & 7, # Dennis et al. (2006)).
    m[ti+1]=theta.remle+(m[ti]+((v[ti]-tausq)/v[ti])*(yt[ti]-m[ti]));
    v[ti+1]=tausq*((v[ti]-tausq)/v[ti])+sigmasq+tausq;
  }

  # The following statement calculates exp{E[X(t) | Y(t), Y(t-1),...,Y(0)]};
  # see equation 54 in Dennis et al. (2006).

  Predict.t=exp(m+((v-tausq)/v)*(yt-m));

  if(plot.it=="TRUE"){
    #  Plot the data & model-fitted values
    #X11()
    plot(tt,exp(yt),xlab="Time",ylab="Population abundance",type="b",cex=1.5,
         main="Predicted (--) and observed (-o-) abundances");#  Population data are circles.
    par(lty="dashed"); #  Predicted abundances are dashed line.
    points(tt,Predict.t, type="l", lwd=1);
  }

  return(data.frame(Time.t = tt, REMLE = Predict.t, Observed.y = exp(yt)))
}
