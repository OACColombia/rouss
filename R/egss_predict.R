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
#' parms1 = egss_mle(yt = yt1, tt = tt1, fguess_egss = guess_egss(yt = yt1, tt = tt1))
#'
#' egss_predict(yt = yt1, tt = tt1, parms = parms1$mles)
#'
egss_predict <- function(yt,tt,parms, plot.it="TRUE"){

  # Time-vector starting in 0.
  t.i     <- tt-tt[1];
  q       <- length(t.i)-1;
  qp1     <- q+1;
  t.s     <- t.i[2:qp1] - t.i[1:q];

  # parameters (only sigmasqr and tausqr - REMLE)
  sigmasq <- parms[1];
  tausq   <- parms[2];

  # Missing t.s
  nmiss   <- t.s-1;
  long.nmiss <- c(0,nmiss);
  Nmiss   <- sum(nmiss)

  vx      <- matrix(0,qp1,qp1);
  for(i in 1:q){
    vx[((i+1):qp1),((i+1):qp1)] <- matrix(1,(qp1-i),(qp1-i))*t.i[(i+1)];
  }
  Sigma.mat    <- sigmasq*vx;
  Itausq       <- matrix(rep(0,(qp1*qp1)),
                         nrow=qp1,
                         ncol=qp1);
  diag(Itausq) <- rep(tausq,qp1);
  V            <- Sigma.mat + Itausq;

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
  x0.remle=(t(j)%*%Vinv%*%(yt-c(theta.remle)*t.i))/(t(j)%*%Vinv%*%j);
  Var_theta.remle=1/(t(j1)%*%V1inv%*%j1); # Variance of theta
  theta_hi.remle=theta.remle+1.96*sqrt(Var_theta.remle); # 95% CI for theta
  theta_lo.remle=theta.remle-1.96*sqrt(Var_theta.remle)

  #Calculate estimated population size for EGSS model

  m=rep(1,qp1); # Will contain Kalman means for Kalman calculations.
  v=rep(1,qp1); # Will contain variances for Kalman calculations.

  m[1]=x0.remle; # Initial mean of Y(t).
  v[1]=tausq; # Initial variance of Y(t).

  for (ti in 1:q) # Loop to generate estimated population abundances
  { # using Kalman filter (see equations 6 & 7, # Dennis et al. (2006)).
    m[ti+1]=theta+(m[ti]+((v[ti]-tausq)/v[ti])*(yt[ti]-m[ti]));
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
