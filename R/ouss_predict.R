#' Predict trajectory for OUSS
#'
#' @description
#' This function provides an easy way to predict the trajectory of a population under OUSS model
#'
#' @param yt A vector of log-abundance observations in the time-series.
#' @param tt A vector of observations times.
#' @param parms Parameters estimates (MLEs or REMLEs) from `ouss_remle()` or `ouss_mle()`
#' @param plot.it Logical argument of plot observed and predicted abundances in the time series
#'
#' @return
#' If "plot.it = TRUE", it will return a simple figure of the population trend.
#'
#' A list of:
#' 1) predicted values for each time point with observed data,
#' 2) predicted values including missing time points (not observed data)
#' @export
#'
#' @examples
#' yt1 = log(c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6))
#' tt1 = c(1966:1995)
#' parms1 = ouss_remle(yt = yt1, tt = tt1, fguess = guess_ouss(yt = yt1, tt = tt1))
#'
#' ouss_predict(yt = yt1, tt = tt1, parms = parms1$remles, plot.it = TRUE)
#'
ouss_predict <- function(yt,tt,parms, plot.it="TRUE"){

  t.i             <- tt-tt[1];
  q               <- length(t.i)-1;
  qp1             <- q+1;

  # parameters
  mu              <- parms[1];
  theta           <- parms[2];
  betasq          <- parms[3];
  tausq           <- parms[4];

  Var.inf         <- betasq/(2*theta);
  t.s             <- t.i[2:qp1] - t.i[1:q];
  t.cols          <- matrix(rep(t.i,each=qp1),nrow=qp1,ncol=qp1, byrow=FALSE);
  t.rows          <- t(t.cols);
  abs.diffs       <- abs(t.rows-t.cols);

  nmiss           <- t.s-1;
  long.nmiss      <- c(0,nmiss);
  Nmiss           <- sum(nmiss)

  long.t          <- t.i[1]:max(t.i)
  where.miss      <- which(is.na(match(x=long.t,table=t.i)),
                           arr.ind=TRUE)
  lt.cols         <- matrix(rep(long.t),
                            nrow=(qp1+Nmiss),
                            ncol=(qp1+Nmiss),
                            byrow=FALSE);
  lt.rows         <- t(lt.cols);
  labs.diffs      <- abs(lt.rows-lt.cols);

  Sigma.mat       <- Var.inf*exp(-theta*abs.diffs);
  Itausq          <- matrix(0,qp1,qp1);
  diag(Itausq)    <- rep(tausq,qp1);
  V               <- Sigma.mat+Itausq;

  long.V          <- Var.inf*exp(-theta*labs.diffs) + diag(rep(tausq,(qp1+Nmiss)))

  Predict.t       <- rep(0,qp1);
  Muvec           <- rep(mu,q);
  miss.predict    <- list()
  Muvec.miss      <- rep(mu,qp1);
  start.miss      <- 1
  stop.miss       <- 0
  for (tj in 1:qp1){
    Y.omitj       <- yt[-tj];    #  Omit observation at time tj.
    V.omitj       <- V[-tj,-tj];  #  Omit row tj and col tj from var-cov matrix.
    V12           <- V[tj,-tj];       #  Submatrix:  row tj without col tj.
    Predict.t[tj] <- mu+V12%*%ginv(V.omitj)%*%(Y.omitj-Muvec);  #  Graybill's 1976 Thm.

    if(long.nmiss[tj]==0){
      miss.predict[[tj]] <- Predict.t[tj]}else
        if(long.nmiss[tj]>0){

          start.miss <- stop.miss+1
          ntjmiss    <- long.nmiss[tj]
          mu.miss    <- rep(mu,ntjmiss);
          ind.tjmiss <- where.miss[start.miss:(start.miss+(ntjmiss-1))]
          stop.miss  <- stop.miss+ntjmiss

          longV12    <- long.V[ind.tjmiss,-where.miss]

          miss.predict[[tj]] <- c(mu.miss + longV12%*%ginv(V)%*%(yt-Muvec.miss), Predict.t[tj])
        }
  }

  Predict.t <- exp(Predict.t);
  LPredict.t <- exp(as.vector(unlist(miss.predict)))

  isinf <- sum(is.infinite(Predict.t))
  if(isinf>0){
    where.infs <- which(is.infinite(Predict.t)==TRUE, arr.ind=TRUE)
    Predict.t[where.infs] <- .Machine$double.xmax
  }

  isinf2 <- sum(is.infinite(LPredict.t))
  if(isinf2>0){
    where.infs <- which(is.infinite(LPredict.t)==TRUE, arr.ind=TRUE)
    LPredict.t[where.infs] <- .Machine$double.xmax
  }

  if(plot.it=="TRUE"){
    #  Plot the data & model-fitted values
    #X11()
    plot(tt,exp(yt),xlab="Time",ylab="Population abundance",type="b",cex=1.5,
         main="Predicted (--) and observed (-o-) abundances");#  Population data are circles.
    par(lty="dashed"); #  Predicted abundances are dashed line.
    points(tt,Predict.t, type="l", lwd=1);
  }

  return(list(cbind(tt,Predict.t,exp(yt)), cbind(long.t,LPredict.t) ))
}
