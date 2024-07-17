#' Estimate parameters, compute predictions and run parametric bootstrapping
#'
#' @description
#' This function integrates several functions at once. It estimates the parameters using two methods (Maximum Likelihood "ML" or Restricted Maximum Likelihood "REML"), compute the predictions of abundance from the OUSS model, and run a parametric bootstrapping to estimate confidence intervals. In addition, it can generates figures of the bootstrapped parameter estimation and/or the trajectories (observed vs predicted).
#'
#' @param yt A vector of log-abundance observations in the time-series.
#' @param tt A vector of observations times.
#' @param pmethod Select one of two methods: ML or REML
#' @param nboot number of bootstrap
#' @param plot.pred Logical argument of plot observed and predicted abundances in the time series
#' @param plot.bootdists Logical argument of plot histograms of parametric bootstrapping of the four parameters
#'
#' @return A list of parameter estimates, log-likelihood estimate, AIC, matrix of bootstrapped values, parametric bootstrapped confidence intervals, and predictions (with CIs) of abundance (REMLEs or MLEs) for times with observations (1) and including gaps (2)
#' @export
#'
#' @examples
#' yt1 = log(c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6))
#' tt1 = c(1966:1995)
#' ouss_calc(yt = yt1, tt = tt1, pmethod = "REML", nboot = 100, plot.pred = TRUE, plot.bootdists = T)
#'
ouss_calc <- function(yt, tt, pmethod="ML", nboot, plot.pred="TRUE", plot.bootdists = "TRUE"){

  # Compute a rough guess of the parameter estimates to initialize the search:
  guesscalc <- guess_ouss(yt = yt, tt=tt)


  # Compute either the ML or the REML estimates, according to what you specified in point 1. above
  if(pmethod=="ML"){
    best.guess <- ouss_mle(yt = yt, tt = tt, fguess = guesscalc);
    AIC <- best.guess[[3]];
    remle.option <- "FALSE"}else if (pmethod=="REML"){
      best.guess <- ouss_remle(yt = yt, tt = tt, fguess = guesscalc);
      remle.option <- "TRUE"}else{
        print("Error: ML and REML are the only options allowed for 'method'")}

  # Parameter estimates and maximized log-likelihood (we will print these at the end)
  parms.est <- best.guess[[1]];
  lnLhat    <- best.guess[[2]];

  # Parametric bootstrap: computing both, parameters and predictions 95 % CI's
  pboot.calcs <- ouss_pboot(B=nboot,
                            yt = yt,
                            tt = tt,
                            REMLE = remle.option,
                            parms = parms.est,
                            plot.it = plot.bootdists);

  if(plot.pred=="TRUE"){
    plot(pboot.calcs$preds.CIs1[,1],
         pboot.calcs$preds.CIs1[,2],
         xlab = "Time", ylab = "Population abundance",
         type = "b", col = "blue", cex = 1.5);
    points(pboot.calcs$preds.CIs2[,1],
           pboot.calcs$preds.CIs2[,3],type="l",col="red", lty = "dashed");
    points(tt,pboot.calcs$preds.CIs1[,4],type="p",col="red");
    points(pboot.calcs$preds.CIs2[,1],
           pboot.calcs$preds.CIs2[,2],type="l",col="gray", lty = "dotted");
    points(pboot.calcs$preds.CIs2[,1],
           pboot.calcs$preds.CIs2[,4],type="l",col="gray", lty = "dotted");
    legend("top", legend=c("Observed", "Predicted", "CI"),
           col=c("blue","red", "gray"), lty=1:3, cex=0.8)
  }

  if(pmethod=="ML"){
    print("AIC score");
    print(AIC);
    out <- list(parms.est = parms.est,
                lnLhat = lnLhat,
                AIC = AIC,
                pbootmat = pboot.calcs[[1]],
                pboot.cis = pboot.calcs[[2]],
                pboot.preds1 = pboot.calcs$preds.CIs1,
                pboot.preds2 = pboot.calcs$preds.CIs2)}else{
                  out <- list(parms.est = parms.est,
                              lnLhat = lnLhat,
                              pbootmat = pboot.calcs[[1]],
                              pboot.cis = pboot.calcs[[2]],
                              pboot.preds1 = pboot.calcs$preds.CIs1,
                              pboot.preds2 = pboot.calcs$preds.CIs2)}

  return(out);
}
