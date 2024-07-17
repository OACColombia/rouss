#' Parametric bootstrapping for OUSS
#'
#' @description
#' This function conduct a parametric bootstrapping to estimate confidence intervals (2.5% to 97.5%) of parameters and estimates (MLEs and REMLEs) for OUSS model.
#'
#' @param B number of bootstraps
#' @param yt A vector of log-abundance observations in the time-series.
#' @param tt A vector of observations times.
#' @param parms Parameters estimates (MLEs or REMLEs) from `ouss_remle()` or `ouss_mle()`
#' @param REMLE Logical argument of use of REMLE
#' @param plot.it Logical argument of plot observed and predicted abundances in the time series
#'
#' @return list of bootsrtrapped values. If "plot.it = TRUE", it will return histograms of the four parameters
#' @export
#'
#' @examples
#' yt1 = log(c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6))
#' tt1 = c(1966:1995)
#' B1 = 10 #recommended for analysis ≥2000
#' parms1 = ouss_remle(yt = yt1, tt = tt1, fguess = guess_ouss(yt = yt1, tt = tt1))
#'
#' ouss_pboot(B = B1, yt = yt1, tt = tt1, parms = parms1$remle, REMLE=TRUE, plot.it = TRUE)
#'
ouss_pboot <- function(B, yt, tt, parms, REMLE="FALSE", plot.it="FALSE"){

  t.i             <- tt-tt[1];
  q               <- length(t.i)-1;
  qp1             <- q+1;

  # parameters
  mu              <- parms[1];
  theta           <- parms[2];
  betasq          <- parms[3];
  tausq           <- parms[4];

  #tt <- Tvec-Tvec[1]
  long.t <- t.i[1]:max(t.i)
  nparms    <- length(parms);
  preds.boot1<- matrix(0,nrow=B,ncol=length(t.i))
  preds.boot2<- matrix(0,nrow=B,ncol=length(long.t))

  if(REMLE=="TRUE"){

    boot.remles <- matrix(0,nrow=B,ncol=nparms+1);
    all.sims    <- ouss_sim(nsims=B,parms=parms,tt=tt);
    all.preds   <- ouss_predict(yt=yt,tt=tt,parms=parms,plot.it="FALSE")
    remle.preds <- all.preds[[1]][,2]
    remle.longpreds <- all.preds[[2]][,2]

    for(b in 1:B ){

      bth.timeseries <- all.sims[,b];
      remles.out <- ouss_remle(yt=bth.timeseries,tt=tt, fguess = parms);
      boot.remles[b,] <- c(remles.out$remles, remles.out$lnLhat);
      all.bootpreds <- ouss_predict( yt=bth.timeseries,
                                     tt=tt,
                                     parms=remles.out$remles,
                                     plot.it="FALSE");
      preds.boot1[b,] <- all.bootpreds[[1]][,2]
      preds.boot2[b,] <- all.bootpreds[[2]][,2]

    }

    CIs.mat <- apply(boot.remles,2,FUN=function(x){
      quantile(x,probs=c(0.025,0.975))});
    CIs.mat <- rbind(CIs.mat[1,1:4],parms,CIs.mat[2,1:4]);
    rownames(CIs.mat) <- c("Lower_CI_2.5","REMLE","Upper_CI_97.5");
    colnames(CIs.mat) <- c("mu", "theta","betasq","tausq");

    preds.CIs1 <- apply(preds.boot1,2,FUN=function(x){quantile(x,probs=c(0.025,0.975))});
    mean.boots <- apply(preds.boot1,2,FUN=function(x){quantile(x,probs=0.50)})
    preds.CIs1 <- t(rbind(tt,
                          exp(yt),
                          remle.preds-(mean.boots-preds.CIs1[1,]),
                          remle.preds,
                          remle.preds+(preds.CIs1[2,]-mean.boots)));
    colnames(preds.CIs1) <- c("Time","Observed","Lower_CI_2.5","REMLE","Upper_CI_97.5");

    preds.CIs2 <- apply(preds.boot2,2,FUN=function(x){quantile(x,probs=c(0.025,0.975))});
    mean.boots2 <- apply(preds.boot2,2,FUN=function(x){quantile(x,probs=0.50)})
    preds.CIs2 <- t(rbind(long.t+tt[1],
                          remle.longpreds-(mean.boots2-preds.CIs2[1,]),
                          remle.longpreds,
                          remle.longpreds+(preds.CIs2[2,]-mean.boots2)));

    #preds.CIs2 <- apply(preds.boot2,2,FUN=function(x){quantile(x,probs=c(0.025,0.5,0.975))});
    #preds.CIs2 <- t(rbind(long.t+Tvec[1],preds.CIs2));
    colnames(preds.CIs2) <- c("Time","Lower_CI_2.5","REMLE","Upper_CI_97.5");
    #pred.CIs2 <- cbind(preds.CIs2[,1], reml.longpreds-(preds.CIs2[,3]-preds.CIs2[,2]),reml.longpreds,reml.longpreds+(preds.CIs2[,4]-preds.CIs2[,3]))

    boot.list <- list(boot.remles = boot.remles, CIs.mat = CIs.mat, preds.CIs1 = preds.CIs1,
                      preds.CIs2=preds.CIs2);

    if(plot.it=="TRUE"){
      par(mfrow=c(2,2));
      hist(boot.remles[,1],main=expression(hat(mu)),xlab="");
      abline(v=parms[1],lwd=2,col="red");
      hist(boot.remles[,2],main=expression(hat(theta)),xlab="");
      abline(v=parms[2],lwd=2,col="red");
      hist(boot.remles[,3],main=expression(hat(beta^2)),xlab="");
      abline(v=parms[3],lwd=2,col="red");
      hist(boot.remles[,4],main=expression(hat(tau^2)),xlab="");
      abline(v=parms[4],lwd=2,col="red");
      par(mfrow=c(1,1));
    }
    return(boot.list)

  }else{

    boot.mles <- matrix(0,nrow=B,ncol=nparms+2);
    all.sims  <- ouss_sim(nsims=B, tt=tt, parms=parms);
    all.preds <- ouss_predict(yt=yt, tt=tt, parms=parms, plot.it="FALSE")
    mle.preds<- all.preds[[1]][,2]
    mle.longpreds <- all.preds[[2]][,2]

    for(b in 1:B ){

      bth.timeseries <- all.sims[,b];
      mles.out <- ouss_mle(yt = bth.timeseries, tt=tt, fguess = parms);
      boot.mles[b,] <- c(mles.out$mles, mles.out$lnL.hat,mles.out$AIC);
      all.bootpreds <- ouss_predict(yt=bth.timeseries, tt=tt, parms=mles.out$mles, plot.it="FALSE");
      preds.boot1[b,] <- all.bootpreds[[1]][,2]
      preds.boot2[b,] <- all.bootpreds[[2]][,2]
    }

    CIs.mat <- apply(boot.mles,2,FUN=function(x){quantile(x,probs=c(0.025,0.975))});
    CIs.mat <- rbind(CIs.mat[1,1:4],
                     parms,
                     CIs.mat[2,1:4]);
    rownames(CIs.mat) <- c("Lower_CI_2.5","MLE","Upper_CI_97.5");
    colnames(CIs.mat) <- c("mu", "theta","betasq","tausq");

    #preds.CIs1 <- apply(preds.boot1,2,FUN=function(x){quantile(x,probs=c(0.025,0.975))});
    #preds.CIs1 <- t(rbind(Tvec,preds.CIs1[1,], ml.preds, preds.CIs1[2,]));
    preds.CIs1 <- apply(preds.boot1,2,FUN=function(x){quantile(x,probs=c(0.025,0.5,0.975))});
    preds.CIs1 <- t(rbind(tt,
                          exp(yt),
                          preds.CIs1));
    colnames(preds.CIs1) <- c("Time", "Observed","Lower_CI_2.5","MLE","Upper_CI_97.5");
    pred.CIs1 <- cbind(preds.CIs1[,1],
                       mle.preds-(preds.CIs1[,3]-preds.CIs1[,2]),
                       mle.preds,
                       mle.preds+(preds.CIs1[,4]-preds.CIs1[,3]))

    #preds.CIs2 <- apply(preds.boot2,2,FUN=function(x){quantile(x,probs=c(0.025,0.975))});
    #preds.CIs2 <- t(rbind(Tvec,preds.CIs2[1,], ml.preds, preds.CIs2[2,]));
    preds.CIs2 <- apply(preds.boot2,2,FUN=function(x){quantile(x,probs=c(0.025,0.5,0.975))});
    preds.CIs2 <- t(rbind(long.t+tt[1],
                          preds.CIs2));
    pred.CIs2 <- cbind(preds.CIs2[,1],
                       mle.longpreds-(preds.CIs2[,3]-preds.CIs2[,2]),
                       mle.longpreds,
                       mle.longpreds+(preds.CIs2[,4]-preds.CIs2[,3]))

    colnames(preds.CIs2) <- c("Time","Lower_CI_2.5","MLE","Upper_CI_97.5");

    boot.list <- list(boot.mles = boot.mles,
                      CIs.mat = CIs.mat,
                      preds.CIs1 = preds.CIs1,
                      preds.CIs2 = preds.CIs2)

    if(plot.it=="TRUE"){
      par(mfrow=c(2,2));
      hist(boot.mles[,1],main=expression(hat(mu)),xlab="");
      abline(v=parms[1],lwd=2,col="red");
      hist(boot.mles[,2],main=expression(hat(theta)),xlab="");
      abline(v=parms[2],lwd=2,col="red");
      hist(boot.mles[,3],main=expression(hat(beta^2)),xlab="");
      abline(v=parms[3],lwd=2,col="red");
      hist(boot.mles[,4],main=expression(hat(tau^2)),xlab="");
      abline(v=parms[4],lwd=2,col="red");
      par(mfrow=c(1,1));
    }

    return(boot.list)

  }# End if/else

}
