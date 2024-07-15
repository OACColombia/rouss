#' Conduct Viable Population Monitoring - a risk-based approach
#'
#' @description
#' This function provides an easy way to conduct a risk-based viable population monitoring for a time-series data
#'
#' @param yt A vector of log-abundance observations in the time-series.
#' @param tt A vector of observations times.
#' @param ntraj Number of trajectories (e.g., 100, 1000) to simulate
#' @param N.critical Critical value of the population (recommended = 3/4*mean(observed in the time series))
#' @param plt Population Long Term simulations. How long you want to simulate the trajectories?
#'
#' @return A dataframe with the vector of observations times, vector of abundance observations, and two computed probability of extinctions (ntraj below N.critical/ntraj; area under the curve in log-normal distribution below N.critical)
#' @export
#'
#' @examples
#' yt1 = log(c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6))
#' tt1 = c(1966:1995)
#' ntraj = 100 #It is better if you simulate over 50000 (the probability stabilizes)
#' N.critical = 3/4*mean(exp(yt2))
#' plt = 5
#'
#' vpm(yt = yt1, tt = tt1, ntraj = ntraj, N.critical = N.critical, plt = plt)

vpm <- function(yt, tt, ntraj, N.critical, plt){
  # A vector of time steps to evaluate
  last.tt <- c(round(1/4*length(tt),0):length(tt))
  # Number of trajectories to simulate
  ntraj = ntraj
  # threshold time steps and length
  thres.times <- 0:plt
  len.sim <- max(thres.times)+1
  # Confirming decimals in N.critical
  N.critical = round(N.critical,3)
  #array to store p.below critical point (last point < N.critical/ntraj)
  probability.ext1 <- rep(0,length(last.tt))
  #array to store p.below critical point - area under the curve
  probability.ext2 <- rep(0,length(last.tt))
  #array to store model in each time
  Model <- rep(0,length(last.tt))

  for(t in last.tt){
    #array to store last point of each trajectory
    last.points <- rep(0,ntraj)
    #Compute Restricted Maximum Likelihood Estimaton for OUSS
    OUSS.partial <- ouss_remle(yt = yt[1:t],
                               tt = tt[1:t],
                               fguess = guess_ouss(yt = yt[1:t],
                                                   tt = tt[1:t]))
    #Select the model to simulate based on density dependence parameter
    model <- if(OUSS.partial$remles[2] < 0.025){"EGSS"}else{"OUSS"}

    Model[t] <- model;

    if(model == "OUSS"){
      for(j in 1:ntraj){
        #simulate the next PLT values (5, 10, 20,...)
        sim.mat <- ouss_sim(1,
                            tt = thres.times,
                            parms = OUSS.partial$remles)

        Pop.sim <- exp(c(yt[t], sim.mat[-1]))
        last.points[j] <- Pop.sim[len.sim]
      }

      probability.extinction1 <- length(which(round(last.points,3) <= N.critical))/ntraj

      probability.ext1[t] <- if(probability.extinction1 <= 0.001){"≤0.001"}else{probability.extinction1}

      kde.sims <- kde1d(x=last.points)
      probability.extinction2 <- round(pkde1d(q=N.critical, obj=kde.sims),3)

      probability.ext2[t] <- if(probability.extinction2 <= 0.001){"≤0.001"}else{probability.extinction2}
    }
    else{
      for(j in 1:ntraj){
        #Compute Maximum Likelihood Estimators for EGSS
        EGSS.partial <- egss_mle(yt = yt[1:t],
                                 tt = tt[1:t],
                                 fguess_egss = guess_egss(yt = yt[1:t],
                                                          tt = tt[1:t]))
        #simulate the next PLT values (5, 10, 20,...)
        sim.mat <- egss_sim(1,
                            tt = thres.times,
                            parms = EGSS.partial$mles)
        Pop.sim <- exp(c(yt[t], sim.mat[-1]))
        last.points[j] <- Pop.sim[len.sim]
      }

      probability.extinction1 <- length(which(round(last.points,3) <= N.critical))/ntraj

      probability.ext1[t] <- if(probability.extinction1 <= 0.001){"≤0.001"}else{probability.extinction1}

      kde.sims <- kde1d(x=last.points)
      probability.extinction2 <- round(pkde1d(q=N.critical, obj=kde.sims),3)

      probability.ext2[t] <- if(probability.extinction2 <= 0.001){"≤0.001"}else{probability.extinction2}
    }
  }

  out <- tibble(Time = tt,
                Model = Model,
                Prob.ext.Risk.1 = probability.ext1,
                Prob.ext.Risk.2 = probability.ext2,
                Observed = exp(yt))
  return(out)
}
