
# rouss

<!-- badges: start -->
<!-- badges: end -->

The goal of `rouss` is to fit the Ornstein-Uhlenbeck State-Space model to time-series data published by Dennis & Ponciano (2014). It also includes fit of Exponential Growth State-Space model from Humbert et al. (2009). We used these functions with community science data of the eBird platform. Finally, we include the risk-based Viable Population Monitoring (or Population Viability Monitoring) to assess probability of (quasi)extinction of the population across time. 

## Installation

You can install the development version of `rouss` like so:

``` r
library(devtools);
install_github("OACColombia/rouss")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(rouss)

#### ~~~ Some data ~~~ ####

#Observed population counts of American Redstart between 1966–1995, from the North American Breeding Bird Survey

#Discrete equal sampling
yt1 = log(c(18, 10, 9, 14, 17, 14, 
            5, 10, 
            9, 5, 11, 11, 4, 5, 4, 8, 2, 3, 9, 2, 
            4, 7, 4, 1, 
            2, 4, 11, 11, 9, 6)) #note, there are not zeros
tt1 = c(1966:1995)

#Arbitrarialy removing values (unequal sampling)
yt2 = log(c(18, 10, 9, 14, 17, 14,
            #5, 10, 
            9, 5, 11, 11, 4, 5, 4, 8, 2, 3, 9, 2, 
            #4, 7, 4, 1, 
            2, 4, 11, 11, 9, 6)) #note, there are not zeros
tt2 = c(1966:1971, 
        1974:1985,
        1990:1995)

#### ~~~ Applying OUSS model ~~~ ####

# Example with American Redstar without gaps in the time series (see Fig. 1 in Dennis et al. 2006)

AmericanRedstarNOGaps <- ouss_calc(yt = yt1, 
                                 tt = tt1, 
                                 pmethod = "REML", 
                                 nboot = 100, 
                                 plot.pred = TRUE, 
                                 plot.bootdists = T)


# Compare with the estimates of the altered dataset with missing data

AmericanRedstarGaps <- ouss_calc(yt = yt2, 
                                 tt = tt2, 
                                 pmethod = "REML", 
                                 nboot = 100, 
                                 plot.pred = TRUE, 
                                 plot.bootdists = T)
```

