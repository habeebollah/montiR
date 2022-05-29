#' @title Schaefer's parameter estimation minimizing function
#'
#' @description
#' Function used in the maximum likelihood minimization to estimate Schaefer's reference points.
#' Observation error is assumed to occur in the relationship between stock biomass and index of abundance,
#' therefore it is estimated assuming lognormal distribution in maximum likelihood (Polacheck et al., 1993).
#'
#' Since fishing effort data collection are not always conducted regularly while catch is likely
#' have a better time series information, this function also allow for some lose of data.
#'
#' This function also consider the different quality data, for instance if the data shows
#' a one way trip pattern which losing increasing rate of increase.
#'
#' @param inpars surplus production parameters which consist of K (carrying capacity), B0
#' (biomass when fishing is started), r (intrinsic growth rate), q (catchability coefficient)
#' @param df dataframe containing three columns; year, catch and effort
#' @param OWT is CPUE plot showing One Way Trip pattern? The default is FALSE, but should be
#' replaced with either "Depletion" or "Biomass" when the plot shows One Way Trip where depletion or biomass
#' is used to complement the lacking contrast information from index of abundance (CPUE)
#' @param currentF Current exploitation rate collected from other survey. The default is 0.7 to say that almost all fish are caught
#' @param weight weight given to the deviation between observed and predicted value in either
#' exploitation rate or biomass level. The default weight for exploitation rate is set at 1, and can be increased to reach the
#' value closest to the current exploitation rate. The default weight for biomass level is set at 0.5 with range
#' between 0-1 (lower accuracy with high variance as closer to 0, constrain the estimation procedure to fit the auxiliary
#' information as closer to 1)
#'
#' @return
#' Using biomass as auxiliary information to complement the one way trip data requires time series data. This
#' information should be added as the fourth columns in the input dataframe.
#'
#' @export
#'
#' @references
#' Hilborn, Ray, and Carl J. Walters, eds. Quantitative fisheries stock assessment: choice,
#' dynamics and uncertainty. Springer Science & Business Media, 1992.
#'
#' Polacheck, T., Hilborn, R., and A.E. Punt. 1993. Fitting surplus production models:
#' Comparing methods and measuring uncertainty. Canadian Journal of Fisheries and Aquatic Sciences, 50: 2597-2607.
#'
#' @examples
#' K <- 1000
#' B0 <- K
#' r <- 0.2
#' q <- 0.00025
#'
#' ### Estimate parameters using optim
#' inpars <- c(log(K), log(B0), log(r), log(q), log(0.1))
#'
#' fit <- optim(par=inpars,
#'              fn=Spar_min,
#'              df=df.goodcontrast,
#'              method="Nelder-Mead",
#'              OWT=FALSE, currentF = 0.7, weight = 0.5)
#'
#' Spar_vals <- data.frame(SPpar = c("K", "B0", "r", "q", "sigma"),
#'                          init_pars = c(K, B0, r, q, 0.1),
#'                          fitted_pars = exp(fit$par))
#'
Spar_min1 <- function(inpars, df, OWT=FALSE, currentF = 0.7, weight = 0.5){
  K <- exp(inpars[1])
  B0 <- exp(inpars[2])
  r <- exp(inpars[3])
  q <- exp(inpars[4])
  sigma <- exp(inpars[5])

  CPUE <- df$catch/df$effort
  EstCatch <- EstCPUE <- EstBt <- vector(length=nrow(df))

  for (i in 1:nrow(df)) {
    if (i == 1) EstBt[i] <- B0
    if (i>1) EstBt[i] <- max(c(0.01,
                               EstBt[i-1] +  EstBt[i-1] * r * (1 - EstBt[i-1]/K) - df$catch[i-1]
    )
    )
  }
  EstCPUE <-  EstBt * q # check Punt & Hilborn (1996) p 17 & 29

  if (OWT==FALSE){
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = sigma, log = TRUE))
    }

  if (OWT=="Depletion"){
    annualFrates <- df[,2]/EstBt # F = catch/biomass = Z-M
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = sigma, log = TRUE)) +
      weight * (tail(annualFrates,1) - currentF)^2
    print(tail(annualFrates,1)) # to check whether the weighting makes the estimates depletion getting closer to the current exploitation rate
  }

  if (OWT=="Biomass"){ # haven't been checked yet
    surveyB <- df[,4]
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = sigma, log = TRUE)) +
      -sum(weight * (surveyB - EstBt)^2)
  }

  return(nll)
}

# the trajectory and stochastic can be found in Punt & Hilborn (1996) p 44
