#' @title Schaefer's reference points minimizing function
#'
#' @description
#' Function used in the maximum likelihood minimization to estimate Schaefer's reference points.
#' Observation error is used to increase the accuracy of data fitting. It is assumed to occur in the relationship
#' between stock biomass and index of abundance and is estimated assuming lognormal distribution in maximum likelihood (Polacheck et al., 1993).
#'
#' Since fishing effort data collection are not always conducted regularly while catch is likely
#' have a better time series information, this function also allow for some lose of data.
#'
#' This function also consider the different quality data, for instance if the data shows
#' a one way trip pattern which losing increasing rate of increase.
#'
#' A continuation step in the following example section can help to estimate the standard error in the reference points value
#'
#' @param inpars reference point parameters which consist of Bmsy (stock biomass at maximum
#' sustainable yield), MSY (maximum sustainable yield), Emsy (effort at maximum sustainable yield), B0 (biomass when fishing is started)
#' and sigma (observation error)
#' @param df dataframe containing three columns; year, catch and effort. A fourth column with biomass should be added
#' if OWT (One Way Trip) option uses "Biomass"
#' @param OWT is CPUE plot showing One Way Trip pattern? The default is FALSE, but should be
#' replaced with either "Biomass" or "Depletion" when the plot shows One Way Trip type of data
#' @param currentF Current exploitation rate collected from other survey. The default is 0.7 to say that almost all fish are caught
#' @param weight weight given to the deviation between observed and predicted value in either
#' biomass or exploitation rate.
#'
#' @return
#' A penalized likelihood is used to fix the lack of contrast in One Way Trip type of data using Depletion or Biomass data.
#'
#' The Biomass option in OWT is used when biomass time series data from acoustic or trawl survey is available and
#' should be added as the fourth columns in the input dataframe. The default weight when Biomass level is set at 0.5 with range
#' between 0-1 (lower accuracy with high variance as closer to 0, constrain the estimation procedure to fit the auxiliary
#' information as closer to 1)
#'
#' The Depletion option in OWT uses current harvest rate from survey or expert knowledge as penalty. The default
#' weight for harvest rate is 10 and can be adjusted so the predicted harvest rate reach a closest value
#' to the current exploitation rate. Predicted harvest rate value in each iteration process will show up when optimization
#' process is being executed.
#'
#' Input for inpars are  kept at initial value without using log() like the other minimization inputs. If
#' the fitted parameters resulting in minus value, use the constrained variables and "L-BFGS-B" optimization method,
#' and produce the standard error from hessian using steps in https://stackoverflow.com/questions/27202395/how-do-i-get-standard-errors-of-maximum-likelihood-estimates-in-stan
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
#' sigma <- 0.1
#'
#' Bmsy <- K/2
#' MSY <- (r*K)/4
#' Emsy <- r/(2*q)
#'
#' ### Estimate parameters using optim
#' inpars <- c(Bmsy, MSY, Emsy, B0, sigma)
#'
#' fit <- optim(par=inpars,
#'              fn=ParRP_min,
#'              df=df.goodcontrast,
#'              method="Nelder-Mead",
#'              OWT=FALSE, currentF = 0.7, weight = 1000,
#'              hessian=TRUE)
#'
#' fitted_pars <- fit$par
#'
#' ### compile the result and calculate the standard error
#' ParRP_vals <- data.frame(RPpar = c("Bmsy", "MSY", "Emsy", "B0", "sigma"),
#'                          init_pars = c(Bmsy, MSY, Emsy, B0, 0.1),
#'                          fitted_pars = fitted_pars,
#'                          std_err = sqrt(abs(diag(solve(-fit$hessian)))))
#'
#'
ParRP_min <- function(inpars, df, OWT=FALSE, currentF = 0.7, weight = 0.5){
  Bmsy <- inpars[1]
  MSY <- inpars[2]
  Emsy <- inpars[3]
  B0 <- inpars[4]
  sigma <- inpars[5]

  K <- 2*Bmsy
  r <- (MSY*4)/K
  q <- r/(2*Emsy)

  CPUE <- df$catch/df$effort
  EstCatch <- EstCPUE <- EstBt <- vector(length=nrow(df))

  for (i in 1:nrow(df)) {
    if (i == 1) EstBt[i] <- B0
    if (i>1) EstBt[i] <- max(c(0.01,
                               EstBt[i-1] +  EstBt[i-1] * r * (1 - EstBt[i-1]/K) - df$catch[i-1]
    )
    )
  }
  EstCPUE <-  EstBt * q

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
