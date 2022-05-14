#' @title Schaefer's reference points minimizing function
#'
#' @description
#' Function used in the maximum likelihood minimization to estimate Schaefer's reference points
#' using lognormal distribution. The use of lognormal in maximum likelihood is important
#' since the index of abundance is assumed to follow lognormal distribution and
#' all the observation errors is the result of the relationships between stock biomass and index of
#' abundance which requires to be estimated (Polacheck et al., 1993).
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
#' sustainable yield), MSY (maximum sustainable yield), Emsy (effort at maximum sustainable yield), and B0 (biomass when fishing is started).
#' @param df dataframe containing three columns; year, catch and effort
#' @param OWT is CPUE plot showing One Way Trip pattern? The default is FALSE, but should be
#' replaced with TRUE when the plot shows One Way Trip
#' @param Frate exploitation rate collected from other survey. The default is 0.7
#' @param weight weight given to the deviation between observed and predicted value in
#' exploitation rate. The default is set at 100 and can be adjusted to create a more make sense result
#'
#' @return input for inpars are  kept at initial value without using log() like the other minimization inputs. If
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
#'
#' Bmsy <- K/2
#' MSY <- (r*K)/4
#' Emsy <- r/(2*q)
#'
#' ### Estimate parameters using optim
#' inpars <- c(Bmsy, MSY, Emsy, B0, 0.1)
#'
#' fit <- optim(par=inpars,
#'              fn=SparRP_min,
#'              df=df.goodcontrast,
#'              method="Nelder-Mead",
#'              OWT=FALSE, Frate = 0.7, weight = 100,
#'              hessian=TRUE)
#'
#' fitted_pars <- fit$par
#'
### compile the result and calculate the standard error
#' SparRP_vals <- data.frame(RPpar = c("Bmsy", "MSY", "Emsy", "B0", "sigma"),
#'                          init_pars = c(Bmsy, MSY, Emsy, B0, 0.1),
#'                          fitted_pars = fitted_pars,
#'                          std_err = sqrt(abs(diag(solve(-fit$hessian)))))
#'
#'
SparRP_min <- function(inpars, df, OWT=FALSE, Frate = 0.7, weight = 100){
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

  annualFrates <- df[,2]/EstBt # F = catch/biomass = Z-M

  if (OWT==FALSE){
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = sigma, log = TRUE))
  }
  else{
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = sigma, log = TRUE)) +
      weight * (tail(annualFrates,1) - Frate)^2
  }
  return(nll)
}
