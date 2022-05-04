#' @title Fox's reference points minimizing function
#'
#' @description
#' Function used in the maximum likelihood minimization to estimate Fox's reference points using lognormal distribution. The use of maximum likelihood is important to remove potential error during data collection.
#' Since fishing effort data collection are not always conducted regularly while catch is likely have a better time series information, this function also allow for some lose of data.
#'
#' This function also consider the different quality data, for instance if the data shows a one way trip pattern which losing increasing rate of increase.
#'
#' A continuation step in the following example section can help to estimate the standard error in the reference points value
#'
#' @param inpars reference point parameters which consist of Bmsy (stock biomass at maximum sustainable yield), MSY (maximum sustainable yield), Emsy (effort at maximum sustainable yield), and B0 (biomass when fishing is started).
#' @param df dataframe containing three columns; year, catch and effort
#' @param OWT is CPUE plot showing One Way Trip pattern? The default is FALSE, but should be replaced with TRUE when the plot shows One Way Trip
#' @param Frate exploitation rate collected from other survey. The default is 0.7
#' @param weight weight given to the deviation between observed and predicted value in exploitation rate. The default is set at 100 and can be adjusted to create a more make sense result
#'
#' @return
#' @export
#'
#' @references
#' Hilborn, Ray, and Carl J. Walters, eds. Quantitative fisheries stock assessment: choice, dynamics and uncertainty. Springer Science & Business Media, 1992.
#'
#' @examples
#' K <- 1000
#' B0 <- K
#' r <- 0.2
#' q <- 0.00025
#'
#' Bmsy <- K/exp(1)
#' MSY <- (r*K)/(exp(1)*log(K))
#' Emsy <- r/q
#'
#' ### Estimate parameters using optim
#' startPars <- c(log(Bmsy), log(MSY), log(Emsy), log(B0), log(0.1))
#'
#' fit <- optim(par=startPars,
#'              fn=RPparF_opt,
#'              df=goodcontrast,
#'              method="Nelder-Mead",
#'              OWT=FALSE, Frate = 0.7, weight = 100,
#'              hessian=TRUE)
#'
#' fitted_pars <- exp(fit$par)
#'
### compile the result and calculate the standard error
#' RPparF_val <- data.frame(RPpar = c("Bmsy", "MSY", "Emsy", "B0", "sigma"),
#'                          init_pars = c(Bmsy, MSY, Emsy, B0, 0.1),
#'                          fitted_pars = fitted_pars,
#'                          std_err = sqrt(abs(diag(solve(-fit$hessian))))) # need to check whether the the hessian should be back transformed using exp
#'
#'
RPparS_opt <- function(inpars, df, OWT=FALSE, Frate = 0.7, weight = 100){
  Bmsy <- exp(inpars[1])
  MSY <- exp(inpars[2])
  Emsy <- exp(inpars[3])
  B0 <- exp(inpars[4])
  sigma <- exp(inpars[5])

  K <- exp(1)*Bmsy
  r <- (MSY*exp(1)*log(K))/K
  q <- r/Emsy

  CPUE <- df$catch/df$effort
  EstCatch <- EstCPUE <- EstBt <- vector(length=nrow(df))

  for (i in 1:nrow(df)) {
    if (i == 1) EstBt[i] <- B0
    if (i>1) EstBt[i] <- max(c(0.01,
                               EstBt[i-1] +  EstBt[i-1] * r * (1 - log(EstBt[i-1])/log(K)) - df$catch[i-1]
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
