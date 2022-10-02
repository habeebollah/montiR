#' @title Schaefer's data fitting and plotting
#'
#' @description
#' This function is used to predict the initial Biomass Dynamic Model parameters as input to do
#' the estimation in the latter process.
#'
#' It also works to check the result of estimated parameters by plotting the data, whether the
#' data and estimated parameter is fitted.
#'
#' @param K surplus production parameters which represents carrying capacity
#' @param B0 surplus production parameters which represents biomass when no fishing
#' activity has started
#' @param r surplus production parameters which represents intrinsic growth rate
#' @param q surplus production parameters which represents catchability coefficient
#' @param df dataframe containing three columns; year, catch and unit of effort
#' @param res option to show estimated data based on estimated surplus production input
#'
#' @export
#'
#' @references
#' Hilborn, Ray, and Carl J. Walters, eds. Quantitative fisheries stock assessment: choice,
#' dynamics and uncertainty. Springer Science & Business Media, 1992.
#'
#' @examples
#'
#' Par.init(K = 1000, B0 = 1000, r = 0.2, q = 0.00025, df = df.goodcontrast)
#'
Par.init <- function(K, B0, r, q, df, res=TRUE){

  CPUE <- df$cpue
  EstEffort <- EstCatch <- EstCPUE <- EstBt <- vector(length=nrow(df))
  EstEffort <- df$catch / df$cpue

  for (i in 1:nrow(df)) {
    if (i == 1) EstBt[i] <- B0
    if (i>1) EstBt[i] <- EstBt[i-1] +  EstBt[i-1] * r * (1 - EstBt[i-1]/K) - df$catch[i-1]
    EstCatch[i] <- EstBt[i] * q * EstEffort[i]
  }

  EstCPUE <- EstBt * q

  plot(CPUE, xlab="Year", ylab="CPUE", type="b", col="blue")
  lines(EstCPUE, type="b", col="red")
  legend("topright", legend=c("Observation", "Estimation"), lty=c(1, 1), col=c("blue", "red"), box.lty=1, cex=0.7, bty="n")

  if(res==TRUE){
    return(data.frame(CPUE=CPUE, EstBt=EstBt, EstCatch=EstCatch, EstCPUE=EstCPUE))
  }
}


#' @title Function used in the minimizing process to calculate Schaefer's management parameters
#'
#' @description
#' This function calculates Schaefer's management parameters using maximum likelihood minimization.
#'
#' @param inpars surplus production parameters which consist of K (carrying capacity), B0
#' (biomass when fishing is started), r (intrinsic growth rate), q (catchability coefficient),
#' s.sigma (observation error)
#' @param df dataframe containing three columns; year, catch and unit of effort. A fourth column with biomass should be added
#' if OWT (One Way Trip) option uses "Biomass"
#' @param OWT is CPUE plot showing One Way Trip pattern? The default is FALSE, but should be
#' replaced with either "Biomass" or "Depletion" when the plot shows One Way Trip type of data
#' @param currentF Current exploitation rate collected from other survey. Only being used
#' when OWT="Depletion" is chosen
#' @param weight weight given to the deviation between observed and predicted value in either
#' biomass or exploitation rate.
#'
#' @return
#' A penalized likelihood is used to fix the lack of contrast in One Way Trip type of data using Depletion or Biomass data.
#'
#' The Biomass option in OWT is used when biomass time series data from acoustic or trawl survey is available and
#' should be added as the fourth columns in the input dataframe. The default weight when Biomass level is set at 0.9 with range
#' between 0-1 (lower accuracy with high variance as closer to 0, constrain the estimation procedure to fit the auxiliary
#' information as closer to 1)
#'
#' The Depletion option in OWT uses current harvest rate from survey or expert knowledge as penalty.
#' Depletion range is between 0 to 1, where higher number represent higher depletion level. The default is 0.7 to
#' say that the depletion is high and many fish were caught. The default weight for harvest rate is 1000 and can be adjusted
#' so the predicted harvest rate reach a closest value to the current exploitation rate. Predicted harvest rate value
#' in each optimization step will show up when optimization process is being executed.
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
#'
Par.min <- function(inpars, df, OWT=FALSE, currentF = 0.7, weight = 1000){
  K <- exp(inpars[1])
  B0 <- exp(inpars[2])
  r <- exp(inpars[3])
  q <- exp(inpars[4])
  s.sigma <- exp(inpars[5])

  CPUE <- df$cpue
  EstCPUE <- EstBt <- vector(length=nrow(df))

  for (i in 1:nrow(df)) {
    if (i == 1) EstBt[i] <- B0
    if (i>1) EstBt[i] <- max(c(0.01,
                               EstBt[i-1] +  EstBt[i-1] * r * (1 - EstBt[i-1]/K) - df$catch[i-1]))
  }
  EstCPUE <-  EstBt * q

  if (OWT==FALSE){
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = s.sigma, log = TRUE))
    }

  if (OWT=="Depletion"){
    annualFrates <- df[,2]/EstBt # F = catch/biomass = Z-M
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = s.sigma, log = TRUE)) +
      weight * (tail(annualFrates,1) - currentF)^2
    print(tail(annualFrates,1)) # to check whether the weighting makes the estimates depletion getting closer to the current exploitation rate
  }

  if (OWT=="Biomass"){
    surveyB <- df[,4]
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = s.sigma, log = TRUE)) +
      sum(weight * (surveyB - EstBt)^2)
  }

  return(nll)
}

#' @title Calculating Schaefer's management parameters
#'
#' @description
#' This function calculates Schaefer's management parameters using maximum likelihood minimization.
#' Observation error is used to increase the accuracy of data fitting. It is assumed to occur in the relationship
#' between stock biomass and index of abundance and is estimated assuming lognormal distribution in
#' maximum likelihood (Polacheck et al., 1993).
#'
#' Since fishing effort data collection are not always conducted regularly while catch is likely
#' have a better time series information, this function also allow for some lose of catch and effort data.
#'
#' This function also consider the different quality data, for instance if the data shows
#' a one way trip pattern which losing rate of catch increase.
#'
#' @param K surplus production parameter which represents carrying capacity
#' @param B0 surplus production parameter which represents biomass when no fishing
#' activity has started
#' @param r surplus production parameter which represents intrinsic growth rate
#' @param q surplus production parameter which represents catchability coefficient
#' @param s.sigma surplus production parameter which represents observation error
#' @param df dataframe containing three columns; year, catch and unit of effort. A fourth column with biomass should be added
#' if OWT (One Way Trip) option uses "Biomass"
#' @param OWT is CPUE plot showing One Way Trip pattern? The default is FALSE, but should be
#' replaced with either "Biomass" or "Depletion" when the plot shows One Way Trip type of data
#' @param currentF Current exploitation rate collected from other survey. Only being used
#' when OWT="Depletion" is chosen
#' @param weight weight given to the deviation between observed and predicted value in either
#' biomass or exploitation rate
#' @param plot option to show the plot as result of fitted estimation to observation data
#'
#' @return
#' A penalized likelihood is used to fix the lack of contrast in One Way Trip type of data using Depletion or Biomass data.
#'
#' The Biomass option in OWT is used when biomass time series data from acoustic or trawl survey is available and
#' should be added as the fourth columns in the input dataframe. The default weight when Biomass level is set at 0.9 with range
#' between 0-1 (lower accuracy with high variance as closer to 0, constrain the estimation procedure to fit the auxiliary
#' information as closer to 1)
#'
#' The Depletion option in OWT uses current harvest rate from survey or expert knowledge as penalty.
#' Depletion range is between 0 to 1, where higher number represent higher depletion level. The default is 0.7 to
#' say that the depletion is high and many fish were caught. The default weight for harvest rate is 1000 and can be adjusted
#' so the predicted harvest rate reach a closest value to the current exploitation rate. Predicted harvest rate value
#' in each optimization step will show up when optimization process is being executed.
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
#'
#' calc.MSY(K=1000, B0=1000, r=0.2, q=0.00025, s.sigma=0.1, df=df.goodcontrast,
#' OWT=FALSE, currentF = 0.7, weight=0.5, plot=TRUE)
#'

calc.MSY <- function(K, B0, r, q, s.sigma, df,
                     OWT = FALSE, currentF = 0.7, weight = 0.9, plot=FALSE) {
  inpars <- c(log(K), log(B0), log(r), log(q), log(s.sigma))

  fit <- optim(par = inpars,
               fn = Par.min,
               df = df,
               method = "Nelder-Mead",
               OWT = OWT,
               currentF = currentF,
               weight = weight
  )

  vals <- exp(fit$par)
  res <- list("Parameter" = data.frame("SPpar" = c("K", "B0", "r", "q", "s.sigma"),
                                       "fitted_pars" = vals
  ),
  "MSY" = data.frame("MSY" = (vals[3] * vals[1]) / 4, # (r*K)/4
                     "Emsy" = vals[3] / (2 * vals[4]), # r/(2*q)
                     "Bmsy" = vals[1] / 2, # K/2
                     "E.rate_MSY" = (tail(df[, 2], 1) / (vals[3] * vals[1] / 4)) * 100, # catch/MSY
                     "E.rate_Emsy" = (tail(df[, 2], 1) / (vals[3] / (2 * vals[4])) * 100 # catch/Emsy
                     )
  )
  )
  if(plot==TRUE){
    Par.init(K=vals[1], B0=vals[2], r=vals[3], q=vals[4], df=df, res=TRUE)
  }
  return(res)
}

#' @title Function used in the minimizing process to calculate standard error in Schaefer's management parameter
#'
#' @description
#' This function calculates the standard error in Schaefer's management parameters using maximum likelihood minimization.
#'
#' @param inpars management parameters which consist of Bmsy (stock biomass at maximum
#' sustainable yield), MSY (maximum sustainable yield), Emsy (effort at maximum sustainable yield), B0 (biomass when fishing is started)
#' and s.sigma (observation error)
#' @param df dataframe containing three columns; year, catch and unit of effort. A fourth column with biomass should be added
#' if OWT (One Way Trip) option uses "Biomass"
#' @param OWT is CPUE plot showing One Way Trip pattern? The default is FALSE, but should be
#' replaced with either "Biomass" or "Depletion" when the plot shows One Way Trip type of data
#' @param currentF Current exploitation rate collected from other survey.
#' @param weight weight given to the deviation between observed and predicted value in either
#' biomass or exploitation rate.
#'
#' @return
#' A penalized likelihood is used to fix the lack of contrast in One Way Trip type of data using Depletion or Biomass data.
#'
#' The Biomass option in OWT is used when biomass time series data from acoustic or trawl survey is available and
#' should be added as the fourth columns in the input dataframe. The default weight when Biomass level is set at 0.9 with range
#' between 0-1 (lower accuracy with high variance as closer to 0, constrain the estimation procedure to fit the auxiliary
#' information as closer to 1)
#'
#' The Depletion option in OWT uses current harvest rate from survey or expert knowledge as penalty.
#' Depletion range is between 0 to 1, where higher number represent higher depletion level. The default is 0.7 to
#' say that the depletion is high and many fish were caught. The default weight for harvest rate is 1000 and can be adjusted
#' so the predicted harvest rate reach a closest value to the current exploitation rate. Predicted harvest rate value
#' in each optimization step will show up when optimization process is being executed.
#'
#' Input are  kept at initial value without using log() like the other minimization inputs. If
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
ParRP.min <- function(inpars, df, OWT=FALSE, currentF = 0.7, weight = 0.5){
  MSY <- inpars[1]
  Emsy <- inpars[2]
  Bmsy <- inpars[3]
  s.sigma <- inpars[4]

  K <- 2*Bmsy
  B0 <- K
  r <- (MSY*4)/K
  q <- r/(2*Emsy)

  CPUE <- df$cpue
  EstCatch <- EstCPUE <- EstBt <- vector(length=nrow(df))

  for (i in 1:nrow(df)) {
    if (i == 1) EstBt[i] <- B0
    if (i>1) EstBt[i] <- max(c(0.01,
                               EstBt[i-1] +  EstBt[i-1] * r * (1 - EstBt[i-1]/K) - df$catch[i-1]))
  }
  EstCPUE <-  EstBt * q

  if (OWT==FALSE){
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = s.sigma, log = TRUE))
  }

  if (OWT=="Depletion"){
    annualFrates <- df[,2]/EstBt # F = catch/biomass = Z-M
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = s.sigma, log = TRUE)) +
      weight * (tail(annualFrates,1) - currentF)^2
    print(tail(annualFrates,1)) # to check whether the weighting makes the estimates depletion getting closer to the current exploitation rate
  }

  if (OWT=="Biomass"){
    surveyB <- df[,4]
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = s.sigma, log = TRUE)) +
      sum(weight * (surveyB - EstBt)^2)
  }

  return(nll)
}

#' @title Calculating standard error in Schaefer's management parameter
#'
#' @description
#' This function calculates the standard error in Schaefer's management parameters using maximum likelihood minimization.
#' Observation error is used to increase the accuracy of data fitting. It is assumed to occur in the relationship
#' between stock biomass and index of abundance and is estimated assuming lognormal distribution in
#' maximum likelihood (Polacheck et al., 1993).
#'
#' Since fishing effort data collection are not always conducted regularly while catch is likely
#' have a better time series information, this function also allow for some lose of catch and effort data.
#'
#' This function also consider the different quality data, for instance if the data shows
#' a one way trip pattern which losing rate of catch increase.
#'
#' @param MSY management parameter representing Maximum Sustainable Yield
#' @param Emsy management parameter representing Effort at maximum sustainable yield
#' @param Bmsy management parameter representing stock Biomass at maximum sustainable yield
#' @param s.sigma surplus production parameter which represents observation error
#' @param df dataframe containing three columns; year, catch and unit of effort. A fourth column with biomass should be added
#' if OWT (One Way Trip) option uses "Biomass"
#' @param OWT is CPUE plot showing One Way Trip pattern? The default is FALSE, but should be
#' replaced with either "Biomass" or "Depletion" when the plot shows One Way Trip type of data
#' @param currentF Current exploitation rate collected from other survey.
#' @param weight weight given to the deviation between observed and predicted value in either
#' biomass or exploitation rate.
#'
#' @return
#' A penalized likelihood is used to fix the lack of contrast in One Way Trip type of data using Depletion or Biomass data.
#'
#' The Biomass option in OWT is used when biomass time series data from acoustic or trawl survey is available and
#' should be added as the fourth columns in the input dataframe. The default weight when Biomass level is set at 0.9 with range
#' between 0-1 (lower accuracy with high variance as closer to 0, constrain the estimation procedure to fit the auxiliary
#' information as closer to 1)
#'
#' The Depletion option in OWT uses current harvest rate from survey or expert knowledge as penalty.
#' Depletion range is between 0 to 1, where higher number represent higher depletion level. The default is 0.7 to
#' say that the depletion is high and many fish were caught. The default weight for harvest rate is 1000 and can be adjusted
#' so the predicted harvest rate reach a closest value to the current exploitation rate. Predicted harvest rate value
#' in each optimization step will show up when optimization process is being executed.
#'
#' Input are  kept at initial value without using log() like the other minimization inputs. If
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
#'
#' calc.SE(MSY=50, Emsy=400, Bmsy=500, s.sigma=0.1, df=df.onewaytrip,
#' OWT = FALSE, currentF = 0.7, weight = 0.9)
#'

calc.SE <- function(MSY, Emsy, Bmsy, s.sigma,
                    df, OWT = FALSE, currentF = 0.7, weight = 0.9) {
  inpars <- c(MSY, Emsy, Bmsy, s.sigma)

  fit <- optim(par = inpars,
               fn = ParRP.min,
               df = df,
               method = "Nelder-Mead",
               OWT = OWT,
               currentF = currentF,
               weight = weight,
               hessian = TRUE
  )

  res <- data.frame(ManagPar = c("Bmsy", "MSY", "Emsy", "sigma"),
                    fitted_pars = fit$par,
                    std_err = sqrt(diag(abs(solve(-fit$hessian)))
                    )
  )
  return(res)
}

#' @title Function used in the minimizing process to calculate likelihood profile in Schaefer's management parameter
#'
#' @description
#' This function calculates the likelihood profile as input to calculate confidence interval
#' in Schaefer's management parameters using maximum likelihood minimization.
#'
#' @param inpar input parameter which uses r (intrinsic growth)
#' @param MSYval consist of MSY (maximum sustainable yield) value
#' @param df dataframe containing three columns; year, catch and unit of effort. A fourth column with biomass should be added
#' if OWT (One Way Trip) option uses "Biomass"
#' @param OWT is CPUE plot showing One Way Trip pattern? The default is FALSE, but should be
#' replaced with either "Biomass" or "Depletion" when the plot shows One Way Trip type of data
#' @param currentF Current exploitation rate collected from other survey.
#' @param weight weight given to the deviation between observed and predicted value in either
#' biomass or exploitation rate.
#'
#' @return
#' A penalized likelihood is used to fix the lack of contrast in One Way Trip type of data using Depletion or Biomass data.
#'
#' The Biomass option in OWT is used when biomass time series data from acoustic or trawl survey is available and
#' should be added as the fourth columns in the input dataframe. The default weight when Biomass level is set at 0.9 with range
#' between 0-1 (lower accuracy with high variance as closer to 0, constrain the estimation procedure to fit the auxiliary
#' information as closer to 1)
#'
#' The Depletion option in OWT uses current harvest rate from survey or expert knowledge as penalty.
#' Depletion range is between 0 to 1, where higher number represent higher depletion level. The default is 0.7 to
#' say that the depletion is high and many fish were caught. The default weight for harvest rate is 1000 and can be adjusted
#' so the predicted harvest rate reach a closest value to the current exploitation rate. Predicted harvest rate value
#' in each optimization step will show up when optimization process is being executed.
#'
#' Input are  kept at initial value without using log() like the other minimization inputs. If
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
#' Punt, A. E., & Hilborn, R. 1996. Biomass dynamic models. FAO Computerized Information Series Fisheries, 10, 1-62.
#'
MSYprofile.min <- function(inpar, df, MSYval, OWT=FALSE, currentF = 0.7, weight = 0.5){

  r <- inpar
  MSY <- MSYval
  K <- (MSY*4)/r
  B0 <- K
  sigma <- 0.1 # should the sigma be used here?

  CPUE <- df$cpue
  EstCatch <- EstCPUE <- EstBt <- vector(length=nrow(df))

  for (i in 1:nrow(df)) {
    if (i == 1) EstBt[i] <- B0
    if (i>1) EstBt[i] <- max(c(0.01,
                               EstBt[i-1] +  EstBt[i-1] * r * (1 - EstBt[i-1]/K) - df$catch[i-1]))
  }

  q <- mean(CPUE) / mean(EstBt)
  EstCPUE <-  EstBt * q

  if (OWT==FALSE){
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = sigma, log = TRUE))
  }

  if (OWT=="Depletion"){
    annualFrates <- df[,2]/EstBt # F = catch/biomass = Z-M
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = sigma, log = TRUE)) +
      weight * (tail(annualFrates,1) - currentF)^2
    #print(tail(annualFrates,1)) # to check whether the weighting makes the estimates depletion getting closer to the current exploitation rate
  }

  if (OWT=="Biomass"){
    surveyB <- df[,4]
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = sigma, log = TRUE)) +
      -sum(weight * (surveyB - EstBt)^2)
  }

  return(nll)
}

#' @title Calculating confident interval for MSY in Schaefer's management parameter
#'
#' @description
#' This function calculates the confident interval for MSY in Schaefer's management parameters
#' using likelihood profile.
#' Observation error is used to increase the accuracy of data fitting. It is assumed to occur in the relationship
#' between stock biomass and index of abundance and is estimated assuming lognormal distribution in
#' maximum likelihood (Polacheck et al., 1993).
#'
#' Since fishing effort data collection are not always conducted regularly while catch is likely
#' have a better time series information, this function also allow for some lose of catch and effort data.
#'
#' This function also consider the different quality data, for instance if the data shows
#' a one way trip pattern which losing rate of catch increase.
#'
#' @param MSYval management parameter representing Maximum Sustainable Yield (MSY)
#' @param rval parameter which uses r (intrinsic growth)
#' @param df dataframe containing three columns; year, catch and unit of effort. A fourth column with biomass should be added
#' if OWT (One Way Trip) option uses "Biomass"
#' @param OWT is CPUE plot showing One Way Trip pattern? The default is FALSE, but should be
#' replaced with either "Biomass" or "Depletion" when the plot shows One Way Trip type of data
#' @param currentF Current exploitation rate collected from other survey.
#' @param weight weight given to the deviation between observed and predicted value in either
#' biomass or exploitation rate.
#' @param plot option to show the plot as result of likelihood profile estimation
#'
#' @return
#' A penalized likelihood is used to fix the lack of contrast in One Way Trip type of data using Depletion or Biomass data.
#'
#' The Biomass option in OWT is used when biomass time series data from acoustic or trawl survey is available and
#' should be added as the fourth columns in the input dataframe. The default weight when Biomass level is set at 0.9 with range
#' between 0-1 (lower accuracy with high variance as closer to 0, constrain the estimation procedure to fit the auxiliary
#' information as closer to 1)
#'
#' The Depletion option in OWT uses current harvest rate from survey or expert knowledge as penalty.
#' Depletion range is between 0 to 1, where higher number represent higher depletion level. The default is 0.7 to
#' say that the depletion is high and many fish were caught. The default weight for harvest rate is 1000 and can be adjusted
#' so the predicted harvest rate reach a closest value to the current exploitation rate. Predicted harvest rate value
#' in each optimization step will show up when optimization process is being executed.
#'
#' Input are  kept at initial value without using log() like the other minimization inputs. If
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
#' Punt, A. E., & Hilborn, R. 1996. Biomass dynamic models. FAO Computerized Information Series Fisheries, 10, 1-62.
#'
#' @examples
#'
#' calc.CI(MSYval= 50, rval= 0.2, df=df.goodcontrast, plot=TRUE)
#'

calc.CI <- function(MSYval, rval, df, OWT=FALSE, currentF = 0.7, weight = 0.5, plot=FALSE){
  MSYvec <- seq(from=MSYval*0.6, to=MSYval*1.4, length.out=1000)

  res <- matrix(NA, nrow=length(MSYvec), ncol=2)
  res[,1] <- MSYvec
  colnames(res) <- c("MSYvec", "NLL")

  for(i in 1:length(MSYvec)){
    x <- optim(par=rval,
               fn=MSYprofile.min,
               df=df,
               MSYval=MSYvec[i],
               method="Brent",
               lower=rval*0.5, upper=rval*1.5,
               OWT=FALSE, currentF = 0.7, weight = 1000)
    res[i,2] <- x$value # NLL
  }


  # calculate the best MSY and CI
  bestNLL <- res[[which.min(res[,2]),2]]
  bestMSY <- res[[which.min(res[,2]),1]]
  lowdf <- res[1:which.min(res[,2]),]
  uppdf <- res[which.min(res[,2]):1000,]

  # use an example from https://stackoverflow.com/questions/30314007/find-nearest-smaller-number
  maxless.low <- max(lowdf[,2][lowdf[,2] <= min(lowdf[,2])+1.92])
  lowCI <- lowdf[[which(lowdf[,2] == maxless.low), 1]]

  maxless.upp <- max(uppdf[,2][uppdf[,2] <= min(uppdf[,2])+1.92])
  uppCI <- uppdf[[which(uppdf[,2] == maxless.upp), 1]]

  finres <- data.frame(MSY=bestMSY, lowerCI=lowCI, upperCI=uppCI)

  if(plot==TRUE){
    plot(x=MSYvec, y=res[,2], type="l", xlab="MSY values", ylab="NLL", col="blue")
    abline(h=min(res[,2]), lty=2, col="red")
    abline(h=min(res[,2])+1.92, lty=2, col="gray") # CI95%
  }

  return(finres)
}

#' @title Schaefer's projection function
#'
#' @description
#' Function to create projection in Schaefer's model based on different reference points (MSY, Emsy) and
#' different Total Allowable Catch setting (TAC). This projection can be used to assist in the development of
#' limit and target reference point.
#'
#' @param inpars fitted surplus production parameters which consist of K (carrying capacity),
#' B0 (biomass when fishing is started), r (intrinsic growth rate), q (catchability coefficient), and
#' sigma (observation error)
#' @param df dataframe containing three columns; year, catch and unit of effort
#' @param nyears number of years the projection for the fishery
#' @param nsims number of iteration performed
#' @param TAC number to drive the management level
#' @param plot whether the projection plot is produced or not
#'
#' @return
#' The TAC is set at 1 in default. It would depend on the fishery to set whether the TAC can be set lower
#' to be more conservative or set higher to increase more catch.
#' For instance, when the TAC is set conservative at 0.8 of reference point, it will return
#' 0.8MSY, and 0.8Emsy. In contrary, the TAC will increase catch when it is set at 1.2
#' (or other value higher than 1) and return 1.2MSY, and 1.2Emsy.
#'
#' @export
#'
#' @examples
#'
#' fit <- calc.MSY(K=1000, B0=1000, r=0.2, q=0.00025, s.sigma=0.1,
#' df=df.goodcontrast, OWT=FALSE, currentF = 0.7, weight=0.5, plot=TRUE)
#'
#' run.Proj(inpars=fit[[1]][,2], df=df.goodcontrast, nyears=30, nsims=100, TAC=1, plot=TRUE)
#'
#'


run.Proj <- function(inpars, df,
                  nyears, nsims,
                  TAC = 1, plot = TRUE) {

  K <- inpars[1]
  B0 <- inpars[2]
  r <- inpars[3]
  q <- inpars[4]
  s.sigma <- inpars[5]

  Emsy <- r / (2 * q) #effort at MSY
  Bmsy <- K / 2 #biomass at MSY
  MSY <- (r * K) / 4 #MSY
  Fmsy <- r / 2 #fishing mortality at MSY

  EBt.msy <- EBt.Emsy  <- EBt.Bmsy <- vector(length = nrow(df))

  Year <- E.msy <- E.Emsy <- E.Bmsy <- C.msy <- C.Emsy <- C.Bmsy <- F.msy <- F.Emsy <- F.Bmsy <- vector(length = nrow(df) + nyears)
  Year <- df$year[1]:(df$year[1] + nrow(df) + nyears - 1)

  # observation df
  EBt.msy[1] <- EBt.Emsy[1] <- EBt.Bmsy[1] <- B0
  for (i in 1:(nrow(df))) {
    EBt.msy[i + 1] <- EBt.msy[i] + EBt.msy[i] * r * (1 - (EBt.msy[i] / K)) - df[i, 2]
    EBt.Emsy[i + 1] <- EBt.Emsy[i] + EBt.Emsy[i] * r * (1 - (EBt.Emsy[i] / K)) - df[i, 2]
    EBt.Bmsy[i + 1] <- EBt.Bmsy[i] + EBt.Bmsy[i] * r * (1 - (EBt.Bmsy[i] / K)) - df[i, 2]
  }

  E.msy[1:(nrow(df))] <- E.Emsy[1:(nrow(df))] <- E.Bmsy[1:(nrow(df))] <- df$catch / df$cpue
  C.msy[1:(nrow(df))] <- C.Emsy[1:(nrow(df))] <- C.Bmsy[1:(nrow(df))] <- df$catch

  for (i in 1:nrow(df)) {
    F.msy[i] <- df$catch[i] / ((EBt.msy[i] + EBt.msy[i + 1]) / 2)
    F.Emsy[i] <- df$catch[i] / ((EBt.Emsy[i] + EBt.Emsy[i + 1]) / 2)
  }

  # projection based on MSY measure
  matr.MSY <- replicate(n = nsims,
                        expr = {
                          for (i in (nrow(df) + 1):((nrow(df)) + nyears)) {
                            C.msy[i] <- TAC * MSY
                            E.msy[i] <- C.msy[i] / (q * EBt.msy[i])
                            EBt.msy[i + 1] <- (EBt.msy[i] + EBt.msy[i] * r * (1 - (EBt.msy[i] / K)) - C.msy[i]) * rlnorm(n = 1, meanlog = (-s.sigma ^ 2 / 2), sdlog = s.sigma)  # stochastic
                            F.msy[i] <- C.msy[i] / ((EBt.msy[i] + EBt.msy[i + 1]) / 2)
                          }

                          B_B.msy <- head(EBt.msy, -1) / Bmsy
                          F_F.msy <- F.msy / Fmsy

                          data.frame(Year = Year, Catch = C.msy, Effort = E.msy, EstBt = head(EBt.msy, -1), B_Bmsy = B_B.msy, F_Fmsy = F_F.msy)
                        },
                        simplify = F)

  temp.MSY <- matrix(NA, nrow=length(Year), ncol=16)
  temp.MSY[,1] <- Year

  ci.fun <- function(x){
    nSample <- length(x)
    meanSample <- mean(x) # mean from each row
    sdSample <- sd(x) # sd from each row

    marginError <- qt(0.975, df=nSample-1) * sdSample / sqrt(nSample)
    lci <- meanSample - marginError
    uci <- meanSample + marginError
    return(c(lci, uci))
  }

  for (rows in 1:length(Year)){
    temp.MSY[rows, 2] <- mean(sapply(matr.MSY, '[[', 2)[rows,]) #mu
    temp.MSY[rows, 3] <- mean(sapply(matr.MSY, '[[', 3)[rows,])
    temp.MSY[rows, 4] <- mean(sapply(matr.MSY, '[[', 4)[rows,])
    temp.MSY[rows, 5] <- mean(sapply(matr.MSY, '[[', 5)[rows,])
    temp.MSY[rows, 6] <- mean(sapply(matr.MSY, '[[', 6)[rows,])
    temp.MSY[rows, 7] <- ci.fun(sapply(matr.MSY, '[[', 2)[rows,])[1] #lci
    temp.MSY[rows, 8] <- ci.fun(sapply(matr.MSY, '[[', 3)[rows,])[1]
    temp.MSY[rows, 9] <- ci.fun(sapply(matr.MSY, '[[', 4)[rows,])[1]
    temp.MSY[rows, 10] <- ci.fun(sapply(matr.MSY, '[[', 5)[rows,])[1]
    temp.MSY[rows, 11] <- ci.fun(sapply(matr.MSY, '[[', 6)[rows,])[1]
    temp.MSY[rows, 12] <- ci.fun(sapply(matr.MSY, '[[', 2)[rows,])[2] #uci
    temp.MSY[rows, 13] <- ci.fun(sapply(matr.MSY, '[[', 3)[rows,])[2]
    temp.MSY[rows, 14] <- ci.fun(sapply(matr.MSY, '[[', 4)[rows,])[2]
    temp.MSY[rows, 15] <- ci.fun(sapply(matr.MSY, '[[', 5)[rows,])[2]
    temp.MSY[rows, 16] <- ci.fun(sapply(matr.MSY, '[[', 6)[rows,])[2]
  }
  colnames(temp.MSY) <- c("Year", "Catch.mu", "Effort.mu", "EstBt.mu", "B_Bmsy.mu", "F_Fmsy.mu",
                          "Catch.lci", "Effort.lci", "EstBt.lci", "B_Bmsy.lci", "F_Fmsy.lci",
                          "Catch.uci", "Effort.uci", "EstBt.uci", "B_Bmsy.uci", "F_Fmsy.uci")

  # projection based on Emsy measure
  matr.EMSY <- replicate(n = nsims,
                         expr = {
                           for (i in (nrow(df) + 1):((nrow(df)) + nyears)) {
                             E.Emsy[i] <- TAC * Emsy
                             C.Emsy[i] <- q * E.Emsy[i] * EBt.Emsy[i] * rlnorm(n = 1, meanlog = (-s.sigma ^ 2 / 2), sdlog = s.sigma)
                             EBt.Emsy[i + 1] <- (EBt.Emsy[i] + EBt.Emsy[i] * r * (1 - (EBt.Emsy[i] / K)) - C.Emsy[i]) * rlnorm(n = 1, meanlog = (-s.sigma ^ 2 / 2), sdlog = s.sigma)  # stochastic
                             F.Emsy[i] <- C.Emsy[i] / ((EBt.Emsy[i] + EBt.Emsy[i + 1]) / 2)
                           }

                           B_B.Emsy <- head(EBt.Emsy, -1) / Bmsy
                           F_F.Emsy <- F.Emsy / Fmsy

                           data.frame(Year = Year, Catch = C.Emsy, Effort = E.Emsy, EstBt = head(EBt.Emsy, -1), B_BEmsy = B_B.Emsy, F_FEmsy = F_F.Emsy)
                         },
                         simplify = F)

  temp.EMSY <- matrix(NA, nrow=length(Year), ncol=16)
  temp.EMSY[,1] <- Year

  for (rows in 1:length(Year)){
    temp.EMSY[rows, 2] <- mean(sapply(matr.EMSY, '[[', 2)[rows,]) #mu
    temp.EMSY[rows, 3] <- mean(sapply(matr.EMSY, '[[', 3)[rows,])
    temp.EMSY[rows, 4] <- mean(sapply(matr.EMSY, '[[', 4)[rows,])
    temp.EMSY[rows, 5] <- mean(sapply(matr.EMSY, '[[', 5)[rows,])
    temp.EMSY[rows, 6] <- mean(sapply(matr.EMSY, '[[', 6)[rows,])
    temp.EMSY[rows, 7] <- ci.fun(sapply(matr.EMSY, '[[', 2)[rows,])[1] #lci
    temp.EMSY[rows, 8] <- ci.fun(sapply(matr.EMSY, '[[', 3)[rows,])[1]
    temp.EMSY[rows, 9] <- ci.fun(sapply(matr.EMSY, '[[', 4)[rows,])[1]
    temp.EMSY[rows, 10] <- ci.fun(sapply(matr.EMSY, '[[', 5)[rows,])[1]
    temp.EMSY[rows, 11] <- ci.fun(sapply(matr.EMSY, '[[', 6)[rows,])[1]
    temp.EMSY[rows, 12] <- ci.fun(sapply(matr.EMSY, '[[', 2)[rows,])[2] #uci
    temp.EMSY[rows, 13] <- ci.fun(sapply(matr.EMSY, '[[', 3)[rows,])[2]
    temp.EMSY[rows, 14] <- ci.fun(sapply(matr.EMSY, '[[', 4)[rows,])[2]
    temp.EMSY[rows, 15] <- ci.fun(sapply(matr.EMSY, '[[', 5)[rows,])[2]
    temp.EMSY[rows, 16] <- ci.fun(sapply(matr.EMSY, '[[', 6)[rows,])[2]
  }
  colnames(temp.EMSY) <- c("Year", "Catch.mu", "Effort.mu", "EstBt.mu", "B_BEmsy.mu", "F_FEmsy.mu",
                           "Catch.lci", "Effort.lci", "EstBt.lci", "B_BEmsy.lci", "F_FEmsy.lci",
                           "Catch.uci", "Effort.uci", "EstBt.uci", "B_BEmsy.uci", "F_FEmsy.uci")

  res <- list(MSY=temp.MSY, Emsy=temp.EMSY)

  if (plot == TRUE){
    par(mfrow=c(2,2), mar=c(0,0,0,0), oma=c(4,4,3,1))

    #MSY
    plot(x=res[[1]][,1], y=res[[1]][,5], type="l", col="blue", ylim=c(0,2.4), xaxt="n")
    lines(x=res[[1]][,1], y=res[[1]][,10], type="l", lty=2, col="red")
    lines(x=res[[1]][,1], y=res[[1]][,15], type="l", lty=2, col="red")
    lines(x=res[[1]][1:nrow(df),1], y=res[[1]][1:nrow(df),10], type="l", col="black")
    mtext(text="B/Bmsy", side=2, line=2.6)
    mtext(text=paste("Trajectory based on", TAC, "* MSY"), side=3, line=1.3)

    #EMSY
    plot(x=res[[2]][,1], y=res[[2]][,5], type="l", col="blue", ylim=c(0,2.4), xaxt="n", yaxt="n")
    lines(x=res[[2]][,1], y=res[[2]][,10], type="l", lty=2, col="red")
    lines(x=res[[2]][,1], y=res[[2]][,15], type="l", lty=2, col="red")
    lines(x=res[[2]][1:nrow(df),1], y=res[[2]][1:nrow(df),10], type="l", col="black")
    mtext(text=paste("Trajectory based on", TAC, "* Emsy"), side=3, line=1.3)
    legend("bottomright", legend=c("Mean projection", "95% CI", "Observation"), col=c("blue", "red", "black"), lty=c(1,2,1), bty="n")

    #MSY
    plot(x=res[[1]][,1], y=res[[1]][,6], type="l", col="blue", ylim=c(0,2.4))
    lines(x=res[[1]][,1], y=res[[1]][,11], type="l", lty=2, col="red")
    lines(x=res[[1]][,1], y=res[[1]][,16], type="l", lty=2, col="red")
    lines(x=res[[1]][1:nrow(df),1], y=res[[1]][1:nrow(df),11], type="l", col="black")
    mtext(text="F/Fmsy", side=2, line=2.6)

    #EMSY
    plot(x=res[[2]][,1], y=res[[2]][,6], type="l", col="blue", ylim=c(0,2.4), yaxt="n")
    lines(x=res[[2]][,1], y=res[[2]][,11], type="l", lty=2, col="red")
    lines(x=res[[2]][,1], y=res[[2]][,16], type="l", lty=2, col="red")
    lines(x=res[[2]][1:nrow(df),1], y=res[[2]][1:nrow(df),11], type="l", col="black")
    mtext(text="Years", side=1, line=2.6, outer=T)
    legend("bottomright", legend=c("Mean projection", "95% CI", "Observation"), col=c("blue", "red", "black"), lty=c(1,2,1), bty="n")
  }
  return(res)
}

