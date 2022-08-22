#' @title Schaefer's data fitting and plotting
#'
#' @description
#' This function is used to predict the initial Biomass Dynamic Model parameters as input to do
#' the estimation in the latter process.
#'
#' It also works to check the result of estimated parameters by plotting the data, whether the
#' data and estimated parameter is fitted.
#'
#' @param inpars surplus production parameters which consist of K (carrying capacity), B0
#' (biomass when fishing is started), r (intrinsic growth rate), q (catchability coefficient)
#' @param df dataframe containing three columns; year, catch and effort
#'
#' @export
#'
#' @references
#' Hilborn, Ray, and Carl J. Walters, eds. Quantitative fisheries stock assessment: choice,
#' dynamics and uncertainty. Springer Science & Business Media, 1992.
#'
#' @examples
#' K <- 1000
#' B0 <- K
#' r <- 0.2
#' q <- 0.00025
#'
#' inpars <- c(K, B0, r, q)
#' Par_init(inpars=inpars, df=df.goodcontrast)
#'
Par_init <- function(inpars, df){
  K <- inpars[1]
  B0 <- inpars[2]
  r <- inpars[3]
  q <- inpars[4]

  CPUE <- df$catch/df$effort
  EstCatch <- EstCPUE <- EstBt <- vector(length=nrow(df))

  for (i in 1:nrow(df)) {
    if (i == 1) EstBt[i] <- B0
    if (i>1) EstBt[i] <- EstBt[i-1] +  EstBt[i-1] * r * (1 - EstBt[i-1]/K) - df$catch[i-1]
    EstCatch[i] <- EstBt[i] * q * df$effort[i]
  }

  EstCPUE <- EstCatch / df$effort

  plot(CPUE, xlab="Year", ylab="CPUE", type="b", col="blue")
  lines(EstCPUE, type="b", col="red")
  legend("topright", legend=c("Observation", "Estimation"), lty=c(1, 1), col=c("blue", "red"), box.lty=1, cex=0.7, bty="n")

  return(data.frame(CPUE=CPUE, EstBt=EstBt, EstCatch=EstCatch, EstCPUE=EstCPUE))
}


#' @title Schaefer's parameter estimation minimizing function
#'
#' @description
#' Function used in the maximum likelihood minimization to estimate Schaefer's parameters.
#' Observation error is used to increase the accuracy of data fitting. It is assumed to occur in the relationship
#' between stock biomass and index of abundance and is estimated assuming lognormal distribution in maximum likelihood (Polacheck et al., 1993).
#'
#' Since fishing effort data collection are not always conducted regularly while catch is likely
#' have a better time series information, this function also allow for some lose of data.
#'
#' This function also consider the different quality data, for instance if the data shows
#' a one way trip pattern which losing increasing rate of increase.
#'
#' @param inpars surplus production parameters which consist of K (carrying capacity), B0
#' (biomass when fishing is started), r (intrinsic growth rate), q (catchability coefficient),
#' s.sigma (observation error)
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
#' s.sigma <- 0.1
#'
#' ### Estimate parameters using optim
#' inpars <- c(log(K), log(B0), log(r), log(q), log(s.sigma))
#'
#' fit <- optim(par=inpars,
#'              fn=Par.min,
#'              df=df.goodcontrast,
#'              method="Nelder-Mead",
#'              OWT=FALSE, currentF = 0.7, weight = 0.5)
#'
#' Par.vals <- data.frame(SPpar = c("K", "B0", "r", "q", "s.sigma"),
#'                          init_pars = c(K, B0, r, q, s.sigma),
#'                          fitted_pars = exp(fit$par))
#'
Par.min <- function(inpars, df, OWT=FALSE, currentF = 0.7, weight = 1000){
  K <- exp(inpars[1])
  B0 <- exp(inpars[2])
  r <- exp(inpars[3])
  q <- exp(inpars[4])
  s.sigma <- exp(inpars[5])

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
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = s.sigma, log = TRUE))
    }

  if (OWT=="Depletion"){
    annualFrates <- df[,2]/EstBt # F = catch/biomass = Z-M
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = s.sigma, log = TRUE)) +
      weight * (tail(annualFrates,1) - currentF)^2
    #print(tail(annualFrates,1)) # to check whether the weighting makes the estimates depletion getting closer to the current exploitation rate
  }

  if (OWT=="Biomass"){
    surveyB <- df[,4]
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = s.sigma, log = TRUE)) +
      -sum(weight * (surveyB - EstBt)^2)
  }

  return(nll)
}

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
#' and s.sigma (observation error)
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
#' s.sigma <- 0.1
#'
#' Bmsy <- K/2
#' MSY <- (r*K)/4
#' Emsy <- r/(2*q)
#'
#' ### Estimate parameters using optim
#' inpars <- c(Bmsy, MSY, Emsy, B0, s.sigma)
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
#' ParRP_vals <- data.frame(RPpar = c("Bmsy", "MSY", "Emsy", "B0", "s.sigma"),
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
  s.sigma <- inpars[5]

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
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = s.sigma, log = TRUE))
  }

  if (OWT=="Depletion"){
    annualFrates <- df[,2]/EstBt # F = catch/biomass = Z-M
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = s.sigma, log = TRUE)) +
      weight * (tail(annualFrates,1) - currentF)^2
    #print(tail(annualFrates,1)) # to check whether the weighting makes the estimates depletion getting closer to the current exploitation rate
  }

  if (OWT=="Biomass"){ # haven't been checked yet
    surveyB <- df[,4]
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = s.sigma, log = TRUE)) +
      -sum(weight * (surveyB - EstBt)^2)
  }

  return(nll)
}


#' @title Schaefer's projection function
#'
#' @description
#' Function to create projection in Schaefer's model based on different reference points (MSY, Emsy) and
#' different Total Allowable Catch setting (TAC). This projection can be used to assist in the development of
#' limit and target reference point.
#'
#' @param inpars fitted surplus production parameters which consist of K (carrying capacity),
#' B0 (biomass when fishing is started), r (intrinsic growth rate), q (catchability coefficient)
#' @param df dataframe containing three columns; year, catch and effort
#' @param nyears number of years the projection for the fishery
#' @param nsims number of iteration performed
#' @param TAC number to drive the management level
#' @param s.sigma standard deviation estimated in the surplus production data fitting process
#' @param graph whether the plot is produced or not
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
#' K <- 1000
#' B0 <- K
#' r <- 0.2
#' q <- 0.00025
#' s.sigma <- 0.1
#'
#' ### Estimate parameters using optim
#' inpars <- c(log(K), log(B0), log(r), log(q), log(s.sigma))
#'
#' fit <- optim(par=inpars,
#'              fn=Par.min,
#'              df=df.goodcontrast,
#'              method="Nelder-Mead",
#'              OWT=FALSE, currentF = 0.7, weight = 0.5)
#'
#' fitted_pars <- exp(fit$par)
#'
#' Sproj(inpars=fitted_pars, df=df.goodcontrast, nyears=30, nsims=100, TAC=1, graph=TRUE)
#'
#'


Sproj <- function(inpars, df,
                  nyears, nsims,
                  TAC = 1, graph = TRUE) {

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

  E.msy[1:(nrow(df))] <- E.Emsy[1:(nrow(df))] <- E.Bmsy[1:(nrow(df))] <- df$effort
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

  if (graph == TRUE){
    par(mfrow=c(2,2), mar=c(0,0,0,0), oma=c(5,5,1,1))

    #MSY
    plot(x=res[[1]][,1], y=res[[1]][,5], type="l", col="blue", ylim=c(0,2.4), xaxt="n")
    lines(x=res[[1]][,1], y=res[[1]][,10], type="l", lty=2, col="red")
    lines(x=res[[1]][,1], y=res[[1]][,15], type="l", lty=2, col="red")
    lines(x=res[[1]][1:nrow(df),1], y=res[[1]][1:nrow(df),10], type="l", col="black")
    mtext(text="B/Bmsy", side=2, line=3)
    mtext(text=paste("Trajectory based on", TAC, "* MSY"), side=3, line=-1.5)

    #EMSY
    plot(x=res[[2]][,1], y=res[[2]][,5], type="l", col="blue", ylim=c(0,2.4), xaxt="n", yaxt="n")
    lines(x=res[[2]][,1], y=res[[2]][,10], type="l", lty=2, col="red")
    lines(x=res[[2]][,1], y=res[[2]][,15], type="l", lty=2, col="red")
    lines(x=res[[2]][1:nrow(df),1], y=res[[2]][1:nrow(df),10], type="l", col="black")
    mtext(text=paste("Trajectory based on", TAC, "* Emsy"), side=3, line=-1.5)
    legend("bottomright", legend=c("Mean projection", "95% CI", "Observation"), col=c("blue", "red", "black"), lty=c(1,2,1), bty="n")

    #MSY
    plot(x=res[[1]][,1], y=res[[1]][,6], type="l", col="blue", ylim=c(0,2.4))
    lines(x=res[[1]][,1], y=res[[1]][,11], type="l", lty=2, col="red")
    lines(x=res[[1]][,1], y=res[[1]][,16], type="l", lty=2, col="red")
    lines(x=res[[1]][1:nrow(df),1], y=res[[1]][1:nrow(df),11], type="l", col="black")
    mtext(text="F/Fmsy", side=2, line=3)

    #EMSY
    plot(x=res[[2]][,1], y=res[[2]][,6], type="l", col="blue", ylim=c(0,2.4), yaxt="n")
    lines(x=res[[2]][,1], y=res[[2]][,11], type="l", lty=2, col="red")
    lines(x=res[[2]][,1], y=res[[2]][,16], type="l", lty=2, col="red")
    lines(x=res[[2]][1:nrow(df),1], y=res[[2]][1:nrow(df),11], type="l", col="black")
    mtext(text="Years", side=1, line=3, outer=T)
    legend("bottomright", legend=c("Mean projection", "95% CI", "Observation"), col=c("blue", "red", "black"), lty=c(1,2,1), bty="n")
  }
  return(res)
}

