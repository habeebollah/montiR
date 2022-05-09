#' @title Schaefer and Fox' parameter estimation using bayesian
#'
#' @description
#' Function to estimate the Schaefer and Fox' parameter using bayesian approach. This approach requires r from the assessed species as input prior in the bayesian,
#' to increase the accuracy of surplus production parameter estimation in Schaefer and Fox model.
#'
#' This function requires JABBA package installed.
#'
#' @param df dataframe containing three columns; year, catch and effort
#' @param K.prior carrying capacity (K) parameter prior. The default is NULL, can also receive input (mu, cv)
#' @param psi.prior ratio between B0/K parameter prior. The default is (0.9, 0.25), can also receive input (mu, cv)
#' @param r.prior intrinsic growth rate (r) parameter prior. The default is (mu, cv) with mu from rGrabber median parameter and cv is set at 0.5
#' @param SPmodel option on Surplus Production model; 1 for Schaefer and 2 for Fox
#'
#' @return
#' @export
#'
#' @references
#' Winker, H., Carvalho, F., Kapur, M. (2018) JABBA: Just Another Bayesian Biomass Assessment. Fisheries Research 204: 275-288.
#'
#' @examples
#' # grab r intrinsic growth parameter and use it for the next analysis
#' r.val <- rGrabber("Hoplostethus atlanticus")
#'
#' library("JABBA")
#' SFpar_bayes(df=df.eastpacCatch, K.prior, psi.prior, r.prior=r.val, SPmodel=1)
#'

SFpar_bayes <- function(df, K.prior=NULL, psi.prior=c(0.9, 0.25), r.prior, SPmodel){
  catch <- df[,c(1,2)]
  cpue <- data.frame(year=df[,1],cpue=df[,2]/df[,3])
  K.prior <- K.prior
  psi.prior <- psi.prior
  r.mean <- r.prior$median_r

  temp <- ifelse(SPmodel==1, type <- "Schaefer",
         ifelse(SPmodel==2, type <- "Fox", type <- print("wrong SPmodel code")))

  jbinput <- build_jabba(catch=df[,c(1,2)], cpue=data.frame(year=df[,1],cpue=df[,2]/df[,3]),
                         K.dist = c("lnorm", "range"), K.prior = K.prior,
                         psi.dist = c("lnorm", "range"), psi.prior = psi.prior,
                         r.dist = c("lnorm", "range"), r.prior=c(r.mean, 0.5), # sd is estimated at 0.5 since there is no way to calculate this value from provided mean and 95% CI
                         model.type = type)
  res.bayes <- fit_jabba(jbinput, init.r=r.mean, quickmcmc=TRUE)
  res <- data.frame(K=res.bayes[["pars"]]$Median[1],
                    B0=res.bayes[["pars"]]$Median[4]*res.bayes[["pars"]]$Median[1],
                    r=res.bayes[["pars"]]$Median[2],
                    q=res.bayes[["pars"]]$Median[3])
  return(res)
}
