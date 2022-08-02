#' @title Fox's data fitting and plotting
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
#' @return
#' @export
#'
#' @references
#' Hilborn, Ray, and Carl J. Walters, eds. Quantitative fisheries stock assessment: choice,
#' dynamics and uncertainty. Springer Science & Business Media, 1992.
#'
#' Bonfil, R. 2005. ‘Fishery stock assessment models and their application to sharks’, in Musick, J. A. and Bonfil,
#' Ramon (eds) Wildlife Conservation. Rome: Food and Agriculture Organizations of The United Nations, pp. 154–181.
#'
#' @examples
#' K <- 1000
#' B0 <- K
#' r <- 0.2
#' q <- 0.00025
#'
#' inpars <- c(log(K), log(B0), log(r), log(q))
#' Fpar(inpars=inpars, df=df.goodcontrast)
#'
Fpar <- function(inpars, df){
  K <- inpars[1]
  B0 <- inpars[2]
  r <- inpars[3]
  q <- inpars[4]

  CPUE <- df$catch/df$effort
  EstCatch <- EstCPUE <- EstBt <- vector(length=nrow(df))

  for (i in 1:nrow(df)) {
    if (i == 1) EstBt[i] <- B0
    if (i>1) EstBt[i] <- EstBt[i-1] +  EstBt[i-1] * r * (1 - log(EstBt[i-1])/log(K)) - df$catch[i-1]
    EstCatch[i] <- EstBt[i] * q * df$effort[i]
  }

  EstCPUE <- EstCatch / df$effort

  plot(CPUE, xlab="Year", ylab="CPUE", type="b", col="blue")
  lines(EstCPUE, type="b", col="red")
  legend("topright", legend=c("Observation", "Estimation"), lty=c(1, 1), col=c("blue", "red"), box.lty=1, cex=0.7)

  return(data.frame(CPUE=CPUE, EstBt=EstBt, EstCatch=EstCatch, EstCPUE=EstCPUE))
}
