#' Surplus production under equilibrium assumption for Schaefer and Fox
#'
#' @param df dataframe containing three columns; year, catch and effort
#'
#' @return
#' This function calculate MSY and fMSY including its confident interval,
#' as well as the r squared and adj R squared, under equilibrium assumption using Schaefer and Fox models
#'
#' @export
#'
#'
#' @examples
#' SP_Eq(goodcontrast)
#'
SP_Eq <- function(df){

  # Schaefer
  CPUE.s <- df[,2]/df[,3]
  LReg.s <- lm(CPUE.s ~ df[,3])

  a.intercept.s <- summary(LReg.s)$coefficient[1]
  b.intercept.s <- summary(LReg.s)$coefficient[2]

  MSY.s <- -0.25 * a.intercept.s^2/b.intercept.s
  fMSY.s <- -0.5 * a.intercept.s/b.intercept.s

  # 2.5 and 97.5% CI is calculated from the slope
  MSY.s_CI <- c(-0.25 * a.intercept.s^2/confint(LReg.s)[2],
                -0.25 * a.intercept.s^2/confint(LReg.s)[4])
  fMSY.s_CI <- c(-0.5 * a.intercept.s/confint(LReg.s)[2],
                 -0.5 * a.intercept.s/confint(LReg.s)[4])

  # Catch[i]/effort[i] = q*B = a+b*effort[i] when f=0 reaches virgin biomass, hence
  # VirginBiomass.s <- a.intercept.s/q

  # Fox
  CPUE.f <- log(df[,2]/df[,3])
  LReg.f <- lm(CPUE.f ~ df[,3])

  a.intercept.f <- summary(LReg.f)$coefficient[1]
  b.intercept.f <- summary(LReg.f)$coefficient[2]

  MSY.f <- -(1/b.intercept.f)*exp(a.intercept.f-1)
  fMSY.f <- -1/b.intercept.f

  # 2.5 and 97.5% CI is calculated from the slope
  MSY.f_CI <- c(-(1/confint(LReg.f)[2])*exp(a.intercept.f-1),
                -(1/confint(LReg.f)[4])*exp(a.intercept.f-1))
  fMSY.f_CI <- c(-1/confint(LReg.f)[2], -1/confint(LReg.f)[4])

  # Catch[i]/effort[i] = q*B = exp(a+b*effort[i]) when f=0 reaches virgin biomass, hence
  # VirginBiomass.f <- exp(a.intercept.f)/q

  res <- list(data=cbind(df, CPUE.s, CPUE.f),
              result=data.frame(analysis=c("MSY", "MSY.CIlow", "MSY.CIupp", "Emsy", "Emsy.CIlow", "Emsy.CIupp", "r2", "adj.r2"),
                                Schaefer=c(MSY.s, MSY.s_CI[1], MSY.s_CI[2], Emsy.s, Emsy.s_CI[1], Emsy.s_CI[2], summary(LReg.s)$r.squared, summary(LReg.s)$adj.r.squared),
                                Fox=c(MSY.f, MSY.f_CI[1], MSY.f_CI[2], Emsy.f, Emsy.f_CI[1], Emsy.f_CI[2], summary(LReg.f)$r.squared, summary(LReg.f)$adj.r.squared)
              )
  )
  return(res)
}
