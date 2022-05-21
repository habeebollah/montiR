#' @title Plotting catch, effort and CPUE
#'
#' @description
#' This function plot the catch, effort and CPUE from a dataframe to check,
#' whether the plot shows one-way trip pattern or show a good contrast as required for the biomass dynamic model. It also
#' gives a histogram on the data, and compare it with the normal and lognormal distribution which
#' will be important on choosing the type of distribution used in the data fitting process.
#'
#' @param df dataframe containing three columns; year, catch and effort
#'
#' @return
#' @export
#'
#' @references
#' Hilborn, Ray, and Carl J. Walters, eds. Quantitative fisheries stock assessment: choice,
#' dynamics and uncertainty. Springer Science & Business Media, 1992.
#'
#' Magnusson, Arni, and Ray Hilborn. "What makes fisheries data informative?." Fish and Fisheries 8, no. 4 (2007): 337-358.
#'
#' @examples
#' plotInit(df.goodcontrast)
#' plotInit(df.onewaytrip)
#'
plotInit <- function(df){
  I <- df[,2]/df[,3]
  lmI <- lm(I ~ df[,1])
  I2 <- seq(min(I), max(I), length=100)
  nline <- dnorm(I2, mean(I), sd(I))
  lnline <- dlnorm(I2, mean(log(I)), sd(log(I)))

  par(mfrow=c(2,2), mar=c(3,4,1,2))
  plot(x=df[,1], y=df[,2], type="l", ylab="Catch")

  plot(x=df[,1], y=df[,3], type="l", ylab="Effort")

  #plot(x=df[,1], y=I, type="l", ylab="CPUE")
  plot(I ~ df[,1], type="l", ylab="CPUE", col="black")
  abline(lmI$coefficients[1], b=lmI$coefficients[2], lty=1, col="red")
  legend("topright", c("This data", "One-way trip"), bty = "n",
         cex=0.8, col = c("black", "red"), lwd = c(1, 1))

  hist(I, prob=T, freq=F, main=NULL)
  lines(density(I), col="black")
  lines(I2, nline, col="red")
  lines(I2, lnline, col="blue")
  legend("topright", c("This data", "Normal", "LogNormal"), bty = "n",
         cex=0.8, col = c("black", "red", "blue"), lwd = c(1, 1, 1))
}
