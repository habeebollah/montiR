#' @title Schaefer's multiple regression method
#'
#' @description This function calculate MSY and Emsy, as well as the r squared and adj R squared under non-equilibrium assumption.
#' The least square in the Walters and Hilborn (1976)'s multiple regression produces bias parameter (Uhler, 1979), which not advised to be used as it would result in further bias when calculating MSY and Emsy.
#' Based on the input, Hilborn and Walters (1992) adjusted the equation as input in the multiple regression.
#'
#' @param df dataframe containing three columns; year, catch and effort
#'
#' @return
#' @export
#'
#' @references
#' Carl J. Walters and Ray Hilborn. 1976. Adaptive Control of Fishing Systems. Journal of the Fisheries Research Board of Canada. 33(1): 145-159
#'
#' Uhler, R. S. 1979. Least squares regression estimates of the Schaefer production model: some Monte Carlo simulation results. Can. J. Fish. Aquat. Sci. 37: 1284-1294.
#'
#' Hilborn, Ray, and Carl J. Walters, eds. Quantitative fisheries stock assessment: choice, dynamics and uncertainty. Springer Science & Business Media, 1992.
#'
#' @examples
#' Smreg(df.eastpacCatch)
#'
Smreg <- function(df){

  df$Ut <- df$catch/df$effort
  for (i in 1:nrow(df)){
    df$Ut1[i] <- (df$Ut[i+1]/df$Ut[i])-1
  }

  # Walters and Hilborn (1976)
  mreg1 <- lm(df$Ut1 ~ df$Ut + df$effort)

  r1 <- mreg1$coefficients[[1]]
  q1 <- -mreg1$coefficients[[3]]
  K1 <- -r1/(mreg1$coefficients[[2]]*q1)

  MSY1 <- (r1*K1)/4
  Emsy1 <- r1/(2*q1)

  # Hilborn and Walters (1992)
  mr2a <- df$Ut1
  mr2b <- df$Ut
  mr2c <- df$Ut^2
  mr2d <- (df$effort*df$Ut)

  mreg2 <- lm(mr2a ~ mr2b + mr2c + mr2d)
  summary(mreg2)$r.squared
  r2 <- mreg2$coefficients[[2]]
  q2 <- -mreg2$coefficients[[4]]
  K2 <- -r2/(mreg2$coefficients[[3]]*q2)

  MSY2 <- (r2*K2)/4
  Emsy2 <- r2/(2*q2)

  res <- data.frame(analysis=c("K", "r", "q", "MSY", "Emsy", "r2", "adj.r2"),
                    WH.1976=c(K1, r1, q1, MSY1, Emsy1, summary(mreg1)$r.squared, summary(mreg1)$adj.r.squared),
                    HW.1992=c(K2, r2, q2, MSY2, Emsy2, summary(mreg2)$r.squared, summary(mreg2)$adj.r.squared))
  return(res)
}

