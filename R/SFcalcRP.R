#' @title Schaefer and Fox' reference point calculator from data fitting parameters
#'
#' @description
#' This function is used to calculate the reference points (Bmsy, MSY, Emsy) under Schaefer and Fox model
#'
#' @param inpars surplus production parameters which consist of K (carrying capacity), B0 (biomass when fishing is started), r (intrinsic growth rate), q (catchability coefficient)
#' @param SPmodel option on Surplus Production model; 1 for Schaefer and 2 for Fox
#'
#' @return
#' @export
#'
#' @references
#' Hilborn, Ray, and Carl J. Walters, eds. Quantitative fisheries stock assessment: choice, dynamics and uncertainty. Springer Science & Business Media, 1992.
#'
#' Bonfil, R. (2005) ‘Fishery stock assessment models and their application to sharks’, in Musick, J. A. and Bonfil, Ramon (eds) Wildlife Conservation. Rome: Food and Agriculture Organizations of The United Nations, pp. 154–181.
#'
#' @examples
#' K <- 1000
#' B0 <- K
#' r <- 0.2
#' q <- 0.00025
#'
#' inpars <- c(K, B0, r, q)
#' SFcalcRP(inpars=inpars, SPmodel=2)
#'

SFcalcRP <- function(inpars, SPmodel){
  K <- inpars[1]
  B0 <- inpars[2]
  r <- inpars[3]
  q <- inpars[4]

  ifelse(SPmodel==1, res <- list(Bmsy = K/2, MSY = (r*K)/4, Emsy = r/(2*q)),
         ifelse(SPmodel==2, res <- list(Bmsy = K/exp(1), MSY = (r*K)/(exp(1)*log(K)), Emsy = r/q), res <- "wrong SPmodel code!"))
return(res)
}
