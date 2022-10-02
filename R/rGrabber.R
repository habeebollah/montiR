#' @title Grab intrinsic growth rate (r) parameter from Fishbase and SeaLifebase
#'
#' @description
#' Function to download intrinsic growth rate (r) parameter for specific species from Fishbase
#' and SeaLifebase. The r growth parameter will serve as prior in the bayesian analysis in the surplus production model.
#'
#' @param SciName insert scientific name for the fish or aquatic biota
#'
#' @return
#' The r value is ranging from 0 to 1, where high r value represents the r strategist species and
#' known to have more resilient characteristics when facing with fishing pressures (Adams, 1980; Kawasaki, 1980, 1983).
#'
#' This function requires rfishbase package installed and internet connection.
#'
#' @export
#'
#' @references
#' Adams, P. (1980) ‘Life history patterns in marine fishes and their consequences for fisheries management’,
#' Fishery Bulletin, 78(1), pp. 1–12.
#'
#' Kawasaki, T. (1980) ‘Fundamental Relations among the Selections of Life History in the Marine Teleosts’,
#' Bulletin of the Japanese Society of Scientific Fisheriess, 463(3), pp. 289–293.
#'
#' Kawasaki, T. (1983) Why do some pelagic fishes have wide fluctuations in their numbers? Biological basis of
#' fluctuation from the viewpoint of evolutionary ecology, FAO Fisheries Report (FAO).
#'
#' https://www.fishbase.in/manual/key%20facts.htm
#'
#' @examples
#' library(rfishbase)
#' #grab r intrinsic growth rate parameter from the database
#' # data grabbed from fishbase
#' rGrabber("Sardinella lemuru") #for finfish species with a high r value
#' rGrabber("Hoplostethus atlanticus") #for finfish species with a low r value
#'
#' # data grabbed from sealifebase
#' rGrabber("Anguilla japonica") # for non finfish species
#'
rGrabber <- function(SciName){
  notes1 <- "These values are generated from the resilience information. The method gives more reliable estimation on r growth parameter (R. Froese pers.comm., March 2022)"
  notes2 <- "These values are generated from stock assessment models. The method is in the evaluation process, sometimes it overestimates on r growth parameter (R. Froese pers.comm., March 2022)"
  notes3 <- "r growth parameter or resilience information for this species has not been populated to this database yet.
  If the information is not available when you open www.fishbase.org or www.sealifebase.org under Estimates of some properties based on models, then you might want to check from another sources"
  if(is.na(estimate(SciName)[[15]])){
    # estimates use resilience value from table 2. https://www.fishbase.de/rfroese/CMSYUserGuideMarch2021.pdf
    temp <- ifelse(is.na(stocks(SciName)[[38]]), res <- notes3,
                   ifelse(stocks(SciName)[[38]]=="High",
                          res <- list(median_r=median(c(0.6,1.5)), lci_r=0.6, uci_r=1.5, resilience="High", notes=notes1),
                          ifelse(stocks(SciName)[[38]]=="Medium",
                                 res <- list(median_r=median(c(0.2,0.8)), lci_r=0.2, uci_r=0.8, resilience="Medium", notes=notes1),
                                 ifelse(stocks(SciName)[[38]]=="Low",
                                        res <- list(median_r=median(c(0.005,0.5)), lci_r=0.005, uci_r=0.5, resilience="Low", notes=notes1),
                                        res <- list(median_r=median(c(0.015,0.1)), lci_r=0.015, uci_r=0.1, resilience="Very low", notes=notes1)))))
  }
  else{
    # directly accessing the fishbase and sealifebase data
    res <- list(median_r=estimate(SciName)[[15]],
                lci_r=estimate(SciName)[[16]],
                ucl_r=estimate(SciName)[[17]],
                resilience=stocks(SciName)[[38]][[1]],
                notes=notes2)
  }
  return(res)
}
