#' @title Grab intrinsic growth rate (r) parameter from Fishbase and SeaLifebase
#'
#' @description
#' Function to download intrinsic growth rate (r) parameter for specific species from Fishbase
#' and SeaLifebase. The r growth parameter will serve as prior in the bayesian analysis in the surplus production model.
#'
#' This function requires rfishpackage installed and internet connection.
#'
#' @param SciName insert scientific name for the fish or aquatic biota
#'
#' @return
#' @export
#'
#' @references
#' https://www.fishbase.in/manual/key%20facts.htm
#'
#' @examples
#' #grab r intrinsic growth rate parameter from the database
#' library("rfishbase")
#' rGrabber("Sardinella lemuru") #for finfish species with a high r value, data grabbed from fishbase
#' rGrabber("Hoplostethus atlanticus") #for finfish species with a low r value, data grabbed from fishbase
#' rGrabber("Anguilla japonica") # for non finfish species, data grabbed from sealifebase
#'


rGrabber <- function(SciName){
  notes1 <- "These values are generated from the resilience information.
  The method gives more reliable estimation on r growth parameter (R. Froese pers.comm.)"
  notes2 <- "These values are generated from stock assessment models.
  The method is in the evaluation process, sometimes it overestimates on r growth parameter (R. Froese pers.comm.)"
  notes3 <- "r growth parameter or resilience information for this species has not been populated to this database yet.
  if the information is not available when you open www.fishbase.org or www.sealifebase.org
  under Estimates of some properties based on models, then you might want to check from another sources"
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
