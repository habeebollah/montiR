#' @title Plotting catch, effort and CPUE
#'
#' @description
#' This function plot the catch, effort and CPUE from a dataframe to check,
#' whether the plot shows one-way trip pattern or show a good contrast as required for the biomass dynamic model.
#'
#' @param df dataframe containing three columns; year, catch and effort
#'
#' @return
#' @export
#'
#' @references
#' Hilborn, Ray, and Carl J. Walters, eds. Quantitative fisheries stock assessment: choice, dynamics and uncertainty. Springer Science & Business Media, 1992.
#'
#' Magnusson, Arni, and Ray Hilborn. "What makes fisheries data informative?." Fish and Fisheries 8, no. 4 (2007): 337-358.
#'
#' @examples
#' plotInit(df.goodcontrast)
#' plotInit(df.onewaytrip)
#'
plotInit <- function(df){
  par(mfrow=c(3,1), mar=c(3,4,1,2))
  plot(x=df[,1], y=df[,2], type="l", ylab="Catch")
  plot(x=df[,1], y=df[,3], type="l", ylab="Effort")
  plot(x=df[,1], y=df[,2]/df[,3], type="l", ylab="CPUE")
}
