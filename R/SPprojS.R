#' @title Schaefer's projection function
#'
#' @description
#' Function to create projection in Schaefer's model based on different reference points.
#'
#' @param inpars fitted surplus production parameters which consist of K (carrying capacity), B0 (biomass when fishing is started), r (intrinsic growth rate), q (catchability coefficient)
#' @param df dataframe containing three columns; year, catch and effort
#' @param nyears number of years the projection for the fishery
#' @param sigma standard deviation estimated in the surplus production data fitting process
#'
#' @return
#' sigma in default is set at a very low value (0.000001) to allow for a deterministic projection. In interest to create stochastic projection, the sigma value produced during maximum likelihood estimation or use another value can be utilized
#'
#' The TAC (Total Allowable Catch) is set at 1 in default. It would depend on the fishery to set whether the TAC can be set lower to be more conservative or set higher to increase more catch.
#'
#' @export
#'
#' @examples
#' SPprojS(inpars=fitted_pars[1:4], df=goodcontrast, nyears=30, TAC=1, sigma=0.000001)
#'
SPprojS <- function(inpars, df, nyears, TAC=1, sigma=0.000001){
  K <- inpars[1]
  B0 <- inpars[2]
  r <- inpars[3]
  q <- inpars[4]

  Emsy <- r/(2*q) #effort at MSY
  Bmsy <- k/2 #biomass at MSY
  MSY <- (r*k)/4 #MSY
  Fmsy <- r/2 #fishing mortality at MSY

  EBt.msy <- EBt.Emsy <- vector(length=nrow(df))

  Year <- E.msy <- E.Emsy <- C.msy <- C.Emsy <- F.msy <- F.Emsy <- vector(length=nrow(df)+nyears)
  Year <- df$year[1]:(df$year[1]+nrow(df)+nyears-1)

  # observation df
  EBt.msy[1] <- EBt.Emsy[1] <- B0
  for (i in 1:(nrow(df))){
    EBt.msy[i+1] <- EBt.msy[i]+EBt.msy[i]*r*(1-(EBt.msy[i]/K))-df[i,2]
    EBt.Emsy[i+1] <- EBt.Emsy[i]+EBt.Emsy[i]*r*(1-(EBt.Emsy[i]/K))-df[i,2]
  }

  E.msy[1:(nrow(df))] <- E.Emsy[1:(nrow(df))] <- df$effort
  C.msy[1:(nrow(df))] <- C.Emsy[1:(nrow(df))] <- df$catch

  for (i in 1:nrow(df)){
    F.msy[i] <- df$catch[i]/((EBt.msy[i]+EBt.msy[i+1])/2)
    F.Emsy[i] <- df$catch[i]/((EBt.Emsy[i]+EBt.Emsy[i+1])/2)
  }

  # projection based on MSY measure
  for (i in (nrow(df)+1):((nrow(df))+nyears)){
    b <- (r/2)*rlnorm(n=1, meanlog=(-sigma^2/2), sdlog=sigma)  # stochastic
    d <- (r/2)
    C.msy[i] <- TAC*MSY
    E.msy[i] <- C.msy[i]/(q*EBt.msy[i])
    EBt.msy[i+1] <- EBt.msy[i]+EBt.msy[i]*(b-d)*(1-(EBt.msy[i]/K))-C.msy[i]
    F.msy[i] <- C.msy[i]/((EBt.msy[i]+EBt.msy[i+1])/2)
  }

  B_B.msy <- head(EBt.msy,-1)/BMSY
  F_F.msy <- F.msy/FMSY

  temp.MSY <- data.frame(Year=Year, Catch=C.msy, Effort=E.msy, EstBt=head(EBt.msy,-1), B_B.08msy=B_B.msy, F_F.msy=F_F.msy)

  # projection based on Emsy measure
  for (i in (nrow(df)+1):((nrow(df))+nyears)){
    b <- (r/2)*rlnorm(n=1, meanlog=(-sigma^2/2), sdlog=sigma)  # stochastic
    d <- (r/2)
    E.Emsy[i] <- TAC*EMSY
    C.Emsy[i] <- q*E.Emsy[i]*EBt.Emsy[i] * rlnorm(n=1, meanlog=(-sigma^2/2), sdlog=sigma)
    EBt.Emsy[i+1] <- EBt.Emsy[i]+EBt.Emsy[i]* (b-d) *(1-(EBt.Emsy[i]/K))-C.Emsy[i]
    F.Emsy[i] <- C.Emsy[i]/((EBt.Emsy[i]+EBt.Emsy[i+1])/2)
  }

  B_B.Emsy <- head(EBt.Emsy,-1)/BMSY
  F_F.Emsy <- F.Emsy/FMSY

  temp.Emsy <- data.frame(Year=Year, Catch=C.Emsy, Effort=E.Emsy, EstBt=head(EBt.Emsy,-1), B_B.Emsy=B_B.Emsy, F_F.Emsy=F_F.Emsy)

  # projection based on 0.2Bmsy measure ?????
  # C <- q*E*B*0.2

  res <- list(temp.MSY, temp.Emsy)
  return(res)
}
