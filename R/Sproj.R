#' @title Schaefer's projection function
#'
#' @description
#' Function to create projection in Schaefer's model based on different reference points (Bmsy, MSY, Emsy) and
#' different Total Allowable Catch setting (TAC). This projection can be used to assist in the development of
#' limit and target reference point.
#'
#' @param inpars fitted surplus production parameters which consist of K (carrying capacity),
#' B0 (biomass when fishing is started), r (intrinsic growth rate), q (catchability coefficient)
#' @param df dataframe containing three columns; year, catch and effort
#' @param nyears number of years the projection for the fishery
#' @param TAC number to drive the management level
#' @param sigma standard deviation estimated in the surplus production data fitting process
#'
#' @return
#' sigma in default is set at a very low value (0.000001) to allow for a deterministic projection. In
#' interest to create stochastic projection, the sigma value produced during maximum likelihood estimation
#' or use another value can be utilized
#'
#' The TAC is set at 1 in default. It would depend on the fishery to set whether the TAC can be set lower
#' to be more conservative or set higher to increase more catch.
#' For instance, when the TAC is set conservative at 0.8 of reference point, it will return (1-0.8)Bmsy,
#' 0.8MSY, and 0.8Emsy. In contrary, the TAC will increase catch when it is set at 1.2
#' (or other value higher than 1) and return (1*1.2)Bmsy, 1.2MSY, and 1.2Emsy.
#'
#' @export
#'
#' @examples
#' res <- Sproj(inpars=fitted_pars[1:4], df=df.goodcontrast, nyears=30, TAC=1, sigma=0.000001)
#'
#' # create simple comparison B/Bmsy and F/Fmsy plotting
#' plot(res[[1]][,5], type="l", col="red", ylab="B/Bmsy", xlab="year")
#' lines(res[[2]][,5], type="l", col="black")
#' lines(res[[3]][,5], type="l", col="blue")
#' legend("bottomleft", col=c("red", "black", "blue"), legend=c("MSY", "Emsy", "Bmsy"), lty=1)
#'
#' plot(res[[1]][,6], type="l", col="red", ylab="F/Fmsy", xlab="year")
#' lines(res[[2]][,6], type="l", col="black")
#' lines(res[[3]][,6], type="l", col="blue")
#' legend("bottomleft", col=c("red", "black", "blue"), legend=c("MSY", "Emsy", "Bmsy"), lty=1)
#'
#'

# The TAC*Bmsy need to have further scrutiny

Sproj <- function(inpars, df, nyears, TAC=1, sigma=0.000001){
  K <- inpars[1]
  B0 <- inpars[2]
  r <- inpars[3]
  q <- inpars[4]

  Emsy <- r/(2*q) #effort at MSY
  Bmsy <- K/2 #biomass at MSY
  MSY <- (r*K)/4 #MSY
  Fmsy <- r/2 #fishing mortality at MSY

  EBt.msy <- EBt.Emsy  <- EBt.Bmsy <- vector(length=nrow(df))

  Year <- E.msy <- E.Emsy <- E.Bmsy <- C.msy <- C.Emsy <- C.Bmsy <- F.msy <- F.Emsy <- F.Bmsy <- vector(length=nrow(df)+nyears)
  Year <- df$year[1]:(df$year[1]+nrow(df)+nyears-1)

  # observation df
  EBt.msy[1] <- EBt.Emsy[1] <- EBt.Bmsy[1] <- B0
  for (i in 1:(nrow(df))){
    EBt.msy[i+1] <- EBt.msy[i]+EBt.msy[i]*r*(1-(EBt.msy[i]/K))-df[i,2]
    EBt.Emsy[i+1] <- EBt.Emsy[i]+EBt.Emsy[i]*r*(1-(EBt.Emsy[i]/K))-df[i,2]
    EBt.Bmsy[i+1] <- EBt.Bmsy[i]+EBt.Bmsy[i]*r*(1-(EBt.Bmsy[i]/K))-df[i,2]
  }

  E.msy[1:(nrow(df))] <- E.Emsy[1:(nrow(df))] <- E.Bmsy[1:(nrow(df))] <- df$effort
  C.msy[1:(nrow(df))] <- C.Emsy[1:(nrow(df))] <- C.Bmsy[1:(nrow(df))] <- df$catch

  for (i in 1:nrow(df)){
    F.msy[i] <- df$catch[i]/((EBt.msy[i]+EBt.msy[i+1])/2)
    F.Emsy[i] <- df$catch[i]/((EBt.Emsy[i]+EBt.Emsy[i+1])/2)
  }

  # projection based on MSY measure
  for (i in (nrow(df)+1):((nrow(df))+nyears)){
    C.msy[i] <- TAC*MSY
    E.msy[i] <- C.msy[i]/(q*EBt.msy[i])
    EBt.msy[i+1] <- (EBt.msy[i]+EBt.msy[i]*r*(1-(EBt.msy[i]/K))-C.msy[i])*rlnorm(n=1, meanlog=(-sigma^2/2), sdlog=sigma)  # stochastic
    F.msy[i] <- C.msy[i]/((EBt.msy[i]+EBt.msy[i+1])/2)
  }

  B_B.msy <- head(EBt.msy,-1)/Bmsy
  F_F.msy <- F.msy/Fmsy

  temp.MSY <- data.frame(Year=Year, Catch=C.msy, Effort=E.msy, EstBt=head(EBt.msy,-1), B_B.08msy=B_B.msy, F_F.msy=F_F.msy)

  # projection based on Emsy measure
  for (i in (nrow(df)+1):((nrow(df))+nyears)){
    E.Emsy[i] <- TAC*Emsy
    C.Emsy[i] <- q*E.Emsy[i]*EBt.Emsy[i] * rlnorm(n=1, meanlog=(-sigma^2/2), sdlog=sigma)
    EBt.Emsy[i+1] <- (EBt.Emsy[i]+EBt.Emsy[i]*r*(1-(EBt.Emsy[i]/K))-C.Emsy[i])*rlnorm(n=1, meanlog=(-sigma^2/2), sdlog=sigma)  # stochastic
    F.Emsy[i] <- C.Emsy[i]/((EBt.Emsy[i]+EBt.Emsy[i+1])/2)
  }

  B_B.Emsy <- head(EBt.Emsy,-1)/Bmsy
  F_F.Emsy <- F.Emsy/Fmsy

  temp.Emsy <- data.frame(Year=Year, Catch=C.Emsy, Effort=E.Emsy, EstBt=head(EBt.Emsy,-1), B_B.Emsy=B_B.Emsy, F_F.Emsy=F_F.Emsy)

  # projection based on Bmsy measure
  for (i in (nrow(df)+1):((nrow(df))+nyears)){
    C.Bmsy[i] <- ifelse(TAC < 1, (1-TAC), TAC)*Bmsy
    EBt.Bmsy[i+1] <- (EBt.msy[i]+EBt.Bmsy[i]*r*(1-(EBt.Bmsy[i]/K))-C.Bmsy[i])*rlnorm(n=1, meanlog=(-sigma^2/2), sdlog=sigma)  # stochastic
    E.Bmsy[i] <- C.Bmsy[i]/(q*EBt.Bmsy[i])
    F.Bmsy[i] <- C.Bmsy[i]/((EBt.Bmsy[i]+EBt.Bmsy[i+1])/2)
  }

  B_B.Bmsy <- head(EBt.Bmsy,-1)/Bmsy
  F_F.Bmsy <- F.Bmsy/Fmsy

  temp.Bmsy <- data.frame(Year=Year, Catch=C.Bmsy, Effort=E.Bmsy, EstBt=head(EBt.Bmsy,-1), B_B.Bmsy=B_B.Bmsy, F_F.Bmsy=F_F.Bmsy)

  res <- list(temp.MSY, temp.Emsy, temp.Bmsy)
  return(res)
}
