#######################################################
## Title: evaluating the management performance of underreported onewaytrip data
## Owner: Abdullah Habibi (abd.habibi@gmail.com)
## Created: May 21st, 2022
#######################################################


rm(list = ls())
# create data
K <- 1000
B0 <- K
r <- 0.2
q <- 0.00025
nYears <- 20
effort.gc <- c(seq(1,500, length.out = nYears/2), rev(seq(1,500, length.out = nYears / 2)))
effort.owt <- seq(1,500, length.out = nYears)


df.creator <- function(K, B0, r, q, nYears, effort){
  B <- CPUE <- C <- rep(NA, nYears)
  procError <- 0.05
  catchError <- 0.05
  for (i in 1:nYears) {
    if (i == 1) B[i] <- B0
    if (i>1) B[i] <- B[i-1] + rlnorm(1, meanlog = -procError^2/2, sdlog = procError) * B[i-1] * r * (1 - B[i-1]/K) - C[i-1]
    C[i] <- q * effort[i] * B[i] * rlnorm(1, meanlog = -catchError^2/2, sdlog = catchError)
    CPUE[i] <- C[i] / effort[i]
  }
  res <- data.frame(year=1:nYears, catch=C, effort=effort, CPUE=CPUE)
  return(res)
}

df.goodcontrast <- df.creator(K=1000, B0=1000, r=0.2, q=0.00025, nYears=20, effort=effort.gc)

# create list of data
nsims <- 500
owt.base <- list()
for(i in 1:nsims){
  owt.base[i] <- list(df.creator(K=1000, B0=1000, r=0.2, q=0.00025, nYears=20, effort=effort.owt))
}
owt.base

# conduct parameter estimation
pars <- c(log(K), log(B0), log(r), log(q), log(0.1))

fit.owt.T01 <- fit.owt.T02 <- fit.owt.T03 <- fit.owt.T04 <- fit.owt.T05 <- fit.owt.T06 <-
  fit.owt.T07 <- fit.owt.T08 <- fit.owt.T09 <- fit.owt.F <- list()
for(i in 1:nsims){
  fit.owt.F[i] <- list(optim(par=pars, fn=Spar_min, df=owt.base[[i]],
                             method="Nelder-Mead", OWT=FALSE, Frate = 0.7, weight = 1000))
  fit.owt.T01[i] <- list(optim(par=pars, fn=Spar_min, df=owt.base[[i]],
                                method="Nelder-Mead", OWT=TRUE, Frate = 0.7, weight = 0.1))
  fit.owt.T02[i] <- list(optim(par=pars, fn=Spar_min, df=owt.base[[i]],
                                method="Nelder-Mead", OWT=TRUE, Frate = 0.7, weight = 0.2))
  fit.owt.T03[i] <- list(optim(par=pars, fn=Spar_min, df=owt.base[[i]],
                               method="Nelder-Mead", OWT=TRUE, Frate = 0.7, weight = 0.3))
  fit.owt.T04[i] <- list(optim(par=pars, fn=Spar_min, df=owt.base[[i]],
                               method="Nelder-Mead", OWT=TRUE, Frate = 0.7, weight = 0.4))
  fit.owt.T05[i] <- list(optim(par=pars, fn=Spar_min, df=owt.base[[i]],
                               method="Nelder-Mead", OWT=TRUE, Frate = 0.7, weight = 0.5))
  fit.owt.T06[i] <- list(optim(par=pars, fn=Spar_min, df=owt.base[[i]],
                               method="Nelder-Mead", OWT=TRUE, Frate = 0.7, weight = 0.6))
  fit.owt.T07[i] <- list(optim(par=pars, fn=Spar_min, df=owt.base[[i]],
                               method="Nelder-Mead", OWT=TRUE, Frate = 0.7, weight = 0.7))
  fit.owt.T08[i] <- list(optim(par=pars, fn=Spar_min, df=owt.base[[i]],
                               method="Nelder-Mead", OWT=TRUE, Frate = 0.7, weight = 0.8))
  fit.owt.T09[i] <- list(optim(par=pars, fn=Spar_min, df=owt.base[[i]],
                               method="Nelder-Mead", OWT=TRUE, Frate = 0.7, weight = 0.9))
}

pars.owt.T01 <- pars.owt.T02 <- pars.owt.T03 <- pars.owt.T04 <- pars.owt.T05 <- pars.owt.T06 <-
  pars.owt.T07 <- pars.owt.T08 <- pars.owt.T09 <- pars.owt.F <- matrix(NA, nrow=nsims, ncol=5)
colnames(pars.owt.T01) <- colnames(pars.owt.T02) <- colnames(pars.owt.T03) <- colnames(pars.owt.T04) <-
  colnames(pars.owt.T05) <- colnames(pars.owt.T06) <- colnames(pars.owt.T07) <- colnames(pars.owt.T08) <-
  colnames(pars.owt.T09) <-colnames(pars.owt.F) <- c("K", "B0", "r", "q", "sigma")
for(i in 1:nsims){
  pars.owt.F[i,] <- exp(fit.owt.F[[i]]$par)
  pars.owt.T01[i,] <- exp(fit.owt.T01[[i]]$par)
  pars.owt.T02[i,] <- exp(fit.owt.T02[[i]]$par)
  pars.owt.T03[i,] <- exp(fit.owt.T03[[i]]$par)
  pars.owt.T04[i,] <- exp(fit.owt.T04[[i]]$par)
  pars.owt.T05[i,] <- exp(fit.owt.T05[[i]]$par)
  pars.owt.T06[i,] <- exp(fit.owt.T06[[i]]$par)
  pars.owt.T07[i,] <- exp(fit.owt.T07[[i]]$par)
  pars.owt.T08[i,] <- exp(fit.owt.T08[[i]]$par)
  pars.owt.T09[i,] <- exp(fit.owt.T09[[i]]$par)
  }

# calculate MSY
RP.owt.T01 <- RP.owt.T02 <- RP.owt.T03 <- RP.owt.T04 <- RP.owt.T05 <- RP.owt.T06 <-
  RP.owt.T07 <- RP.owt.T08 <- RP.owt.T09 <- RP.owt.F <- matrix(NA, nrow=nsims, ncol=3)
colnames(RP.owt.T01) <- colnames(RP.owt.T02) <- colnames(RP.owt.T03) <- colnames(RP.owt.T04) <- colnames(RP.owt.T05) <-
  colnames(RP.owt.T06) <- colnames(RP.owt.T07) <- colnames(RP.owt.T08) <- colnames(RP.owt.T09) <- colnames(RP.owt.F) <- c("Bmsy", "MSY", "Emsy")
for(i in 1:nsims){
  RP.owt.F[i,] <- unlist(SFcalcRP(inpars=pars.owt.F[i,], SPmodel=1))
  RP.owt.T01[i,] <- unlist(SFcalcRP(inpars=pars.owt.T01[i,], SPmodel=1))
  RP.owt.T02[i,] <- unlist(SFcalcRP(inpars=pars.owt.T02[i,], SPmodel=1))
  RP.owt.T03[i,] <- unlist(SFcalcRP(inpars=pars.owt.T03[i,], SPmodel=1))
  RP.owt.T04[i,] <- unlist(SFcalcRP(inpars=pars.owt.T04[i,], SPmodel=1))
  RP.owt.T05[i,] <- unlist(SFcalcRP(inpars=pars.owt.T05[i,], SPmodel=1))
  RP.owt.T06[i,] <- unlist(SFcalcRP(inpars=pars.owt.T06[i,], SPmodel=1))
  RP.owt.T07[i,] <- unlist(SFcalcRP(inpars=pars.owt.T07[i,], SPmodel=1))
  RP.owt.T08[i,] <- unlist(SFcalcRP(inpars=pars.owt.T08[i,], SPmodel=1))
  RP.owt.T09[i,] <- unlist(SFcalcRP(inpars=pars.owt.T09[i,], SPmodel=1))
}

plot(RP.owt.F[,1], col="red")
points(RP.owt.T01[,1], col="blue1")

# see different effect of weighting to the OWT
# higher weight, we are effectively constraining the estimation procedure to fit the auxiliary information exactly.
# lower weight, assuming that an auxiliary observation was measured with low accuracy (high variance)
Sparbase <- data.frame(OWT=rep(c("F", "T01", "T02", "T03", "T04", "T05", "T06", "T07", "T08", "T09"), each=nsims),
                       K=c(pars.owt.F[,1], pars.owt.T01[,1], pars.owt.T02[,1], pars.owt.T03[,1], pars.owt.T04[,1],
                           pars.owt.T05[,1], pars.owt.T06[,1], pars.owt.T07[,1], pars.owt.T08[,1], pars.owt.T09[,1]),
                       B0=c(pars.owt.F[,2], pars.owt.T01[,2], pars.owt.T02[,2], pars.owt.T03[,2], pars.owt.T04[,2],
                            pars.owt.T05[,2], pars.owt.T06[,2], pars.owt.T07[,2], pars.owt.T08[,2], pars.owt.T09[,2]),
                       r=c(pars.owt.F[,3], pars.owt.T01[,3], pars.owt.T02[,3], pars.owt.T03[,3], pars.owt.T04[,3],
                           pars.owt.T05[,3], pars.owt.T06[,3], pars.owt.T07[,3], pars.owt.T08[,3], pars.owt.T09[,3]),
                       q=c(pars.owt.F[,4], pars.owt.T01[,4], pars.owt.T02[,4], pars.owt.T03[,4], pars.owt.T04[,4],
                           pars.owt.T05[,4], pars.owt.T06[,4], pars.owt.T07[,4], pars.owt.T08[,4], pars.owt.T09[,4]),
                       sigma=c(pars.owt.F[,5], pars.owt.T01[,5], pars.owt.T02[,5], pars.owt.T03[,5], pars.owt.T04[,5],
                               pars.owt.T05[,5], pars.owt.T06[,5], pars.owt.T07[,5], pars.owt.T08[,5], pars.owt.T09[,5]))


par(mfrow=c(1, 5))
boxplot(Sparbase[,2]~Sparbase[,1], xlab="Weighting", ylab="K", ylim=c(0,2000))
boxplot(Sparbase[,3]~Sparbase[,1], xlab="Weighting", ylab="B0", ylim=c(0,2000))
boxplot(Sparbase[,4]~Sparbase[,1], xlab="Weighting", ylab="r", ylim=c(0,1))
boxplot(Sparbase[,5]~Sparbase[,1], xlab="Weighting", ylab="q")
boxplot(Sparbase[,6]~Sparbase[,1], xlab="Weighting", ylab="sigma")

RPbase <- data.frame(OWT=rep(c("F", "T01", "T02", "T03", "T04", "T05", "T06", "T07", "T08", "T09"), each=nsims),
                     Bmsy=c(RP.owt.F[,1], RP.owt.T01[,1], RP.owt.T02[,1], RP.owt.T03[,1], RP.owt.T04[,1],
                            RP.owt.T05[,1], RP.owt.T06[,1], RP.owt.T07[,1], RP.owt.T08[,1], RP.owt.T09[,1]),
                     MSY=c(RP.owt.F[,2], RP.owt.T01[,2], RP.owt.T02[,2], RP.owt.T03[,2], RP.owt.T04[,2],
                           RP.owt.T05[,2], RP.owt.T06[,2], RP.owt.T07[,2], RP.owt.T08[,2], RP.owt.T09[,2]),
                     Emsy=c(RP.owt.F[,3], RP.owt.T01[,3], RP.owt.T02[,3], RP.owt.T03[,3], RP.owt.T04[,3],
                            RP.owt.T05[,3], RP.owt.T06[,3], RP.owt.T07[,3], RP.owt.T08[,3], RP.owt.T09[,3]))

par(mfrow=c(1, 3))
boxplot(RPbase[,2]~RPbase[,1], xlab="Weighting", ylab="Bmsy", ylim=c(0, 1000))
boxplot(RPbase[,3]~RPbase[,1], xlab="Weighting", ylab="MSY", ylim=c(0, 100))
boxplot(RPbase[,4]~RPbase[,1], xlab="Weighting", ylab="Emsy", ylim=c(0, 1000))

# create projection
msyproj.gc.base <- emsyproj.gc.base <- list()
for(i in 1:nsims){
  msyproj.gc.base[i] <- list(Sproj(inpars=pars.gc.base[i,], df=gc.base[[i]], nyears=30, TAC=1, sigma=0.000001)[[1]])
  emsyproj.gc.base[i] <- list(Sproj(inpars=pars.gc.base[i,], df=gc.base[[i]], nyears=30, TAC=1, sigma=0.000001)[[2]])
}
msyproj.gc.base
emsyproj.gc.base

# plot the data and projection
sapply(gc.base, "[[", 3)[,1] # grabbing only first simulation on effort data from all simulations

plotted.col <- 5 # colnames = c(Year, Catch, Effort, EstBt, B_B.msy, E_E.msy)
yval_msyproj.gc.base <- sapply(msyproj.gc.base, "[[", plotted.col)
yval_emsyproj.gc.base <- sapply(emsyproj.gc.base, "[[", plotted.col)

par(mfrow=c(2,1))
plot(x=1:50, y=rowMeans(yval_msyproj.gc.base), col="red", lwd=2, xlim=c(1,50), ylim=c(0, 2.7), type="l",
     ylab="B/Bmsy", xlab="Year") # mean value from all simulations
lines(x=1:50, y=t(apply(yval_msyproj.gc.base,1,quantile,c(0.025,0.975)))[,1], lty=2, col="red", lwd=2) # lower confident interval
lines(x=1:50, y=t(apply(yval_msyproj.gc.base,1,quantile,c(0.025,0.975)))[,2], lty=2, col="red", lwd=2) # upper confident interval
lines(x=1:50, y=rowMeans(yval_emsyproj.gc.base), col="blue", lwd=2) # mean value from all simulations
lines(x=1:50, y=t(apply(yval_emsyproj.gc.base,1,quantile,c(0.025,0.975)))[,1], lty=2, col="blue", lwd=2) # lower confident interval
lines(x=1:50, y=t(apply(yval_emsyproj.gc.base,1,quantile,c(0.025,0.975)))[,2], lty=2, col="blue", lwd=2) # upper confident interval

## need to be added to the minimization with normal distribution version
#' @param dt type of data distribution in the abundance index. The default is "log-normal",
#' and can be replaced with "normal" when the abundance index show a closer pattern to normal distribution.
#'
#' dt=c("log-normal", "normal")
#'
if (OWT==FALSE){
ifelse(dt=="log-normal", nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = sigma, log = TRUE)),
       ifelse(dt=="normal", nll <- -sum(dnorm(x= na.omit(CPUE), mean = na.omit(EstCPUE), sd = sigma, log = TRUE)),
              nll <- "wrong dt code!"))
}
else{
  ifelse(dt=="log-normal", nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = sigma, log = TRUE)) +
           weight * (tail(annualFrates,1) - Frate)^2,
         ifelse(dt=="normal", nll <- -sum(dnorm(x= na.omit(CPUE), mean = na.omit(EstCPUE), sd = sigma, log = TRUE)) +
                  weight * (tail(annualFrates,1) - Frate)^2,
                nll <- "wrong dt code!"))
}

