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

df.creator(K=1000, B0=1000, r=0.2, q=0.00025, nYears=20, effort=effort.owt)

# create list of data
nsims <- 500
owt.base <- list()
for(i in 1:nsims){
  owt.base[i] <- list(df.creator(K=1000, B0=1000, r=0.2, q=0.00025, nYears=20, effort=effort.owt))
}
owt.base

# conduct parameter estimation
pars <- c(log(K), log(B0), log(r), log(q), log(0.1))

fit.owt.T250 <- fit.owt.T200 <- fit.owt.T300 <- fit.owt.T350 <- fit.owt.F <- list()
for(i in 1:nsims){
  fit.owt.F[i] <- list(optim(par=pars, fn=Spar_min, df=owt.base[[i]],
                             method="Nelder-Mead", OWT=FALSE, Frate = 0.7, weight = 1000))
  fit.owt.T250[i] <- list(optim(par=pars, fn=Spar_min, df=owt.base[[i]],
                                method="Nelder-Mead", OWT=TRUE, Frate = 0.7, weight = 250))
  fit.owt.T200[i] <- list(optim(par=pars, fn=Spar_min, df=owt.base[[i]],
                                method="Nelder-Mead", OWT=TRUE, Frate = 0.7, weight = 200))
  fit.owt.T300[i] <- list(optim(par=pars, fn=Spar_min, df=owt.base[[i]],
                                method="Nelder-Mead", OWT=TRUE, Frate = 0.7, weight = 300))
  fit.owt.T350[i] <- list(optim(par=pars, fn=Spar_min, df=owt.base[[i]],
                                method="Nelder-Mead", OWT=TRUE, Frate = 0.7, weight = 350))
}

pars.owt.T250 <- pars.owt.T200 <- pars.owt.T300 <- pars.owt.T350 <- pars.owt.F <- matrix(NA, nrow=nsims, ncol=5)
colnames(pars.owt.T250) <- colnames(pars.owt.T200) <- colnames(pars.owt.T300) <-
  colnames(pars.owt.T350) <-colnames(pars.owt.F) <- c("K", "B0", "r", "q", "sigma")
for(i in 1:nsims){
  pars.owt.F[i,] <- exp(fit.owt.F[[i]]$par)
  pars.owt.T250[i,] <- exp(fit.owt.T250[[i]]$par)
  pars.owt.T200[i,] <- exp(fit.owt.T200[[i]]$par)
  pars.owt.T300[i,] <- exp(fit.owt.T300[[i]]$par)
  pars.owt.T350[i,] <- exp(fit.owt.T350[[i]]$par)
}

# calculate MSY
RP.owt.T250 <- RP.owt.T200 <- RP.owt.T300 <- RP.owt.T350 <- RP.owt.F <- matrix(NA, nrow=nsims, ncol=3)
colnames(RP.owt.T250) <- colnames(RP.owt.T200) <- colnames(RP.owt.T300) <-
  colnames(RP.owt.T350) <- colnames(RP.owt.F) <- c("Bmsy", "MSY", "Emsy")
for(i in 1:nsims){
  RP.owt.F[i,] <- unlist(SFcalcRP(inpars=pars.owt.F[i,], SPmodel=1))
  RP.owt.T250[i,] <- unlist(SFcalcRP(inpars=pars.owt.T250[i,], SPmodel=1))
  RP.owt.T200[i,] <- unlist(SFcalcRP(inpars=pars.owt.T200[i,], SPmodel=1))
  RP.owt.T300[i,] <- unlist(SFcalcRP(inpars=pars.owt.T300[i,], SPmodel=1))
  RP.owt.T350[i,] <- unlist(SFcalcRP(inpars=pars.owt.T350[i,], SPmodel=1))
}

plot(RP.owt.F[,1], col="red")
points(RP.owt.T250[,1], col="blue1")
points(RP.owt.T200[,1], col="blue2")
points(RP.owt.T300[,1], col="blue3")
points(RP.owt.T350[,1], col="blue4")

# see different effect of weighting to the OWT
# how do we know ideal weighing value, because T250 gives median r 0.2 ? but K and q doesn't give a close value like r!
Sparbase <- data.frame(OWT=rep(c("F", "T250", "T200", "T300", "T350"), each=nsims),
                       K=c(pars.owt.F[,1], pars.owt.T250[,1], pars.owt.T200[,1], pars.owt.T300[,1], pars.owt.T350[,1]),
                       B0=c(pars.owt.F[,2], pars.owt.T250[,2], pars.owt.T200[,2], pars.owt.T300[,2], pars.owt.T350[,2]),
                       r=c(pars.owt.F[,3], pars.owt.T250[,3], pars.owt.T200[,3], pars.owt.T300[,3], pars.owt.T350[,3]),
                       q=c(pars.owt.F[,4], pars.owt.T250[,4], pars.owt.T200[,4], pars.owt.T300[,4], pars.owt.T350[,4]),
                       sigma=c(pars.owt.F[,5], pars.owt.T250[,5], pars.owt.T200[,5], pars.owt.T300[,5], pars.owt.T350[,5]))

par(mfrow=c(1, 5))
boxplot(Sparbase[,2]~Sparbase[,1], xlab="Weighting", ylab="K", ylim=c(0,3000))
boxplot(Sparbase[,3]~Sparbase[,1], xlab="Weighting", ylab="B0", ylim=c(0,3000))
boxplot(Sparbase[,4]~Sparbase[,1], xlab="Weighting", ylab="r", ylim=c(0,3))
boxplot(Sparbase[,5]~Sparbase[,1], xlab="Weighting", ylab="q")
boxplot(Sparbase[,6]~Sparbase[,1], xlab="Weighting", ylab="sigma")

RPbase <- data.frame(OWT=rep(c("F", "T250", "T200", "T300", "T350"), each=nsims),
                     Bmsy=c(RP.owt.F[,1], RP.owt.T250[,1], RP.owt.T200[,1], RP.owt.T300[,1], RP.owt.T350[,1]),
                     MSY=c(RP.owt.F[,2], RP.owt.T250[,2], RP.owt.T200[,2], RP.owt.T300[,2], RP.owt.T350[,2]),
                     Emsy=c(RP.owt.F[,3], RP.owt.T250[,3], RP.owt.T200[,3], RP.owt.T300[,3], RP.owt.T350[,3]))

par(mfrow=c(1, 3))
boxplot(RPbase[,2]~RPbase[,1], xlab="Weighting", ylab="Bmsy", ylim=c(0, 2000))
boxplot(RPbase[,3]~RPbase[,1], xlab="Weighting", ylab="MSY", ylim=c(0, 200))
boxplot(RPbase[,4]~RPbase[,1], xlab="Weighting", ylab="Emsy")

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


