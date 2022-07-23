# Title: Data number and type effects
# Written by: Abdullah Habibi
# Coded since: July 4th, 2022

rm(list=ls())


devtools::install_github("habeebollah/montiR")
library(montiR)

# set surplus production parameters
r <- 0.2
K <- 1000
q <- 0.00025
B0 <- K
procError <- 0.01
catchError <- 0.01
sigma <- 0.1
inpars <- c(log(K), log(B0), log(r), log(q), log(sigma))

# functions to create fake data and calculate the surplus production parameters

repsims <- function(nsims=1, nYears){
    B <- CPUE <- C <- rep(NA, nYears)
    effort.gc <- c(seq(1,500, length.out = nYears/2), rev(seq(1,500, length.out = nYears/2))) # good contrast
    effort.owt <- seq(1,500, length.out = nYears) # one way trip
    effort.f <- seq(1,1, length.out = nYears) # flat data

    dat.gc <- replicate(n=nsims,
                        expr={for (i in 1:nYears) {
                          if (i == 1) B[i] <- B0
                          if (i>1) B[i] <- B[i-1] + rlnorm(1, meanlog = -procError^2/2, sdlog = procError) * B[i-1] * r * (1 - B[i-1]/K) - C[i-1]
                          C[i] <- q * effort.gc[i] * B[i] * rlnorm(1, meanlog = -catchError^2/2, sdlog = catchError)
                          CPUE[i] <- C[i] / effort.gc[i]
                        }
                          data.frame(year=1:nYears, catch = C, effort = effort.gc, CPUE = CPUE)},
                        simplify=F)

    dat.f <- replicate(n=nsims,
                       expr={for (i in 1:nYears) {
                         if (i == 1) B[i] <- B0
                         if (i>1) B[i] <- B[i-1] + rlnorm(1, meanlog = -procError^2/2, sdlog = procError) * B[i-1] * r * (1 - B[i-1]/K) - C[i-1]
                         C[i] <- q * effort.f[i] * B[i] * rlnorm(1, meanlog = -catchError^2/2, sdlog = catchError)
                         CPUE[i] <- C[i] / effort.f[i]
                       }
                         data.frame(year=1:nYears, catch = C, effort = effort.f, CPUE = CPUE)},
                       simplify=F)

    dat.owt <- replicate(n=nsims,
                         expr={for (i in 1:nYears) {
                           if (i == 1) B[i] <- B0
                           if (i>1) B[i] <- B[i-1] + rlnorm(1, meanlog = -procError^2/2, sdlog = procError) * B[i-1] * r * (1 - B[i-1]/K) - C[i-1]
                           C[i] <- q * effort.owt[i] * B[i] * rlnorm(1, meanlog = -catchError^2/2, sdlog = catchError)
                           CPUE[i] <- C[i] / effort.owt[i]
                         }
                           data.frame(year=1:nYears, catch = C, effort = effort.owt, CPUE = CPUE)},
                         simplify=F)

    par(mfrow=c(3,1))
    plot(dat.gc[[1]]$CPUE, type='l', ylab="CPUE", xlab="Year", col="grey", main=paste('good-contrast data |', nsims, "sims"))
    for(sims in 2:nsims) lines(dat.gc[[sims]]$CPUE, type='l', col='grey')
    plot(dat.f[[1]]$CPUE, type='l', ylab="CPUE", xlab="Year", col="grey", main=paste('flat data |', nsims, "sims"))
    for(sims in 2:nsims) lines(dat.f[[sims]]$CPUE, type='l', col='grey')
    plot(dat.owt[[1]]$CPUE, type='l', ylab="CPUE", xlab="Year", col="grey", main=paste('one-way-trip data |', nsims, "sims"))
    for(sims in 2:nsims) lines(dat.owt[[sims]]$CPUE, type='l', col='grey')

    ### Estimate parameters using optim
    res.owtpl <- res.owtnpl <- res.fpl <- res.fnpl <- res.gc <- data.frame(matrix(NA, nrow=nsims, ncol=5))
    res <- NULL

    for (i in 1:nsims){
      temp1 <- optim(par=inpars, fn=Par.min, df=dat.gc[[i]], method="Nelder-Mead")
      temp2 <- optim(par=inpars, fn=Par.min, df=dat.f[[i]], method="Nelder-Mead")
      temp3 <- optim(par=inpars, fn=Par.min, df=dat.f[[i]],
                     method="Nelder-Mead", OWT="Depletion", currentF = 0.3, weight = 1000)
      temp4 <- optim(par=inpars, fn=Par.min, df=dat.owt[[i]], method="Nelder-Mead")
      temp5 <- optim(par=inpars, fn=Par.min, df=dat.owt[[i]],
                     method="Nelder-Mead", OWT="Depletion", currentF = 0.7, weight = 1000)

      for (j in 1:5){
        res.gc[i,j] <- exp(temp1$par[j]) # for good contrast data
        res.fnpl[i,j] <- exp(temp2$par[j]) # for flat data, analysed without penalized likelihood
        res.fpl[i,j] <- exp(temp3$par[j]) # for flat data, analysed with penalized likelihood
        res.owtnpl[i,j] <- exp(temp4$par[j]) # for one way trip data, analysed without penalized likelihood
        res.owtpl[i,j] <- exp(temp5$par[j]) # for one way trip data, analysed with penalized likelihood

      }

    }
    colnames(res.owtpl) <- colnames(res.owtnpl) <- colnames(res.fpl) <-colnames(res.fnpl) <- colnames(res.gc) <- c("K", "B0", "r", "q", "obs.err")

    MSY.calc <- function(df){
      return(cbind(df, MSY=df$r*df$K/4, fMSY=df$r/(2*df$q))) # adding MSY and fMSY columns
    }
    res <- MSY.calc(rbind(res.gc, res.fnpl, res.fpl, res.owtnpl, res.owtpl))
    res <- cbind(res, tipe=rep(c("good-contrast", "flat", "flat_fix", "one-way-trip", "one-way-trip_fix"), each=nsims))
    return(res)
}

nsims=10
nYears=10
df10 <- repsims(nsims=nsims, nYears=nYears)

df10
df20
df30

# combine data in one dataframe
dat_all <- cbind(rbind(df10, df20, df30),
                 ndata=rep(c('10', '20', '30'), each=nsims*5))

head(dat_all)
unique(dat_all$ndata)

# conduct quick analysis
library(dplyr)

sims2000_02 <- dat_all %>% group_by(ndata, tipe) %>%
  summarize_each(funs(mean=mean(., na.rm=T),min=min(., na.rm=T), max=max(., na.rm=T),sd=sd(., na.rm=T), se=sd(., na.rm=T)/sqrt(sum(!is.na(.)))))
write.csv(sims2000_02, 'sims2000_02.csv')

# plot the estimated parameters
library(ggplot2)
library(tidyr)

dat_all[,-c(1:5)] %>%
  pivot_longer(-c(tipe,ndata),'parameter','value') %>%
  ggplot(aes(ndata, value)) + geom_violin() +
  facet_wrap(tipe~parameter, scales = 'free_y', nrow=5) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.title.y = element_blank()) +
  labs(x = "Number of annual time series data")


