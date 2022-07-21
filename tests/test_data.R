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
repsims <- function(nYears, nsims){

  for (year in 1:length(nYears)){
    #nYears=nYears[year]

    B <- CPUE <- C <- rep(NA, nYears[year])
    effort.gc <- c(seq(1,500, length.out = nYears[year]/2), rev(seq(1,500, length.out = nYears[year]/2))) # good contrast
    effort.owt <- seq(1,500, length.out = nYears[year]) # one way trip

    dat.gc <- replicate(n=nsims,
                        expr={for (i in 1:nYears[year]) {
                          if (i == 1) B[i] <- B0
                          if (i>1) B[i] <- B[i-1] + rlnorm(1, meanlog = -procError^2/2, sdlog = procError) * B[i-1] * r * (1 - B[i-1]/K) - C[i-1]
                          C[i] <- q * effort.gc[i] * B[i] * rlnorm(1, meanlog = -catchError^2/2, sdlog = catchError)
                          CPUE[i] <- C[i] / effort.gc[i]
                        }
                          data.frame(year=1:nYears[year], catch = C, effort = effort.gc, CPUE = CPUE)},
                        simplify=F)

    dat.owt <- replicate(n=nsims,
                         expr={for (i in 1:nYears[year]) {
                           if (i == 1) B[i] <- B0
                           if (i>1) B[i] <- B[i-1] + rlnorm(1, meanlog = -procError^2/2, sdlog = procError) * B[i-1] * r * (1 - B[i-1]/K) - C[i-1]
                           C[i] <- q * effort.owt[i] * B[i] * rlnorm(1, meanlog = -catchError^2/2, sdlog = catchError)
                           CPUE[i] <- C[i] / effort.owt[i]
                         }
                           data.frame(year=1:nYears[year], catch = C, effort = effort.owt, CPUE = CPUE)},
                         simplify=F)

    #par(mfrow=c(1,2))
    #plot(dat.gc[[1]]$CPUE, type='l', ylab="CPUE", xlab="Year", col="grey", main=paste('good-contrast data |', nsims, "sims"))
    #for(sims in 2:nsims) lines(dat.gc[[sims]]$CPUE, type='l', col='grey')
    #plot(dat.owt[[1]]$CPUE, type='l', ylab="CPUE", xlab="Year", col="grey", main=paste('one-way-trip data |', nsims, "sims"))
    #for(sims in 2:nsims) lines(dat.owt[[sims]]$CPUE, type='l', col='grey')


    ### Estimate parameters using optim
    res.owtpl <- res.owtnpl <- res.gc <- data.frame(matrix(NA, nrow=nsims, ncol=5))

    for (i in 1:nsims){ # for good contrast data
      temp <- optim(par=inpars, fn=Par.min, df=dat.gc[[i]],
                    method="Nelder-Mead")#, OWT="Depletion", currentF = 0.7, weight = 1000)
      for (j in 1:5){ res.gc[i,j] <- exp(temp$par[j])}
    }

    for (i in 1:nsims){ # for one way trip data, analysed without penalized likelihood
      temp <- optim(par=inpars, fn=Par.min, df=dat.owt[[i]],
                    method="Nelder-Mead")#, OWT="Depletion", currentF = 0.7, weight = 1000)
      for (j in 1:5){ res.owtnpl[i,j] <- exp(temp$par[j])}
    }

    for (i in 1:nsims){ # for one way trip data, analysed with penalized likelihood
      temp <- optim(par=inpars, fn=Par.min, df=dat.owt[[i]],
                    method="Nelder-Mead", OWT="Depletion", currentF = 0.7, weight = 1000)
      for (j in 1:5){ res.owtpl[i,j] <- exp(temp$par[j])}
    }

    colnames(res.owtpl) <- colnames(res.owtnpl) <- colnames(res.gc) <- c("K", "B0", "r", "q", "obs.err")

    MSY.calc <- function(df){
      return(cbind(df, MSY=df$r*df$K/4, fMSY=df$r/(2*df$q))) # adding MSY and fMSY columns
    }

    res <- rbind(res.gc, res.owtnpl, res.owtpl)
    #res <- MSY.calc(res)
    res <- cbind(res, tipe=rep(c("good-contrast", "one-way-trip", "one-way-trip_fix"), each=nsims))
    return(res)
  }
  all_res <- rbind(res[year])
  return(all_res)
}



inputyears=c(10, 20, 30)
inputsims=10
repsims(nYears=inputyears, nsims=inputsims)























df10
df20
df30
yvec <- c(10,20,30)
replicate(n=yvec, repsims, simplify = F)

# combine data in one dataframe
dat_all <- rbind(df10, df20, df30)
dat_all <- cbind(dat_all, ndata=rep(c("10", "20", "30"), each=nsims*3))
head(dat_all)

# conduct quick analysis
library(dplyr)

dat_all %>% group_by(ndata, tipe) %>%
  summarize_each(funs(sd=sd(., na.rm=T), se=sd(., na.rm=T)/sqrt(sum(!is.na(.)))))

# plot the estimated parameters
library(ggplot2)
library(tidyr)

dat_all[,-c(1:5)] %>%
  pivot_longer(-c(tipe,ndata),'parameter','value') %>%
  ggplot(aes(ndata, value)) + geom_violin() +
  facet_wrap(tipe~parameter, scales = 'free_y', nrow=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.title.y = element_blank()) +
  labs(x = "Number of annual time series data")
