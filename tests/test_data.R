# Title: Data number and type effects
# Written by: Abdullah Habibi
# Coded since: July 4th, 2022

rm(list=ls())


#devtools::install_github("habeebollah/montiR")
library(montiR)

# set surplus production parameters
r <- 0.2
K <- 1000
q <- 0.00025
B0 <- K
procError <- 0.05
catchError <- 0.05
sigma <- 0.1
inpars <- c(log(K), log(B0), log(r), log(q), log(sigma))

nYears <- 10
effort <- c(seq(1,500, length.out = nYears/2), rev(seq(1,500, length.out = nYears / 2))) # good contrast
#effort <- seq(1,500, length.out = nYears) # one way trip

# create fake data
B <- CPUE <- C <- rep(NA, nYears)
nsims <- 500

dat <- replicate(n=nsims,
                     expr={for (i in 1:nYears) {
                       if (i == 1) B[i] <- B0
                       if (i>1) B[i] <- B[i-1] + rlnorm(1, meanlog = -procError^2/2, sdlog = procError) * B[i-1] * r * (1 - B[i-1]/K) - C[i-1]
                       C[i] <- q * effort[i] * B[i] * rlnorm(1, meanlog = -catchError^2/2, sdlog = catchError)
                       CPUE[i] <- C[i] / effort[i]
                     }
                       data.frame(year=1:nYears, catch = C, effort = effort)},
                     simplify=F)

### Estimate parameters using optim
res <- data.frame(matrix(NA, nrow=nsims, ncol=5))
for (i in 1:nsims){
  temp <- optim(par=inpars,
                       fn=Par.min,
                       df=dat[[i]],
                       method="Nelder-Mead")#, OWT="Depletion", currentF = 0.7, weight = 1000)
  for (j in 1:5){
    res[i,j] <- exp(temp$par[j])
  }
}
colnames(res) <- c("K", "B0", "r", "q", "obs.err")

MSY.calc <- function(df){
  return(cbind(df, MSY=df$r*df$K/4, fMSY=df$r/(2*df$q))) # adding MSY and fMSY columns
  }

col.analysis <- function(df){
  return(data.frame(SE=sd(df)/sqrt(length(df)),
                    min=min(df), max=max(df), mean=mean(df))) # calculate the standard error
}

res <- MSY.calc(res)
sapply(res, col.analysis)[c(1,4),]

gc_res10 <- res # replicate from line 21 to create data with 10, 20, 30 annual data

# combine data in one dataframe
gc_all <- rbind(gc_res10,gc_res20,gc_res30)
gc_all <- cbind(gc_all,tipe=rep(c("10", "20", "30"), each=nsims))

owt_all <- rbind(owt_res10,owt_res20,owt_res30)
owt_all <- cbind(owt_all,tipe=rep(c("10", "20", "30"), each=nsims))

owtfix_all <- rbind(owtfix_res10,owtfix_res20,owtfix_res30)
owtfix_all <- cbind(owtfix_all,tipe=rep(c("10", "20", "30"), each=nsims))

dat_all <- rbind(gc_all, owt_all, owtfix_all)
dat_all <- cbind(dat_all, data=rep(c("good-contrast", "one-way-trip", "one-way-trip_fix"), each=nsims*3))

# plot the estimated parameters
library(tidyr)
library(ggplot2)

gc_all[,-c(5:7)] %>%
  pivot_longer(-tipe,'parameter','value') %>%
  ggplot(aes(tipe, value, color = parameter)) +
  geom_violin() +
  facet_wrap(~parameter, scales = 'free_y', nrow=1)

owt_all[,-5] %>%
  pivot_longer(-tipe,'parameter','value') %>%
  ggplot(aes(tipe, value, color=parameter)) +
  geom_violin() +
  facet_wrap(~parameter, scales = 'free_y', nrow=1)

owtfix_all[,-5] %>%
  pivot_longer(-tipe,'parameter','value') %>%
  ggplot(aes(tipe, value)) +
  geom_violin() + facet_wrap(~parameter, scales = 'free_y', nrow=1)

dat_all[,-c(5:7)] %>%
  pivot_longer(-c(tipe,data),'parameter','value') %>%
  ggplot(aes(tipe, value)) + geom_violin() +
  facet_wrap(data~parameter, scales = 'free_y', nrow=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.title.y = element_blank()) +
  labs(x = "Number of annual time series data")
