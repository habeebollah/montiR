rm(list=ls())


library(rstan)
library(tidyverse)
library(ggdist)
library(here)
library(montiR)
library(rfishbase)
library(ggplot2)
rstan_options(auto_write = T) #To avoid recompilation of unchanged Stan programs
options(mc.cores = parallel::detectCores()) #For execution on a local, multicore CPU with excess RAM

setwd("~/GitHub/montiR/tests") # PC
setwd("~/Documents/GitHub/montiR/tests") # iMac

#load in data
dat <- readr::read_csv(here::here("tests","namibian_hake_data.csv"))
#dat <- montiR::df.namibianCatch
#dat <- montiR::df.eastpacCatch

# clean up data names
#dat[,3] <- dat[,2]/dat[,3]
#colnames(dat) <- c("year", "catch", "cpue")

# quick exploratory plot
dat %>%
  pivot_longer(-year,'variable','value') %>%
  ggplot(aes(year, value, color = variable)) +
  geom_line() +
  facet_wrap(~variable, scales = 'free_y')


# montiR::rGrabber('Merluccius capensis')
# montiR::rGrabber('Merluccius paradoxus')
# montiR::rGrabber("Thunnus albacares")

r <- montiR::rGrabber("Merluccius capensis")[[1]] # grab r median parameter
warmup <- 4000 # number of warmup iterations
iter <- 10000 # total iterations
chains <- 1 # number of chains

dat.list <- list(n_years = nrow(dat),
                    years = dat$year,
                    harvest = dat$catch,
                    cpue = dat$cpue,
                    logr_mean=log(r))

# fit all parameters
fit <- rstan::stan(file ="schaefer_noinput.stan",
                          data = dat.list,
                          iter = iter, warmup = warmup,
                          chains = chains, cores = chains)

fit.r <- rstan::stan(file = "schaefer_rinput.stan",
                   data = dat.list,
                   iter = iter, warmup = warmup,
                   chains = chains, cores = chains)


summary(fit)$summary[c(3:7,122:123),c(1,2,5,7)] # namibianHake [c(2:6,122,123),]; eastpacCatch [c(2:6,117:118),]
summary(fit.r)$summary[c(2:6,121:122),c(1,2,5,7)] # namibianHake [c(2:6,121,122),]; eastpacCatch [c(2:6,116:117),]

cpue_fit <-
  tidybayes::spread_draws(fit, EstCPUE[year], pp_CPUE[year] , EstBt[year], ndraws = 200) # use the tidybayes package to pull out the things I want to plot, in particular the estimated abundance index and the posterior predictive abundance index


# plot fit to cpue
cpue_fit %>%
  ggplot(aes(year, pp_CPUE)) +
  stat_lineribbon(aes(fill = stat(.width)), .width = ppoints(50), alpha = 0.75) +
  geom_point(data = dat,
             aes(year, cpue),
             size = 4,
             color = "tomato")

# plot posterior predictive fits
cpue_fit %>%
  ggplot(aes(year, EstCPUE)) +
  stat_lineribbon(aes(fill = stat(.width)), .width = ppoints(50), alpha = 0.75) +
  geom_point(data = dat,
             aes(year, cpue),
             size = 4,
             color = "tomato")

# plot fit to cpue and posterior predictive. still an unfinished code !
cpue_fit %>% pivot_longer(-c(3,4,5,7), 'parameter', 'value') %>%
  ggplot(aes(tipe, value)) +
  geom_point(data = dat, aes(year, cpue), size = 4, color = "tomato") +
  facet_wrap(pp_CPUE~parameter, scales = 'fixed', nrow=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.title.y = element_blank()) +
  labs(x = "Year")

