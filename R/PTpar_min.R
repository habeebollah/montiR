
PTpar_min <- function(inpars, df, mtype=1, OWT=FALSE, currentF = 0.7, weight = 0.5){

  temp <- ifelse(mtype==1, inpars <- inpars[-6], ifelse(mtype==2, inpars <- inpars[-6], inpars))
  p.par <- ifelse(mtype==1, 1, ifelse(mtype==2, 1e-08, inpars[6]))

  K <- exp(inpars[1])
  B0 <- exp(inpars[2])
  r <- exp(inpars[3])
  q <- exp(inpars[4])
  sigma <- exp(inpars[5])
  #p.par <- ifelse(inpars[6]==3, 2, ifelse(inpars[6]==1, 1, 0)) # 1= Schaefer (p=1), 2= Fox (p=0), 3= Pella-Tomlinson (p0-3)
  #p.par <- ifelse(mtype==1, 1, ifelse(mtype==2, 0, 0.5)) # 1= Schaefer, 2= Fox, 3= Pella-Tomlinson

  CPUE <- df$catch/df$effort
  EstCatch <- EstCPUE <- EstBt <- vector(length=nrow(df))

  for (i in 1:nrow(df)) {
    if (i == 1) EstBt[i] <- B0
    if (i>1) EstBt[i] <- max(c(0.01,
                               EstBt[i-1] +  EstBt[i-1] * (r/p.par) * (1 - (EstBt[i-1]/K)^p.par) - df$catch[i-1]
    )
    )
  }
  EstCPUE <-  EstBt * q

  if (OWT==FALSE){
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = sigma, log = TRUE))
  }

  if (OWT=="Depletion"){
    annualFrates <- df[,2]/EstBt # F = catch/biomass = Z-M
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = sigma, log = TRUE)) +
      weight * (tail(annualFrates,1) - currentF)^2
    print(tail(annualFrates,1)) # to check whether the weighting makes the estimates depletion getting closer to the current exploitation rate
  }

  if (OWT=="Biomass"){ # haven't been checked yet
    surveyB <- df[,4]
    nll <- -sum(dlnorm(x= na.omit(CPUE), meanlog = log(na.omit(EstCPUE)), sdlog = sigma, log = TRUE)) +
      -sum(weight * (na.omit(surveyB) - na.omit(EstBt))^2)
  }

  return(nll)
}

# the trajectory and stochastic can be found in Punt & Hilborn (1996) p 44



K <- 1000
B0 <- K
r <- 0.2
q <- 0.00025

### Estimate parameters using optim
inpars <- c(log(K), log(B0), log(r), log(q), 0.1)
mtype=1

temp <- ifelse(mtype==1, inpars <- inpars[-6], ifelse(mtype==2, inpars <- inpars[-6], inpars))
p.par <- ifelse(mtype==1, 1, ifelse(mtype==2, 1e-08, inpars[6]))

inpars
p.par

#fit <- optim(par=inpars,
#             fn=biodyn_min,
#             df=df.goodcontrast, mtype=1,
#             method="Nelder-Mead",
#             OWT=FALSE, currentF = 0.7, weight = 0.5)

fit <- optim(par=inpars,
             fn=Par.min,
             df=df.goodcontrast,
             method="Nelder-Mead",
             OWT=FALSE, currentF = 0.7, weight = 0.5)

#Spar_vals <- data.frame(SPpar = c("K", "B0", "r", "q", "sigma", "p.par"),
#                        init_pars = c(exp(inpars[1:5]), inpars[6]),
#                        fitted_pars = c(exp(fit$par[1:5]),fit$par[6]))

Par_vals <- data.frame(SPpar = c("K", "B0", "r", "q", "sigma"),
                        init_pars = exp(inpars),
                        fitted_pars = exp(fit$par[1:5]))

Par_vals



#plot productivity and density-dependence functions Fig7.4
prodfun <- function(r,Bt,K,p) return((r*Bt/p)*(1-(Bt/K)^p))
densdep <- function(Bt,K,p) return((1/p)*(1-(Bt/K)^p))
r <- 0.75; K <- 1000.0; Bt <- 1:1000
sp <- prodfun(r,Bt,K,1.0)  # Schaefer equivalent
sp0 <- prodfun(r,Bt,K,p=1e-08)  # Fox equivalent
sp3 <- prodfun(r,Bt,K,3) #left skewed production, marine mammal?
parset(plots=c(2,1),margin=c(0.35,0.4,0.1,0.05))
plot1(Bt,sp,type="l",lwd=2,xlab="Stock Size",
      ylab="Surplus Production",maxy=200,defpar=FALSE)
lines(Bt,sp0 * (max(sp)/max(sp0)),lwd=2,col=2,lty=2) # rescale
lines(Bt,sp3*(max(sp)/max(sp3)),lwd=3,col=3,lty=3)   # production
legend(275,100,cex=1.1,lty=1:3,c("p = 1.0 Schaefer","p = 1e-08 Fox",
                                 "p = 3 LeftSkewed"),col=c(1,2,3),lwd=3,bty="n")
plot1(Bt,densdep(Bt,K,p=1),xlab="Stock Size",defpar=FALSE,
      ylab="Density-Dependence",maxy=2.5,lwd=2)
lines(Bt,densdep(Bt,K,1e-08),lwd=2,col=2,lty=2)
lines(Bt,densdep(Bt,K,3),lwd=3,col=3,lty=3)
