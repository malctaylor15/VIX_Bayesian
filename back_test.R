# Trading Strategy 

setwd("C:/Users/board/Desktop/Kaggle/VIX_Bayesian")
library(mvtnorm)

vix<-read.csv("VIX-weekly.csv")
vixdates<-as.Date(vix$Date,"%m/%d/%y")
numb_obs <- 151
nsim  <-40000

# Create Log Returns
vixw<-vix$VIX[1:numb_obs]
lv<-log(vixw/100)

# Prepare linear regression parts 
TT<-numb_obs - 3
xmat  <-cbind(rep(1,TT),lv[3:(numb_obs-1)],lv[2:(numb_obs-2)],lv[1:(numb_obs-3)])
xpxinv<-solve(t(xmat)%*%xmat)
bhat  <-xpxinv%*%t(xmat)%*%lv[4:numb_obs]
fitted<-xmat%*%bhat
resid <-lv[4:numb_obs]-fitted

nureg <- TT - 4
nusq  <- sum(resid^2)

# Start the Monte Carlo Simulation 
sigdraws<-sqrt(nusq/rchisq(nsim,nureg))		# Draw all the sigmas
bdraws <-rmvnorm(nsim,sigma=xpxinv)  # Draw the betas with the
bdraws <-bdraws*sigdraws+rep(bhat,each=nsim)


efdraws <- matrix(rnorm(nsim*3),ncol=3) # Draw the noise? 
efdraws <- efdraws * sigdraws  	# Same trick, as for betas


lvf     <- matrix(0,ncol=3,nrow=nsim)	# draws of lv
lvf[,1] <- bdraws %*% c(1,lv[numb_obs],lv[(numb_obs-1)],lv[(numb_obs-2)]) + efdraws[,1] 

lvf[,2] <- bdraws[,c(1,3,4)]%*%c(1,lv[numb_obs],lv[(numb_obs-1)])+bdraws[,2]*lvf[,1]
lvf[,2] <- lvf[,2] + efdraws[,2] 

lvf[,3] <- bdraws[,c(1,4)]%*%c(1,lv[numb_obs])+bdraws[,2]*lvf[,2]+ bdraws[,3]* lvf[,1]
lvf[,3]  <- lvf[,3] + efdraws[,3] 

evix <-exp(lvf) 	# Make the draws of vix itself


# See where real values fall on empirical distribution 


ecdf(evix[,1])(vix$VIX[152]/100)
ecdf(evix[,2])(vix$VIX[153]/100)
ecdf(evix[,3])(vix$VIX[154]/100)

t(round(apply(evix,2,quantile,c(0.05,0.25,0.5,0.75,0.95)),3))
apply(evix,2,mean)
apply(evix,2,ecdf)

