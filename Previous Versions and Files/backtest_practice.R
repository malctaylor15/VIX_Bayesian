# Trading Strategy 

setwd("C:/Users/board/Desktop/Kaggle/VIX_Bayesian")
library(mvtnorm)

vix1<-read.csv("VIX-weekly-26yr.csv")
vixdates<-as.Date(vix1$Date,"%m/%d/%y")
numb_obs <- 200
nsim  <-50000

distibution_func <- function(vix, numb_obs, nsim ){

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
sigdraws<-sqrt(nusq/rchisq(nsim,nureg))  	# Draw all the sigmas
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

emp1 <- ecdf(evix[,1])(vix$VIX[(numb_obs+1)]/100)
emp2 <- ecdf(evix[,2])(vix$VIX[(numb_obs+2)]/100)
emp3 <- ecdf(evix[,3])(vix$VIX[(numb_obs+3)]/100)
quant1 <- quantile(evix[,1], c(0.05,0.25,0.5,0.75,0.95))
quant2 <- quantile(evix[,2], c(0.05,0.25,0.5,0.75,0.95))
quant3 <- quantile(evix[,3], c(0.05,0.25,0.5,0.75,0.95))
mean1 <- apply(evix,2,mean)
final_return <- c(emp1, emp2, emp3, quant1, quant2, quant3, mean1)
return(final_return)

} 

# Pre allocate 
results <- data.frame(matrix(0,ncol = 21, nrow = 1188))


for (i in 1:1188){
dataset <- vix1[i:(numb_obs+2+i),]
results[i,] <- distibution_func(dataset, numb_obs, nsim)
}


# Diagonistics - change for each column 
hist(results[,1], breaks = 40)
vixdates1 <- vixdates[204:1391]
plot(vixdates1, results[,2])

dist <- results[,21] - results[,1]
hist(dist, breaks= 40)
makeHist <- function(x, color = "blue", title = "Histogram"){
  h<-hist(x,main=title, breaks = 40) 
  xfit<-seq(min(x),max(x),length=40) 
  yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
  yfit <- yfit*diff(h$mids[1:2])*length(x) 
  lines(xfit, yfit, col=color, lwd=2)
}
makeHist(dist)

