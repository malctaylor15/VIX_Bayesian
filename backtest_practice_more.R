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
  sigdraws<-sqrt(nusq/rchisq(nsim,nureg))    # Draw all the sigmas
  bdraws <-rmvnorm(nsim,sigma=xpxinv)  # Draw the betas with the
  bdraws <-bdraws*sigdraws+rep(bhat,each=nsim)
  
  
  efdraws <- rnorm(nsim,mean=0,sd=1) # Draw the noise? 
  efdraws <- efdraws * sigdraws  	# Same trick, as for betas
  
  
  lvf     <- matrix(0,ncol=3,nrow=nsim)	# draws of lv
  lvf[,1] <- bdraws %*% c(1,lv[numb_obs],lv[(numb_obs-1)],lv[(numb_obs-2)]) + efdraws 
  
  evix <-exp(lvf) 	# Make the draws of vix itself
  
  
  # See where real values fall on empirical distribution 
  
  emp1 <- ecdf(evix[,1])(vix$VIX[(numb_obs+1)]/100)
  quant1 <- quantile(evix[,1], c(0.05,0.25,0.5,0.75,0.95))
  mean1 <- mean(evix)
  final_return <- c(emp1, quant1, mean1)
  return(final_return)
 
} 

# Pre allocate 
results <- data.frame(matrix(0,ncol = 7, nrow = 1188))
names(results)[c(1,4,7)] <- c("Empirical location","Median", "Mean")

for (i in 1:1188){
  dataset <- vix1[i:(numb_obs+2+i),]
  results[i,] <- distibution_func(dataset, numb_obs, nsim)
}


# Diagonistics - change for each column 
hist(results[,1], breaks = 40)
vixdates1 <- vixdates[204:1391]
plot(vixdates1, results[,1])

dist1 <- results[,7] - results[,1]
hist(dist1, breaks = 40)
ecdf(dist1)(0) # How often model mean under estimates historical returns 

# Look at the distribution of the results that are greater than the mean 
dist2 <- dist1[dist1>0]
mean(dist2)
median(dist2)
hist(dist2, breaks = 30)

# Look at the distribution of the results that are LOWER than the mean 
dist3 <- dist1[dist1<0]
mean(dist3)
median(dist3)
hist(dist3, breaks = 30)

# Squared difference => Chi Square 
dist <- dist1^2
hist(dist, breaks= 40)

makeHistChi <- function(x, df1, color = "blue", title = "Histogram"){
  h<-hist(x,main=title, breaks = 40) 
  xfit<-seq(min(x),max(x),length=40) 
  yfit<-dchisq(xfit,df = df1) 
  yfit <- yfit*diff(h$mids[1:2])*length(x) 
  lines(xfit, yfit, col=color, lwd=2)
}
# Fit Chi Sq distribution to data, not sure how to pick k 
# Check wiki and try different things 

# Try the mean = k 
makeHistChi(dist, df1 = mean(dist))

# Try the variance 
makeHistChi(dist, df1 = (sd(dist))^2)

# Try equation for median = f(k) 
# Found solution on wolfram = 0.5507 
makeHistChi(dist, df1 = 0.5507)

# Next Step- add chi sq goodness of fit 


