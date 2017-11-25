# Bayesian Analysis of the VIX

Most recent copy is backtest_practice_more.R which uses VIX-weekly-26yr.csv . 
VIX-weekly-26yr.csv is the date and weekly historical prices of the VIX for the past 26 years- since the inception in 1990. 

The project was a Bayesian Monte Carlo simulation using VIX data. We used an AR(3) model with weekly VIX data to find the posterior distribution of the next VIX return. 

The future returns were distributed given the following equations. 
r_(t+1)=β_0+β_1 r_t+β_2 r_(t-1)+β_3 r_(t-2)+ϵ
Where the β~ N(β_MLE,σ) and σ~1/(χ^2 (#obs-1) ) and ϵ~N(0,1)

In English, this means that the future return is based on the previous 3 time periods. The coefficients for the returns are normally dsitributed with the Maximum likelihood estimate is the mean and the standard deviation of each of the coefficients has an inverted chi-squared distribution.  

More information is available in the Ideas for a trading strategy document. 

Many thanks to Eric Jacquier for my induction to my first Bayesian Monte Carlo simulation and his help constructing the core of this algorithm. 
