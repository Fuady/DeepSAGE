########################################################################################################
##	                Statistical Method for Modelling Sequencing Data from                             ## 
##                   Different Technologies in Longitudinal Studies with                              ##
##                             Application to Huntington Disease                                      ##
##   						                          2020 - 04 - 28			                               		      ##
########################################################################################################

# libraries --------------------------------------------------------------------------------------------
library(dplyr)

# Simulation of DeepSAGE prediction from RNA-Seq

set.seed(12345)


results <- c()
n = 100
ntrain_percent = 0.2
theta=3.5 
mean1=2.15
mean2 = 1.95
sd1 = 0.22
sd2 = 0.17

data.sim.log.nb <- sim.log.nb(N = n, theta = theta, mean1 = mean1, sd1 = sd1, mean2 = mean2, sd2 = sd2)
pred.value <- prediction.lm(n = n, ntrain = n*ntrain_percent, data.log.nb = data.sim.log.nb)

qqplot(pred.value$predicted,pred.value$actual, col=adjustcolor("black"), 
       pch=3, cex=0.7, xlab="Predicted", ylab ="Observed")
abline(coef = c(0,1))


