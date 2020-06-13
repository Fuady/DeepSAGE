########################################################################################################
##                                          Second simulation:                                        ##
##                       Evaluate the null distribution of the test statistics                        ##
##                          for the hypothesis of no association between                              ##
##                                 gene expression and motor score                                    ##
##                                             31 -05- 2020                                           ##
########################################################################################################

# Libraries
# --------------------------------------------------------------------------------------------
library(dplyr)
library(MASS)
library(doParallel)
library(foreach)

# Functions
count_nb_gen <- function(xb, theta) {
  # Function to generate count data that follow negative binomial distribution
  # Input:
  # xb: design matrix times regression coefficient vector
  # theta: dispersion parameter
  # Output:
  # C1: Count data
  XB.tilde = (1/theta) * exp(xb)
  
  C1 <- matrix(NA, nrow = length(XB.tilde), ncol = 1)
  
  for (i in 1:length(XB.tilde)) {
    C1[i, ] <- rnbinom(1, size = XB.tilde[i], prob = (1/theta)/(1 + (1/theta)))
  }
  return(as.numeric(C1))
}

simul_log_nb <- function(n1, n2, n3, sim_sigmau, sim_sigmae, sim_sigmao, sim_beta0, 
                         sim_beta1, sim_beta2, sim_beta3, sim_theta) {
  # Function to generate log of multivariate count data 
  # Input: n1, n2, n3: totalsample of S1, S2, and S3, respectively 
  # sim_sigmau, sim_sigmae, sim_sigmao: Initial value of sigmau, sigmae, sigmao 
  # sim_beta0, sim_beta1, sim_beta2, sim_beta3: initial value of regression coefficients 
  # sim_theta: Initial value for dispersion parameter 
  # Output: 
  # X1,X2: motor score of S1 at time point 1 and 2, respectively 
  # y1,y2: DeepSAGE of S1 at time point 1 and 2, respectively
  # Xunpair1, Xunpair2: Motor score of S2 and S3, respectively 
  # yunpair1, yunpair2: DeepSAGE of S2 and S3, respectively
  
  N = n1 + n2 + n3
  sim_ms <- mvtnorm::rmvnorm(n = N, mean = c(0, 0), sigma = matrix(c(1, 0.9, 0.9, 
                                                                     1), ncol = 2))
  
  X1 <- as.matrix(sim_ms[1:n1, 1], ncol = 1)
  X2 <- as.matrix(sim_ms[1:n1, 2], ncol = 1)
  
  big_s_L1 <- matrix(c(sim_sigmau^2 + sim_sigmae^2, 0.8 * sqrt(sim_sigmau^2 + sim_sigmae^2) * 
                         sqrt(sim_sigmau^2 + sim_sigmae^2 + sim_sigmao^2), 0.8 * sqrt(sim_sigmau^2 + 
                                                                                        sim_sigmae^2) * sqrt(sim_sigmau^2 + sim_sigmae^2 + sim_sigmao^2), sim_sigmau^2 + 
                         sim_sigmae^2 + sim_sigmao^2), ncol = 2)
  
  U_L1 <- mvtnorm::rmvnorm(n = n1, mean = c(0, 0), sigma = big_s_L1)
  xb_datx1 <- sim_beta0 + sim_beta2 * X1 + U_L1[, 1]
  xb_datx2 <- (sim_beta0 + sim_beta1) + (sim_beta2 + sim_beta3) * X2 + U_L1[, 2]
  count_t1_L1 <- count_nb_gen(xb = xb_datx1, theta = sim_theta)
  count_t2_L1 <- count_nb_gen(xb = xb_datx2, theta = sim_theta)
  y1 <- as.matrix(log(count_t1_L1 + 1), ncol = 1)
  y2 <- as.matrix(log(count_t2_L1 + 1), ncol = 1)
  
  Xunpair1 <- as.matrix(sim_ms[(n1 + 1):(n1 + n2), 1], ncol = 1)
  U_L2 <- rnorm(n = n2, mean = 0, sd = sqrt(sim_sigmau^2 + sim_sigmae^2))
  xb_datx11 <- sim_beta0 + sim_beta2 * Xunpair1 + U_L2
  count_t1_L2 <- count_nb_gen(xb = xb_datx11, theta = sim_theta)
  yunpair1 <- as.matrix(log(count_t1_L2 + 1), ncol = 1)
  
  Xunpair2 <- as.matrix(sim_ms[(n1 + n2 + 1):(n1 + n2 + n3), 2], ncol = 1)
  U_L3 <- rnorm(n = n3, mean = 0, sd = sqrt(sim_sigmau^2 + sim_sigmae^2 + sim_sigmao^2))
  xb_datx22 <- (sim_beta0 + sim_beta1) + (sim_beta2 + sim_beta3) * Xunpair2 + U_L3
  count_t2_L3 <- count_nb_gen(xb = xb_datx22, theta = sim_theta)
  yunpair2 <- as.matrix(log(count_t2_L3 + 1), ncol = 1)
  
  return(list(X1 = as.matrix(X1, ncol = 1), X2 = as.matrix(X2, ncol = 1), y1 = as.matrix(y1, 
                                                                                         ncol = 1), y2 = as.matrix(y2, ncol = 1), Xunpair1 = as.matrix(Xunpair1, ncol = 1), 
              Xunpair2 = as.matrix(Xunpair2, ncol = 1), yunpair1 = as.matrix(yunpair1, 
                                                                             ncol = 1), yunpair2 = as.matrix(yunpair2, ncol = 1)))
}

loglik_fullmodel <- function(theta, y1, y2, X1, X2, yunpair1, yunpair2, Xunpair1, 
                             Xunpair2, sigmaomega) {
  # Function to compute full model log-likelihood 
  # Input: 
  # theta: initial parameters
  # y1: log-transformed count data of pair deepSAGE at time point 1 (subsample S1)
  # y2: log-transformed count data of pair deepSAGE at time point 2 (subsample S1)
  # yunpair1: log-transformed count data of unpair deepSAGE at time point 1 (subsample S2) 
  # yunpair2: log-transformed count data of unpair deepSAGE at time point 2 (subsample S3) 
  # X1: the motor score of pair DeepSAGE at time point 1 
  # X2: the motor score of pair DeepSAGE at time point 2 
  # Xunpair1: the motor score of unpair DeepSAGE at time point 1 
  # Xunpair2: the motor score of unpair DeepSAGE at time point 2 
  # sigmaomega: estimated standard deviation of measurement error
  # Output: 
  # a: log-likelihood of full model
  print(theta)
  betat1 <- c(theta[1], theta[2])
  betat2 <- c(theta[3], theta[4])
  sigu <- exp(theta[5])
  sigeps <- exp(theta[6])
  sigom2 <- sigmaomega^2
  
  if (sigu^2 <= 0 | sigeps^2 <= 0) 
    return(NA)
  n <- length(y1)
  e1 <- y1 - cbind(1, X1) %*% betat1
  e2 <- y2 - cbind(1, X2) %*% betat2
  sig1 <- sqrt(sigu^2 + sigeps^2)
  sig2 <- sqrt(sigu^2 + sigeps^2 + sigom2)
  rho <- sigu^2/(sig1 * sig2)
  logl <- -n * (log(sig1) + log(sig2) + 0.5 * log(1 - rho^2)) - (0.5/(1 - rho^2)) * 
    (sum((t(e1) %*% e1))/sig1^2 + sum((t(e2) %*% e2))/sig2^2 - 2 * rho * sum((t(e2) %*% 
                                                                                e1))/(sig1 * sig2))
  
  loglunpair1 <- sum(dnorm(yunpair1, mean = cbind(1, Xunpair1) %*% betat1, sd = sig1, 
                           log = TRUE))
  loglunpair2 <- sum(dnorm(yunpair2, mean = cbind(1, Xunpair2) %*% betat2, sd = sig2, 
                           log = TRUE))
  
  a <- logl + loglunpair1 + loglunpair2
  return(a)
}

loglik_redmodel <- function(theta, y1, y2, X1, X2, yunpair1, yunpair2, Xunpair1, Xunpair2, 
                          sigmaomega) {
  # Function to compute reduced model log-likelihood (full model without beta2)
  # Input: 
  # theta: initial parameters
  # y1: log-transformed count data of pair deepSAGE at time point 1 (subsample S1)
  # y2: log-transformed count data of pair deepSAGE at time point 2 (subsample S1)
  # yunpair1: log-transformed count data of unpair deepSAGE at time point 1 (subsample S2) 
  # yunpair2: log-transformed count data of unpair deepSAGE at time point 2 (subsample S3) 
  # X1: the motor score of pair DeepSAGE at time point 1 
  # X2: the motor score of pair DeepSAGE at time point 2 
  # Xunpair1: the motor score of unpair DeepSAGE at time point 1 
  # Xunpair2: the motor score of unpair DeepSAGE at time point 2 
  # sigmaomega: estimated standard deviation of measurement error
  # Output: 
  # a: log-likelihood of full model
  print(theta)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- 0
  beta3 <- theta[3]
  betat1 <- c(beta0, beta2)
  betat2 <- c((beta0 + beta1), (beta2 + beta3))
  sigu <- exp(theta[4])
  sigeps <- exp(theta[5])
  sigom2 <- sigmaomega^2
  
  if (sigu^2 <= 0 | sigeps^2 <= 0) 
    return(NA)
  n <- length(y1) 
  e1 <- y1 - cbind(1, X1) %*% betat1
  e2 <- y2 - cbind(1, X2) %*% betat2
  sig1 <- sqrt(sigu^2 + sigeps^2)
  sig2 <- sqrt(sigu^2 + sigeps^2 + sigom2)
  rho <- sigu^2/(sig1 * sig2)
  logl <- -n * (log(sig1) + log(sig2) + 0.5 * log(1 - rho^2)) - (0.5/(1 - rho^2)) * 
    (sum((t(e1) %*% e1))/sig1^2 + sum((t(e2) %*% e2))/sig2^2 - 2 * rho * sum((t(e2) %*% 
                                                                                e1))/(sig1 * sig2))
  
  loglunpair1 <- sum(dnorm(yunpair1, mean = cbind(1, Xunpair1) %*% betat1, sd = sig1, 
                           log = TRUE))
  loglunpair2 <- sum(dnorm(yunpair2, mean = cbind(1, Xunpair2) %*% betat2, sd = sig2, 
                           log = TRUE))
  
  a <- logl + loglunpair1 + loglunpair2
  return(a) 
}

