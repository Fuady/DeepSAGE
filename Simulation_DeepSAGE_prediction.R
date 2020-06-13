########################################################################################################
##                                           First simulation:                                        ##
##                                      DeepSAGE prediction from RNA-Seq                              ##
##                                             31 -05- 2020                                           ##
########################################################################################################

# Libraries --------------------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)

# Functions --------------------------------------------------------------------------------------------
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

sim.log.nb <- function(N, theta, mean1, mean2) {
  # Function to generate log-count data
  # Input:
  # N: Number of sample
  # theta: dispersion parameter
  # mean1: Mean of count data 1
  # mean2: Mean of count data 2
  # Output:
  # data.log.nb: data frame of log-count DeepSAGE and RNA-Seq
  
  U <- rnorm(n = N, mean = 0, sd = 0.9)
  mu1 = mean1
  mu2 = mean2
  
  count_t1_L1 <- count_nb_gen(xb = mu1 + U, theta = theta)
  count_t2_L1 <- count_nb_gen(xb = mu2 + U, theta = theta)
  deepsage.sim <- log(count_t1_L1 + 1)
  rnaseq.sim <- log(count_t2_L1 + 1)
  
  data.log.nb <- data.frame(log.SAGE = deepsage.sim, log.RNA = rnaseq.sim)
  
  return(data.log.nb)
}

prediction.lm <- function(N, ntrain, data.log.nb) {
  # Function to predict DeepSAGE from RNA-Seq
  # Input:
  # N: number of sample
  # ntrain: number of sample for training set
  # data.log.nb: log transform count data of DeepSAGE and RNA-Seq
  # Output:
  # res: data frame of actual and predicted DeepSAGE 
  
  split.sample <- sample(1:N, ntrain, replace = FALSE)
  
  train_log.data <- data.log.nb[split.sample, ]
  test_log.data <- data.log.nb[-split.sample, ]
  
  fit <- lm(log.SAGE ~ log.RNA, data = train_log.data)
  actual.SAGE <- test_log.data$log.SAGE
  test_log.data$log.SAGE <- NA
  predicted.SAGE <- predict(fit, newdata = test_log.data)
  
  res <- data.frame(actual = actual.SAGE, predicted = predicted.SAGE)
  return(res)
  
}
