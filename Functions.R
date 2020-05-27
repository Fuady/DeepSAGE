


rm(list=ls())


simul_log_nb <- function(n1, n2, n3, sim_sigmau, sim_sigmae, sim_sigmao, sim_beta0, sim_beta1, sim_beta2, sim_beta3, sim_theta ){
  #   Function to generate log of multivariate count data
  # Input:
  #  n1, n2, n3: total sample of S1, S2, and S3, respectively
  #  sim_sigmau, sim_sigmae, sim_sigmao: Initial value of sigmau, sigmae, sigmao
  #  sim_beta0, sim_beta1, sim_beta2, sim_beta3: initial value of regression coefficients
  #  sim_theta: Initial value for dispersion parameter  
  # Output:
  # X1,X2: motor score of S1 at time point 1 and 2, respectively
  # y1,y2: DeepSAGE of S1 at time point 1 and 2, respectively 
  # Xunpair1, Xunpair2: motor score of S2 and S3, respectively
  # yunpair1, yunpair2: DeepSAGE of S2 and S3, respectively
  
  N = n1 + n2 + n3
  sim_ms <- mvtnorm::rmvnorm(n = N, mean = c(0,0), sigma = matrix(c(1,0.9,0.9,1), ncol=2))
  
  X1 <- as.matrix(sim_ms[1:n1,1], ncol=1)
  X2 <- as.matrix(sim_ms[1:n1,2], ncol=1)
  
  big_s_L1 <- matrix(c(sim_sigmau^2 + sim_sigmae^2, 
                       0.8*sqrt(sim_sigmau^2 + sim_sigmae^2)*sqrt(sim_sigmau^2 + sim_sigmae^2+ sim_sigmao^2), 
                       0.8*sqrt(sim_sigmau^2 + sim_sigmae^2)*sqrt(sim_sigmau^2 + sim_sigmae^2+ sim_sigmao^2), 
                       sim_sigmau^2 + sim_sigmae^2 + sim_sigmao^2),ncol=2)
  
  U_L1 <- mvtnorm::rmvnorm(n = n1, mean = c(0,0), sigma = big_s_L1)
  xb_datx1 <- sim_beta0 + sim_beta2 * X1 + U_L1[,1]
  xb_datx2 <- (sim_beta0 + sim_beta1) + (sim_beta2 + sim_beta3) * X2 + U_L1[,2]
  count_t1_L1 <- count_nb_gen(xb = xb_datx1, theta = sim_theta)
  count_t2_L1 <- count_nb_gen(xb = xb_datx2, theta = sim_theta)
  y1 <- as.matrix(log(count_t1_L1 + 1),ncol=1)
  y2 <- as.matrix(log(count_t2_L1 + 1),ncol=1)
  
  Xunpair1 <- as.matrix(sim_ms[(n1+1):(n1+n2),1],ncol=1)
  U_L2 <- rnorm(n=n2, mean = 0, sd = sqrt(sim_sigmau^2 + sim_sigmae^2))
  xb_datx11 <- sim_beta0 + sim_beta2 * Xunpair1 + U_L2
  count_t1_L2 <- count_nb_gen(xb = xb_datx11, theta = sim_theta)
  yunpair1 <- as.matrix(log(count_t1_L2 + 1),ncol=1)
  
  Xunpair2 <- as.matrix(sim_ms[(n1+n2+1):(n1+n2+n3),2],ncol=1)
  U_L3 <- rnorm(n = n3, mean = 0, sd = sqrt(sim_sigmau^2 + sim_sigmae^2 + sim_sigmao^2))
  xb_datx22 <- (sim_beta0 + sim_beta1) + (sim_beta2 + sim_beta3) * Xunpair2 + U_L3
  count_t2_L3 <- count_nb_gen(xb = xb_datx22, theta = sim_theta)
  yunpair2 <- as.matrix(log(count_t2_L3 + 1),ncol=1)
  
  return(list(X1 = as.matrix(X1,ncol=1), X2 = as.matrix(X2,ncol=1), 
              y1 = as.matrix(y1,ncol=1), y2 = as.matrix(y2,ncol=1), 
              Xunpair1 = as.matrix(Xunpair1,ncol=1), Xunpair2 = as.matrix(Xunpair2,ncol=1),
              yunpair1 = as.matrix(yunpair1,ncol=1), yunpair2 = as.matrix(yunpair2,ncol=1)))
}

count_nb_gen <- function(xb, theta){
  # Function to generate count data using negative binomial (NB)
  # Input: 
  #       xb : mean of NB
  #       theta: dispersion parameter
  # Output:
  #       C1: count data
  
  XB.tilde = (1/theta)*exp(xb)
  
  C <- matrix(NA,nrow=nrow(XB.tilde),ncol=1)
  
  for (i in 1:nrow(XB.tilde)){
    C[i,] <- rnbinom(1,size = XB.tilde[i,],prob = (1/theta)/(1+(1/theta)))
  }
  return(as.numeric(C))
}

loglik_fullmodel <- function(theta, y1,y2, X1,X2,yunpair1,yunpair2,Xunpair1,Xunpair2,sigmaomega){
  # Function to compute full model log-likelihood
  # Input:
  #       theta: initial parameters
  #       y1: log-transformed count data of pair deepSAGE at time point 1 (subsample S1)
  #       y2: log-transformed count data of pair deepSAGE at time point 2 (subsample S1)
  #       yunpair1: log-transformed count data of unpair deepSAGE at time point 1 (subsample S2)
  #       yunpair2: log-transformed count data of unpair deepSAGE at time point 2 (subsample S3)
  #       X1: the motor score of pair DeepSAGE at time point 1
  #       X2: the motor score of pair DeepSAGE at time point 2
  #       Xunpair1: the motor score of unpair DeepSAGE at time point 1
  #       Xunpair2: the motor score of unpair DeepSAGE at time point 2
  #       sigmaomega: estimated standard deviation of measurement error
  # Output:
  #       a: log-likelihood of full model      
  print(theta) 
  betat1  <- c(theta[1], theta[2])
  betat2  <- c(theta[3], theta[4])
  sigu   <- exp(theta[5])
  sigeps <- exp(theta[6])
  sigom2  <- sigmaomega^2
  
  if (sigu^2 <= 0 | sigeps^2 <= 0) return(NA)
  n <- length(y1) 
  e1 <- y1 - cbind(1,X1)%*%betat1                                  
  e2 <- y2 - cbind(1,X2)%*%betat2
  sig1 <- sqrt(sigu^2 + sigeps^2)
  sig2 <- sqrt(sigu^2 + sigeps^2 + sigom2)
  rho <- sigu^2/(sig1*sig2)  
  logl <-  -n*(log(sig1) + log(sig2) + 0.5*log(1-rho^2)) - (0.5/(1-rho^2))*(
    sum((t(e1)%*%e1))/sig1^2 + sum((t(e2)%*%e2))/sig2^2 - 2*rho*sum((t(e2)%*%e1))/(sig1*sig2) )
  
  loglunpair1 <- sum(dnorm(yunpair1, mean = cbind(1,Xunpair1) %*% betat1, sd = sig1, log = TRUE))
  loglunpair2 <- sum(dnorm(yunpair2, mean = cbind(1,Xunpair2) %*% betat2, sd = sig2, log = TRUE))
  
  a <- logl + loglunpair1 + loglunpair2
  return(a) 
}


count_nb_gen <- function(xb, theta){
  XB.tilde = (1/theta)*exp(xb)

  C1 <- matrix(NA,nrow=length(XB.tilde),ncol=1)

  for (i in 1:length(XB.tilde)){
    C1[i,] <- rnbinom(1,size = XB.tilde[i],prob = (1/theta)/(1+(1/theta)))
  }
  return(as.numeric(C1))
}

sim.log.nb <- function(N, theta, mean1, sd1, mean2, sd2){
  U <- rnorm(n = N, mean = 0, sd = 0.9)
  mu1 = rnorm(n = N, mean = mean1, sd = sd1)
  mu2 = rnorm(n = N, mean = mean2, sd = sd2)
  
  count_t1_L1 <- count_nb_gen(xb = mu1+U, theta = theta)
  count_t2_L1 <- count_nb_gen(xb = mu2+U, theta = theta)
  deepsage.sim <- log(count_t1_L1 + 1)
  rnaseq.sim <- log(count_t2_L1 + 1)
  
  data.log.nb <- data.frame(log.SAGE = deepsage.sim, log.RNA = rnaseq.sim)
  
  return(data.log.nb)
}

prediction.lm <- function(n, ntrain, data.log.nb){
  split.sample <- sample(1:n, ntrain, replace = FALSE)
  
  train_log.data <- data.log.nb[split.sample,]
  test_log.data <- data.log.nb[-split.sample,]
  
  fit <- lm(log.SAGE~ log.RNA, data = train_log.data)
  actual.SAGE <- test_log.data$log.SAGE
  test_log.data$log.SAGE <- NA
  predicted.SAGE <- predict(fit, newdata = test_log.data)
  
  cor.actual <- cor(actual.SAGE, test_log.data$log.RNA )
  cor.predicted <- cor(predicted.SAGE, test_log.data$log.RNA )
  
  res <- data.frame(actual = actual.SAGE, predicted = predicted.SAGE, 
                    cor.actual = cor.actual, cor.predicted = cor.predicted,
                    cor.begin = cor(data.log.nb)[1,2])
  return(res)
  
}
#metrics <- function(ntrain, actual.value, predicted.value){
#   A <- actual.value
#   F <- predicted.value
#   e <- A - F
#   
#   R2   = caret::R2(pred = F, obs = A)
#   MAE  = caret::MAE(pred = F, obs = A)        # the mean absolute error (MAE)
#   RMSE = caret::RMSE(pred = F, obs = A)     # the root mean squared error (RMSE)
#   MAPE = mean(abs(e) / A) * 100   # the mean absolute percentage error (MAPE)
#   
#   res <- list(MAE = MAE, RMSE = RMSE, MAPE = MAPE, R2 = R2, ntrain = ntrain)
#   return(res)
# }

bivar.fun5.ms <- function(theta, y1,y2, X1,X2,yunpair1,yunpair2,Xunpair1,Xunpair2,sigmaomega) {
  print(theta) 
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- 0
  beta3 <- theta[3]
  betat1  <- c(beta0,beta2)
  betat2  <- c((beta0+beta1),(beta2+beta3))  
  sigu   <- exp(theta[4])
  sigeps <- exp(theta[5])
  sigom2  <- sigmaomega^2
  
  if (sigu^2 <= 0 | sigeps^2 <= 0) return(NA)
  n <- length(y1) #nrow(X1)
  e1 <- y1 - cbind(1,X1)%*%betat1                                  
  e2 <- y2 - cbind(1,X2)%*%betat2
  sig1 <- sqrt(sigu^2 + sigeps^2)
  sig2 <- sqrt(sigu^2 + sigeps^2 + sigom2)
  rho <- sigu^2/(sig1*sig2)  
  logl <-  -n*(log(sig1) + log(sig2) + 0.5*log(1-rho^2)) - (0.5/(1-rho^2))*(
    sum((t(e1)%*%e1))/sig1^2 + sum((t(e2)%*%e2))/sig2^2 - 2*rho*sum((t(e2)%*%e1))/(sig1*sig2) )
  
  loglunpair1 <- sum(dnorm(yunpair1, mean = cbind(1,Xunpair1) %*% betat1, sd = sig1, log = TRUE))
  loglunpair2 <- sum(dnorm(yunpair2, mean = cbind(1,Xunpair2) %*% betat2, sd = sig2, log = TRUE))
  
  a <- logl + loglunpair1 + loglunpair2
  return(a) # since optim() does minimisation by default.
}
