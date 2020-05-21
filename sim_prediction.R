rm(list=ls())

########################################################################################################
##	                Statistical Method for Modelling Sequencing Data from                             ## 
##                   Different Technologies in Longitudinal Studies with                              ##
##                             Application to Huntington Disease                                      ##
##   						                          2020 - 04 - 28			                               		      ##
########################################################################################################

# libraries --------------------------------------------------------------------------------------------
library(dplyr)
library(MASS)

# functions --------------------------------------------------------------------------------------------
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


n1 = 53 
ns = 19
n2 = 68
n3 = 17

sim_beta0 <- 5.308087
sim_beta1 <- -0.9917194
sim_beta2 <- 0
sim_beta3 <- 0

sim_sigmau <- 0.06893209
sim_sigmae <- 0.3738808
sim_sigmao <- 0.2216958

sim_theta <- 0.15

# 5300, 1900, 6800, 1700
N = n1 + n2 + n3



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


foreach(r = 1:10000) %dopar% {
  print(r)
  set.seed(r)
  
  dat_sim <- simul_log_nb(n1 = n1, n2 = n2, n3 = n3, 
                          sim_sigmau = sim_sigmau, sim_sigmae = sim_sigmae, sim_sigmao = sim_sigmao, 
                          sim_beta0 =  sim_beta0, sim_beta1 = sim_beta1, sim_beta2 = sim_beta2, sim_beta3 = sim_beta3,   
                          sim_theta = sim_theta)
  
  k <- 1
  
  eps <- .Machine$double.eps
  
  optimbivar <- try(optim(rep(0.1,6), loglik_fullmodel, method="BFGS",   
                          y1=dat_sim$y1,
                          y2=dat_sim$y2,
                          X1=dat_sim$X1,
                          X2=dat_sim$X2,
                          yunpair1 = dat_sim$datx11, 
                          yunpair2 = dat_sim$datx22, 
                          Xunpair1 = dat_sim$Xunpair1, 
                          Xunpair2 = dat_sim$Xunpair2,
                          sigmaomega = dat_sim$sim_sigmao,
                          control = list(fnscale=-1,parscale=rep(0.1,6),maxit=500), hessian = TRUE))
  
  optimbivar_ms5 <- try(optim(rep(0.1,5), bivar.fun5.ms, method="BFGS",   
                              y1=y1[,k],
                              y2=y2[,k],
                              X1=X1[,k],
                              X2=X2[,k],
                              yunpair1 = datx11[,k], 
                              yunpair2 = datx22[,k] , 
                              Xunpair1 = Xunpair1[,k], 
                              Xunpair2 = Xunpair2[,k],
                              sigmaomega = sim_sigmao,
                              control = list(fnscale=-1,parscale=rep(0.1,5),maxit=500), hessian = TRUE))
  
  
  beta0 <- optimbivar$par[1]
  beta2 <- optimbivar$par[2]
  beta1 <- optimbivar$par[3]-optimbivar$par[1]
  beta3 <- optimbivar$par[4]-optimbivar$par[2]
  beta23 <- optimbivar$par[4]
  sigmau <- exp(optimbivar$par[5])
  sigmaeps <- exp(optimbivar$par[6])
  converge <- optimbivar$convergence
  
  hessi <- optimbivar$hessian
  
  if(length(which(eigen(-hessi)$values<="1e-10"))>0){  se <- sqrt(diag(MASS::ginv(-hessi + 0.005*diag(ncol(hessi)),tol = 1e-50)))  
  } else { se <- sqrt(diag(MASS::ginv(-hessi,tol = 1e-50))) }
  
  if(length(which(eigen(-hessi)$values<="1e-10"))>0){  sebeta1 <- sqrt(MASS::ginv(-hessi+ 0.005*diag(ncol(hessi)),tol=1e-50)[1,1] + MASS::ginv( -hessi+ 0.005*diag(ncol(hessi)),tol=1e-50)[3,3] - 2*MASS::ginv( -hessi+ 0.005*diag(ncol(hessi)),tol=1e-50)[1,3])  
  } else { sebeta1 <- sqrt(MASS::ginv(-hessi,tol=1e-50)[1,1] + MASS::ginv( -hessi,tol=1e-50)[3,3] - 2*MASS::ginv( -hessi,tol=1e-50)[1,3]) }
  
  if(length(which(eigen(-hessi)$values<="1e-10"))>0){  sebeta3 <- sqrt(MASS::ginv(-hessi+ 0.005*diag(ncol(hessi)),tol=1e-50)[2,2] + MASS::ginv( -hessi+ 0.005*diag(ncol(hessi)),tol=1e-50)[4,4] - 2*MASS::ginv( -hessi+ 0.005*diag(ncol(hessi)),tol=1e-50)[2,4])  
  } else { sebeta3 <- sqrt(MASS::ginv(-hessi,tol=1e-50)[2,2] + MASS::ginv( -hessi,tol=1e-50)[4,4] + 2*MASS::ginv( -hessi,tol=1e-50)[2,4]) }
  
  sebeta0 <- se[1]
  sebeta2 <- se[2]
  sebeta23 <- se[4]
  
  
  pval_b20 <- pchisq(-2*(optimbivar_ms5$value-optimbivar$value),df=1,lower.tail=FALSE)
  #cor_L1 <- cor(sim_dat_L1$RNAseq_t1, sim_dat_L1$RNAseq_t2)
  #cor_L2 <- cor(sim_dat_L2$RNAseq_t1, sim_dat_L2$RNAseq_t2)
  
  tempo <- data.frame(gene= c("gene_1"), position = 1, beta0 = beta0, beta1 = beta1, beta2 = beta2, beta3= beta3,
                      sigmau = sigmau, sigmae = sigmaeps,
                      pvalb20 = pval_b20, llik_ms6 = optimbivar$value,
                      llik_ms5 = optimbivar_ms5$value) 
  write.csv(x=tempo, file=paste("sim_deepsage_ccr2_0_15_",r,".csv", sep=""))
  
}