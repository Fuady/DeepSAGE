########################################################################################################
##	                Statistical Method for Modelling Sequencing Data from                             ## 
##                   Different Technologies in Longitudinal Studies with                              ##
##                             Application to Huntington Disease                                      ##
##   						                          2020 - 04 - 28			                               		      ##
########################################################################################################

# libraries --------------------------------------------------------------------------------------------
library(dplyr)
library(MASS)
library(doParallel)
library(foreach)

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

N = n1 + n2 + n3


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
#  write.csv(x=tempo, file=paste("sim_deepsage_ccr2_0_15_",r,".csv", sep=""))
  
}
