##########################################################################################

### Simulation code for the Additive Hazards Mixed Model for Cluster Randomized Trials ###

##########################################################################################

###########################################################################################################################
# Main Code

# Dependencies: pracma, Rfast - libraries
#               survMix.R, survPerm.R, Integrals.R - Files containing supporting macros and functions

# Author: Ondrej Blaha, Fan Li

# Creation Date: November 3, 2021 (R version 4.0.5)

# Purpose: This code uses a list of initial parameters described below and 
#          (a) calculates required sample size (number of clusters) for given input parameters
#          (b) generates simulation and parameters specific data
#          (c) estimates treatment effect and other relevant parameters and quantities
#          (d) evaluates type I error and power for both, z-test and t_(n-1)-test.

# Required Initial Parameters (all scalars):
#          lambda0: Baseline hazard
#          beta0: Treatment effect
#          theta0: Clustering effect
#          CV: Coefficient of variation
#          tau: Duration of the study
#          m_bar: Average cluster size
#          nsim_init: Desired number of simulations
#          alpha: Desired type I error for the underlying test of parameters (z-test, t-test)
#          beta: Desired power for the underlying test of parameters  (z-test, t-test)

# Output:
#          Sample size for both types of underlying tests  (z-test, t-test)
#          Treatment effect estimate
#          Empirical type I error and predicted and empirical power for both underlying tests and all
#          proposed bias-correction methods, including our proposed randomization test.

# Example: Supplying the macro the below specified initial values:
#          lambda0=4
#          beta0<-beta_init=1.5
#          theta0=0.5
#          CV=0
#          m_bar=10
#          tau=1
#          nsim_init=100
#          alpha=0.050
#          beta=0.2
#          yields the following results:
#
#          ######################### 
#          ###Initial Parameters:### 
#          ######################### 
#          Baseline Hazard (lambda0):       4 
#          Treatment Effect (beta0)  :      0 
#          Clustering Effect (theta0):      0.5 
#          Coefficient of Variation(CV):    0 
#          Study Duration (tau):            1 
#          Cluster Size (m):                10 
#          Number of Simulations:           100 
#          ############################ 
#          ###Estimated Sample Size:### 
#          ############################ 
#          Sample Size for z-test:          48 
#          Sample Size for t(n-1)-test:     50 
#          ################## 
#          ###z-test Power### 
#          ################## 
#          Negative Hazard Occurance:       0 
#          Number of Simulations Excluded:  0 
#          Parameter Estimate:              0.05395097 
#          Empirical Variance:              0.1753433 
#
#          Predicted:                       0.814 
#          Cai and Zeng:                    0.830 
#          Kauermann and Carroll:           0.820 
#          Mancl and DeRouen:               0.820 
#          Morel, Bokossa, and Neerchal:    0.820 
#          ######################### 
#          ###z-test type-I error### 
#          ######################### 
#          Negative Hazard Occurance:       0 
#          Number of Simulations Excluded:  0 
#          Parameter Estimate:              0.05395097 
#          Empirical Variance:              0.1753433 
#          
#          Cai and Zeng:                    0.080 
#          Kauermann and Carroll:           0.080 
#          Mancl and DeRouen:               0.080 
#          Morel, Bokossa, and Neerchal:    0.080 
#          ################## 
#          ###t-test Power### 
#          ################## 
#          Negative Hazard Occurance:       0 
#          Number of Simulations Excluded:  0 
#          Parameter Estimate:              0.05395097 
#          Empirical Variance:              0.1753433 
#          
#          Predicted:                       0.814 
#          Cai and Zeng:                    0.910 
#          Kauermann and Carroll:           0.890 
#          Mancl and DeRouen:               0.890 
#          Morel, Bokossa, and Neerchal:    0.890 
#          Permutation:                     0.890 
#          ######################### 
#          ###t-test type-I error### 
#          ######################### 
#          Negative Hazard Occurance:       0 
#          Number of Simulations Excluded:  0 
#          Parameter Estimate:              0.05395097 
#          Empirical Variance:              0.1753433 
#          
#          Cai and Zeng:                    0.040 
#          Kauermann and Carroll:           0.030 
#          Mancl and DeRouen:               0.030 
#          Morel, Bokossa, and Neerchal:    0.030 
#          Permutation:                     0.030 
#          
#          Simulation Time:  64.11 (seconds)

###########################################################################################################################

######################
#Libraries and Macros#
######################
require(pracma)
require(Rfast)
source("survMix.R")
source("survPerm.R")
source("Integrals.R")
##########################
#Fixed Initial Parameters#
##########################
lambda0<-4
beta0<-beta_init<-1.5
theta0<-0.5
CV<-0
m_bar <- 10
tau<-1
nsim_init <- 100
alpha<-0.050
beta<-0.2
###################################################
#Individual quantities for sample size calculation#
###################################################
# Meat of VBV "B"
B=(integral(one,0,tau) + 
     (m_bar*(1+CV^2)-1)*(dblquad(two,0,tau,0,tau,dim=2) + 
                         dblquad(three,0,tau,0,tau,dim=2) +
                         dblquad(four,0,tau,0,tau,dim=2) +
                         dblquad(five,0,tau,0,tau,dim=2)))
# Bread of VBV "V"
V=sqrt(m_bar)*(integral(six,0,tau))

# Sandwich Variance
var_beta <- solve(V)%*%B%*%solve(V)
#########################################################
#Sample Size - for given power, beta0, and var(beta_hat)#
#########################################################
z<-((qnorm(1-beta)+qnorm(1-alpha/2))^2/beta0^2)*(var_beta)
N1<-ceiling(as.numeric(z))+mod(ceiling(as.numeric(z)),2)
powerz<-pnorm(sqrt(N1)*abs(beta0)/sqrt(var_beta)-qnorm(1-alpha/2),mean=0,sd=1)

n_init<-N1
while (n_init<((qt(1-beta,df=n_init-1)-qt(alpha/2,df=n_init-1))^2/beta0^2)*(var_beta))
{n_init<-n_init+1}
N2<-n_init+mod(n_init,2)
powert1<-pt(sqrt(N2)*abs(beta0)/sqrt(var_beta)-qt(1-alpha/2,df=n_init-1),df=N2-1)


start<-proc.time()[[3]]
#############################
#Simulations part I - z-test#
#############################
n <- N1
nsim <- nsim_init
nbeta<-length(beta0)
cl <- 1:n

BETA=matrix(NA,nsim,nbeta)
VARCZ=VARMD=VARKC=VARFG=VARMBN=matrix(NA,nsim,nbeta)
n_cl_power=rep(NA,nsim)
neg_haz_power=0
PVAL = rep(NA, nsim)

# Load permutation matrix
if(choose(n, n/2)<=50000){
  # Enumeration of All Possible Schemes
  emr <- t(combn(n, n / 2))
  R <- dim(emr)[1]
  pmt <- matrix(0, R, n)
  for (r in 1:R){
    pmt[r, emr[r,]] <- 1
  }
} else{
  S=50000
  # Simulation of Permutations
  pmt <- matrix(0, S, n) 
  for (s in 1:S){
    trt <- sample(cl, n / 2)
    pmt[s, trt] <- 1
  }
  pmt <- unique(pmt) # Indicator matrix
  R <- dim(pmt)[1]
}

nk <- matrix(NA,n,nsim)
for(s in 1:nsim){
  set.seed(100+s)
  skip<-0
  
  # Cluster level information
  xi <- rnorm(n,mean=0,sd=sqrt(theta0))
  trt <- rep(0, n)
  trt[sample(1:n,n/2,replace=F)] <- 1
  data <- matrix(nrow = 0, ncol = (3+nbeta))
  
  # Check if Hazard is negative here
  for (k in 1:n){
    if (trt[k]*beta0+xi[k]+lambda0<1e-5) 
    {neg_haz_power=neg_haz_power+1
    skip<-1}
    else{
      if (CV==0) {temp<-m_bar}
      else {
        temp<-ceiling(rgamma(1,shape=CV^-2,rate=m_bar^-1*CV^-2))}
      nk[k,s] <- max(2,temp)
      Xk <- rep(trt[k],nk[k,s])
      Tk <- -log(runif(nk[k,s]))/c(Xk*beta0+xi[k]+lambda0)
      Ck <- runif(nk[k,s])*tau
      data <- rbind(data,cbind(repmat(k,nk[k,s],1),Xk,pmin(Tk,Ck),Tk<=Ck))
    }
  }
  
  ID <- data[,1]
  X <- data[,1+(1:nbeta),drop=FALSE]
  Y <- data[,2+nbeta]
  Delta <- data[,ncol(data)]
  
  res=survMix(ID, X, Y, Delta)
  if (skip==0) {PVAL[s]=survPerm(ID, X_cl=trt, Y, Delta, pmt)} else {PVAL[s]=NA}
  BETA[s,]=res$outbeta[,1]
  VARCZ[s,]=res$outbeta[,2]
  VARKC[s,]=res$outbeta[,3]
  VARMD[s,]=res$outbeta[,4]
  VARFG[s,]=res$outbeta[,5]
  VARMBN[s,]=res$outbeta[,6]
  n_cl_power[s]=res$outbeta[,7]
}

nsim<-nsim-sum(n_cl_power!=n)
BETA<-BETA[n_cl_power==n]
VARCZ_power<-VARCZ[n_cl_power==n]
VARKC_power<-VARKC[n_cl_power==n]
VARMD_power<-VARMD[n_cl_power==n]
VARFG_power<-VARFG[n_cl_power==n]
VARMBN_power<-VARMBN[n_cl_power==n]
################
#Power - z-test#
################
##Cai and Zeng
t_hat<-abs(BETA-0)/sqrt(VARCZ_power)
power_CZ<-sum(t_hat>qnorm(1-alpha/2))/nsim

##Kauermann and Carroll
t_hat<-abs(BETA-0)/sqrt(VARKC_power)
power_KC<-sum(t_hat>qnorm(1-alpha/2))/nsim

##Mancl and DeRouen
t_hat<-abs(BETA-0)/sqrt(VARMD_power)
power_MD<-sum(t_hat>qnorm(1-alpha/2))/nsim

##Fay and Graubard
t_hat<-abs(BETA-0)/sqrt(VARFG_power)
power_FG<-sum(t_hat>qnorm(1-alpha/2))/nsim

##Morel, Bokossa, and Neerchal
t_hat<-abs(BETA-0)/sqrt(VARMBN_power)
power_MBN<-sum(t_hat>qnorm(1-alpha/2))/nsim

##Permutation Test
power_perm <- mean(PVAL[n_cl_power==n]<0.05)
##############################
#Simulations part II - z-test#
##############################
beta0<-0
nsim <- nsim_init
nbeta<-length(beta0)
cl <- 1:n

BETA=matrix(NA,nsim,nbeta)
VARCZ=VARMD=VARKC=VARFG=VARMBN=matrix(NA,nsim,nbeta)
n_cl_typeI=rep(NA,nsim)
neg_haz_typeI=0
PVAL = rep(NA, nsim)

# Load permutation matrix
if(choose(n, n/2)<=50000){
  # Enumeration of All Possible Schemes
  emr <- t(combn(n, n / 2))
  R <- dim(emr)[1]
  pmt <- matrix(0, R, n)
  for (r in 1:R){
    pmt[r, emr[r,]] <- 1
  }
} else{
  S=50000
  # Simulation of Permutations
  pmt <- matrix(0, S, n) 
  for (s in 1:S){
    trt <- sample(cl, n / 2)
    pmt[s, trt] <- 1
  }
  pmt <- unique(pmt) # Indicator matrix
  R <- dim(pmt)[1]
}

nk <- matrix(NA,n,nsim)
for(s in 1:nsim){
  set.seed(100+s)
  skip<-0
  
  # Cluster level information
  xi <- rnorm(n,mean=0,sd=sqrt(theta0))
  trt <- rep(0, n)
  trt[sample(1:n,n/2,replace=F)] <- 1
  data <- matrix(nrow = 0, ncol = (3+nbeta))
  
  # Check if Hazard is negative here
  for (k in 1:n){
    if (trt[k]*beta0+xi[k]+lambda0<1e-5) 
    {neg_haz_typeI=neg_haz_typeI+1
    skip<-1}
    else{
      if (CV==0) {temp<-m_bar}
      else {
        temp<-ceiling(rgamma(1,shape=CV^-2,rate=m_bar^-1*CV^-2))}
      nk[k,s] <- max(2,temp)
      Xk <- rep(trt[k],nk[k,s])
      Tk <- -log(runif(nk[k,s]))/c(Xk*beta0+xi[k]+lambda0)
      Ck <- runif(nk[k,s])*tau
      data <- rbind(data,cbind(repmat(k,nk[k,s],1),Xk,pmin(Tk,Ck),Tk<=Ck))
    }
  }
  
  ID <- data[,1]
  X <- data[,1+(1:nbeta),drop=FALSE]
  Y <- data[,2+nbeta]
  Delta <- data[,ncol(data)]
  
  res=survMix(ID, X, Y, Delta)
  if (skip==0) {PVAL[s]=survPerm(ID, X_cl=trt, Y, Delta, pmt)} else {PVAL[s]=NA}
  BETA[s,]=res$outbeta[,1]
  VARCZ[s,]=res$outbeta[,2]
  VARKC[s,]=res$outbeta[,3]
  VARMD[s,]=res$outbeta[,4]
  VARFG[s,]=res$outbeta[,5]
  VARMBN[s,]=res$outbeta[,6]
  n_cl_typeI[s]=res$outbeta[,7]
}

nsim<-nsim-sum(n_cl_typeI!=n)
BETA<-BETA[n_cl_typeI==n]
VARCZ_typeI<-VARCZ[n_cl_typeI==n]
VARKC_typeI<-VARKC[n_cl_typeI==n]
VARMD_typeI<-VARMD[n_cl_typeI==n]
VARFG_typeI<-VARFG[n_cl_typeI==n]
VARMBN_typeI<-VARMBN[n_cl_typeI==n]
#######################
#Type I-error - z-test#
#######################
##Cai and Zeng
t_hat<-abs(BETA-0)/sqrt(VARCZ_typeI)
typeI_CZ<-sum(t_hat>qnorm(1-alpha/2))/nsim

##Kauermann and Carroll
t_hat<-abs(BETA-0)/sqrt(VARKC_typeI)
typeI_KC<-sum(t_hat>qnorm(1-alpha/2))/nsim

##Mancl and DeRouen
t_hat<-abs(BETA-0)/sqrt(VARMD_typeI)
typeI_MD<-sum(t_hat>qnorm(1-alpha/2))/nsim

##Fay and Graubard
t_hat<-abs(BETA-0)/sqrt(VARFG_typeI)
typeI_FG<-sum(t_hat>qnorm(1-alpha/2))/nsim

##Morel, Bokossa, and Neerchal
t_hat<-abs(BETA-0)/sqrt(VARMBN_typeI)
typeI_MBN<-sum(t_hat>qnorm(1-alpha/2))/nsim

##Permutation Test
typeI_perm <- mean(PVAL[n_cl_typeI==n]<0.05)
#############################
#Simulations part I - t-test#
#############################
beta0<-beta_init
n <- N2
nsim <- nsim_init
nbeta<-length(beta0)
cl <- 1:n

BETA=matrix(NA,nsim,nbeta)
VARCZ=VARMD=VARKC=VARFG=VARMBN=matrix(NA,nsim,nbeta)
n_cl_power=rep(NA,nsim)
neg_haz_power=0
PVAL = rep(NA, nsim)

# Load permutation matrix
if(choose(n, n/2)<=50000){
  # Enumeration of All Possible Schemes
  emr <- t(combn(n, n / 2))
  R <- dim(emr)[1]
  pmt <- matrix(0, R, n)
  for (r in 1:R){
    pmt[r, emr[r,]] <- 1
  }
} else{
  S=50000
  # Simulation of Permutations
  pmt <- matrix(0, S, n) 
  for (s in 1:S){
    trt <- sample(cl, n / 2)
    pmt[s, trt] <- 1
  }
  pmt <- unique(pmt) # Indicator matrix
  R <- dim(pmt)[1]
}

ind <- matrix(NA,n,nsim)
nk <- matrix(NA,n,nsim)
for(s in 1:nsim){
  set.seed(100+s)
  skip<-0
  
  # Cluster level information
  xi <- rnorm(n,mean=0,sd=sqrt(theta0))
  trt <- rep(0, n)
  trt[sample(1:n,n/2,replace=F)] <- 1
  data <- matrix(nrow = 0, ncol = (3+nbeta))
  
  # Check if Hazard is negative here
  for (k in 1:n){
    if (trt[k]*beta0+xi[k]+lambda0<1e-5) 
    {neg_haz_power=neg_haz_power+1
    skip<-1}
    else{
      if (CV==0) {temp<-m_bar}
      else {
        temp<-ceiling(rgamma(1,shape=CV^-2,rate=m_bar^-1*CV^-2))}
      nk[k,s] <- max(2,temp)
      ind[k,s] <- as.numeric(2>temp)
      Xk <- rep(trt[k],nk[k,s])
      Tk <- -log(runif(nk[k,s]))/c(Xk*beta0+xi[k]+lambda0)
      Ck <- runif(nk[k,s])*tau
      data <- rbind(data,cbind(repmat(k,nk[k,s],1),Xk,pmin(Tk,Ck),Tk<=Ck))
    }
  }
  
  ID <- data[,1]
  X <- data[,1+(1:nbeta),drop=FALSE]
  Y <- data[,2+nbeta]
  Delta <- data[,ncol(data)]
  
  res=survMix(ID, X, Y, Delta)
  if (skip==0) {PVAL[s]=survPerm(ID, X_cl=trt, Y, Delta, pmt)} else {PVAL[s]=NA}
  BETA[s,]=res$outbeta[,1]
  VARCZ[s,]=res$outbeta[,2]
  VARKC[s,]=res$outbeta[,3]
  VARMD[s,]=res$outbeta[,4]
  VARFG[s,]=res$outbeta[,5]
  VARMBN[s,]=res$outbeta[,6]
  n_cl_power[s]=res$outbeta[,7]
}

nsim<-nsim-sum(n_cl_power!=n)
BETA<-BETA[n_cl_power==n]
VARCZ_power<-VARCZ[n_cl_power==n]
VARKC_power<-VARKC[n_cl_power==n]
VARMD_power<-VARMD[n_cl_power==n]
VARFG_power<-VARFG[n_cl_power==n]
VARMBN_power<-VARMBN[n_cl_power==n]
################
#Power - t-test#
################
##Cai and Zeng
t_hat<-abs(BETA-0)/sqrt(VARCZ_power)
power1_CZ<-sum(t_hat>qt(1-alpha/2,n-1))/nsim
power2_CZ<-sum(t_hat>qt(1-alpha/2,n-2))/nsim

##Kauermann and Carroll
t_hat<-abs(BETA-0)/sqrt(VARKC_power)
power1_KC<-sum(t_hat>qt(1-alpha/2,n-1))/nsim
power2_KC<-sum(t_hat>qt(1-alpha/2,n-2))/nsim

##Mancl and DeRouen
t_hat<-abs(BETA-0)/sqrt(VARMD_power)
power1_MD<-sum(t_hat>qt(1-alpha/2,n-1))/nsim
power2_MD<-sum(t_hat>qt(1-alpha/2,n-2))/nsim

##Fay and Graubard
t_hat<-abs(BETA-0)/sqrt(VARFG_power)
power1_FG<-sum(t_hat>qt(1-alpha/2,n-1))/nsim
power2_FG<-sum(t_hat>qt(1-alpha/2,n-2))/nsim

##Morel, Bokossa, and Neerchal
t_hat<-abs(BETA-0)/sqrt(VARMBN_power)
power1_MBN<-sum(t_hat>qt(1-alpha/2,n-1))/nsim
power2_MBN<-sum(t_hat>qt(1-alpha/2,n-2))/nsim

##Permutation Test
power_perm <- mean(PVAL[n_cl_power==n]<0.05)
##############################
#Simulations part II - t-test#
##############################
beta0<-0
nsim <- nsim_init
nbeta<-length(beta0)
cl <- 1:n

BETA=matrix(NA,nsim,nbeta)
VARCZ=VARMD=VARKC=VARFG=VARMBN=matrix(NA,nsim,nbeta)
n_cl_typeI=rep(NA,nsim)
neg_haz_typeI=0
PVAL = rep(NA, nsim)

# Load permutation matrix
if(choose(n, n/2)<=50000){
  # Enumeration of All Possible Schemes
  emr <- t(combn(n, n / 2))
  R <- dim(emr)[1]
  pmt <- matrix(0, R, n)
  for (r in 1:R){
    pmt[r, emr[r,]] <- 1
  }
} else{
  S=50000
  # Simulation of Permutations
  pmt <- matrix(0, S, n) 
  for (s in 1:S){
    trt <- sample(cl, n / 2)
    pmt[s, trt] <- 1
  }
  pmt <- unique(pmt) # Indicator matrix
  R <- dim(pmt)[1]
}

ind <- matrix(NA,n,nsim)
nk <- matrix(NA,n,nsim)
for(s in 1:nsim){
  set.seed(100+s)
  skip<-0
  
  # Cluster level information
  xi <- rnorm(n,mean=0,sd=sqrt(theta0))
  trt <- rep(0, n)
  trt[sample(1:n,n/2,replace=F)] <- 1
  data <- matrix(nrow = 0, ncol = (3+nbeta))
  
  # Check if Hazard is negative here
  for (k in 1:n){
    if (trt[k]*beta0+xi[k]+lambda0<1e-5) 
    {neg_haz_typeI=neg_haz_typeI+1
    skip<-1}
    else{
      if (CV==0) {temp<-m_bar}
      else {
        temp<-ceiling(rgamma(1,shape=CV^-2,rate=m_bar^-1*CV^-2))}
      nk[k,s] <- max(2,temp)
      ind[k,s] <- as.numeric(2>temp)
      Xk <- rep(trt[k],nk[k,s])
      Tk <- -log(runif(nk[k,s]))/c(Xk*beta0+xi[k]+lambda0)
      Ck <- runif(nk[k,s])*tau
      data <- rbind(data,cbind(repmat(k,nk[k,s],1),Xk,pmin(Tk,Ck),Tk<=Ck))
    }
  }
  
  ID <- data[,1]
  X <- data[,1+(1:nbeta),drop=FALSE]
  Y <- data[,2+nbeta]
  Delta <- data[,ncol(data)]
  
  res=survMix(ID, X, Y, Delta)
  if (skip==0) {PVAL[s]=survPerm(ID, X_cl=trt, Y, Delta, pmt)} else {PVAL[s]=0}
  BETA[s,]=res$outbeta[,1]
  VARCZ[s,]=res$outbeta[,2]
  VARKC[s,]=res$outbeta[,3]
  VARMD[s,]=res$outbeta[,4]
  VARFG[s,]=res$outbeta[,5]
  VARMBN[s,]=res$outbeta[,6]
  n_cl_typeI[s]=res$outbeta[,7]
}

nsim<-nsim-sum(n_cl_typeI!=n)
BETA<-BETA[n_cl_typeI==n]
VARCZ_typeI<-VARCZ[n_cl_typeI==n]
VARKC_typeI<-VARKC[n_cl_typeI==n]
VARMD_typeI<-VARMD[n_cl_typeI==n]
VARFG_typeI<-VARFG[n_cl_typeI==n]
VARMBN_typeI<-VARMBN[n_cl_typeI==n]
######################
#TypeI-error - t-test#
######################
##Cai and Zeng
t_hat<-abs(BETA-0)/sqrt(VARCZ_typeI)
typeI1_CZ<-sum(t_hat>qt(1-alpha/2,n-1))/nsim
typeI2_CZ<-sum(t_hat>qt(1-alpha/2,n-2))/nsim

##Kauermann and Carroll
t_hat<-abs(BETA-0)/sqrt(VARKC_typeI)
typeI1_KC<-sum(t_hat>qt(1-alpha/2,n-1))/nsim
typeI2_KC<-sum(t_hat>qt(1-alpha/2,n-2))/nsim

##Mancl and DeRouen
t_hat<-abs(BETA-0)/sqrt(VARMD_typeI)
typeI1_MD<-sum(t_hat>qt(1-alpha/2,n-1))/nsim
typeI2_MD<-sum(t_hat>qt(1-alpha/2,n-2))/nsim

##Fay and Graubard
t_hat<-abs(BETA-0)/sqrt(VARFG_typeI)
typeI1_FG<-sum(t_hat>qt(1-alpha/2,n-1))/nsim
typeI2_FG<-sum(t_hat>qt(1-alpha/2,n-2))/nsim

##Morel, Bokossa, and Neerchal
t_hat<-abs(BETA-0)/sqrt(VARMBN_typeI)
typeI1_MBN<-sum(t_hat>qt(1-alpha/2,n-1))/nsim
typeI2_MBN<-sum(t_hat>qt(1-alpha/2,n-2))/nsim

##Permutation Test
typeI_perm <- mean(PVAL[n_cl_typeI==n]<0.05)

##Cluster-Level Summary
neg_haz_typeI2<-neg_haz_typeI
neg_haz_power2<-neg_haz_power
n_cl_typeI2<-n_cl_typeI
n_cl_power2<-n_cl_power

end<-proc.time()[[3]]
Elapsed<-end-start
#########
#Results#
#########
cat(
"#########################","\n",
"###Initial Parameters:###","\n",
"#########################","\n",
"Baseline Hazard (lambda0):      ",lambda0,"\n",
"Treatment Effect (beta0)  :     ",beta0,"\n",
"Clustering Effect (theta0):     ",theta0,"\n",
"Coefficient of Variation(CV):   ",CV,"\n",
"Study Duration (tau):           ",tau,"\n",
"Cluster Size (m):               ",m_bar,"\n",
"Number of Simulations:          ",nsim_init,"\n",
"############################","\n",
"###Estimated Sample Size:###","\n",
"############################","\n",
"Sample Size for z-test:         ",N1,"\n",
"Sample Size for t(n-1)-test:    ",N2,"\n",
"##################","\n",
"###z-test Power###","\n",
"##################","\n",
"Negative Hazard Occurance:      ",neg_haz_power,"\n",
"Number of Simulations Excluded: ",sum(n_cl_power!=n),"\n",
"Parameter Estimate:             ",mean(BETA),"\n",
"Empirical Variance:             ",var(BETA),"\n",
"\n",
"Predicted:                      ",format(round(powerz, digits = 3), nsmall=3),"\n",
"Cai and Zeng:                   ",format(round(power_CZ, digits = 3), nsmall=3),"\n",
"Kauermann and Carroll:          ",format(round(power_KC, digits = 3), nsmall=3),"\n",
"Mancl and DeRouen:              ",format(round(power_MD, digits = 3), nsmall=3),"\n",
"Morel, Bokossa, and Neerchal:   ",format(round(power_MBN, digits = 3), nsmall=3),"\n",
"#########################","\n",
"###z-test type-I error###","\n",
"#########################","\n",
"Negative Hazard Occurance:      ",neg_haz_typeI,"\n",
"Number of Simulations Excluded: ",sum(n_cl_typeI!=n),"\n",
"Parameter Estimate:             ",mean(BETA),"\n",
"Empirical Variance:             ",var(BETA),"\n",
"\n",
"Cai and Zeng:                   ",format(round(typeI_CZ, digits = 3), nsmall=3),"\n",
"Kauermann and Carroll:          ",format(round(typeI_KC, digits = 3), nsmall=3),"\n",
"Mancl and DeRouen:              ",format(round(typeI_MD, digits = 3), nsmall=3),"\n",
"Morel, Bokossa, and Neerchal:   ",format(round(typeI_MBN, digits = 3), nsmall=3),"\n",
"##################","\n",
"###t-test Power###","\n",
"##################","\n",
"Negative Hazard Occurance:      ",neg_haz_power,"\n",
"Number of Simulations Excluded: ",sum(n_cl_power!=n),"\n",
"Parameter Estimate:             ",mean(BETA),"\n",
"Empirical Variance:             ",var(BETA),"\n",
"\n",
"Predicted:                      ",format(round(powert1, digits = 3), nsmall=3),"\n",
"Cai and Zeng:                   ",format(round(power1_CZ, digits = 3), nsmall=3),"\n",
"Kauermann and Carroll:          ",format(round(power1_KC, digits = 3), nsmall=3),"\n",
"Mancl and DeRouen:              ",format(round(power1_MD, digits = 3), nsmall=3),"\n",
"Morel, Bokossa, and Neerchal:   ",format(round(power1_MBN, digits = 3), nsmall=3),"\n",
"Permutation:                    ",format(round(power_perm, digits = 3), nsmall=3),"\n",
"#########################","\n",
"###t-test type-I error###","\n",
"#########################","\n",
"Negative Hazard Occurance:      ",neg_haz_typeI,"\n",
"Number of Simulations Excluded: ",sum(n_cl_typeI!=n),"\n",
"Parameter Estimate:             ",mean(BETA),"\n",
"Empirical Variance:             ",var(BETA),"\n",
"\n",
"Cai and Zeng:                   ",format(round(typeI1_CZ, digits = 3), nsmall=3),"\n",
"Kauermann and Carroll:          ",format(round(typeI1_KC, digits = 3), nsmall=3),"\n",
"Mancl and DeRouen:              ",format(round(typeI1_MD, digits = 3), nsmall=3),"\n",
"Morel, Bokossa, and Neerchal:   ",format(round(typeI1_MBN, digits = 3), nsmall=3),"\n",
"Permutation:                    ",format(round(typeI_perm, digits = 3), nsmall=3),"\n",
"\n",
"Simulation Time: ",Elapsed,"(seconds)")