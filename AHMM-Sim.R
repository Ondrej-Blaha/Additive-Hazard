##########################################################################################

### Simulation code for the Additive Hazards Mixed Model for Cluster Randomized Trials ###

##########################################################################################

###########################################################################################################################
# Main Code

# Dependencies: pracma, Rfast - libraries
#               survMix.R, survPerm.R, Integrals.R - Files containing supporting macros and functions

# Author: Ondrej Blaha, Fan Li

# Creation Date: Deceber 13, 2021 (R version 4.0.5)

# Purpose: This code uses a list of initial parameters described below and 
#          (a) calculates required sample size (number of clusters) for given input parameters
#          (b) generates simulation and parameters specific data
#          (c) estimates treatment effect and other relevant parameters and quantities
#          (d) evaluates type I error and power for both, z-test and t_(n-1)-test.

# Required Initial Parameters (all scalars):
#          lambda0: Pre-specified baseline hazard
#          delta: Pre-specified treatment effect
#          sigma2: Pre-specified clustering effect
#          CV: Coefficient of variation
#          tau: Duration of the study
#          m_bar: Average cluster size
#          nsim: Desired number of simulations
#          alpha: Desired type I error for the underlying test of parameters (z-test, t-test)
#          beta: Desired power for the underlying test of parameters  (z-test, t-test)

# Output:
#          Sample size for both types of underlying tests  (z-test, t-test)
#          Treatment effect estimate (under the power and type I error settings for both types of test)
#          Empirical type I error and predicted and empirical power for both underlying tests and all
#          proposed bias-correction methods, including our proposed randomization test.

# Example: Supplying the macro the below specified initial values:
#          lambda0=4
#          delta=1.5
#          sigma2=0.5
#          CV=0
#          m_bar=10
#          tau=1
#          nsim=100
#          alpha=0.050
#          beta=0.2
#          yields the following results:
#
#          ######################### 
#          ###Initial Parameters:### 
#          ######################### 
#          Baseline Hazard (lambda0):       4 
#          Treatment Effect (delta)  :      0 
#          Clustering Effect (sigma2):      0.5 
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
#          Parameter Estimate:              1.547017 
#          Empirical Variance:              0.3196896 
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
#          Parameter Estimate:              0.0322313 
#          Empirical Variance:              0.2290252 
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
#          Parameter Estimate:              1.571302 
#          Empirical Variance:              0.2264048 
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
delta<-1.5
sigma2<-0.5
pi<-1/2
CV<-0
m_bar<-10
tau<-1
nsim<-100
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
                         dblquad(five,0,tau,0,tau,dim=2) +
                         dblquad(six,0,tau,0,tau,dim=2)))
# Bread of VBV "V"
V=sqrt(m_bar)*(integral(h,0,tau))

# Sandwich Variance
sigma2_delta <- solve(V)%*%B%*%solve(V)
########################################################
#Sample Size - for given power, delta, and sigma2_delta#
########################################################
z<-((qnorm(1-beta)+qnorm(1-alpha/2))^2/delta^2)*(sigma2_delta)
N1<-ceiling(as.numeric(z))+mod(ceiling(as.numeric(z)),2)
powerz<-pnorm(sqrt(N1)*abs(delta)/sqrt(sigma2_delta)-qnorm(1-alpha/2),mean=0,sd=1)

n_init<-N1
while (n_init<((qt(1-beta,df=n_init-1)-qt(alpha/2,df=n_init-1))^2/delta^2)*(sigma2_delta))
{n_init<-n_init+1}
N2<-n_init+mod(n_init,2)
powert<-pt(sqrt(N2)*abs(delta)/sqrt(sigma2_delta)-qt(1-alpha/2,df=n_init-1),df=N2-1)

start<-proc.time()[[3]]
#############################
#Simulations part I - z-test#
#############################
delta_init<-delta
nsim_init <- nsim
n <- N1
cl <- 1:n

ndelta<-length(delta)
DELTA=matrix(NA,nsim,ndelta)
VARCZ=VARMD=VARKC=VARFG=VARMBN=matrix(NA,nsim,ndelta)
neg_haz_powerz=0
n_cl_powerz=rep(NA,nsim)
PVAL = rep(NA, nsim)

# Load permutation matrix
if(choose(n, n/2)<=50000){
  # Enumeration of all possible schemes
  emr <- t(combn(n, n / 2))
  R <- dim(emr)[1]
  pmt <- matrix(0, R, n)
  for (r in 1:R){
    pmt[r, emr[r,]] <- 1
  }
} else{
  S=50000
  # Simulation of permutations
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
  xi <- rnorm(n,mean=0,sd=sqrt(sigma2))
  trt <- rep(0, n)
  trt[sample(1:n,n/2,replace=F)] <- 1
  data <- matrix(nrow = 0, ncol = (3+ndelta))
  
  # Check if hazard is negative
  for (k in 1:n){
    if (trt[k]*delta+xi[k]+lambda0<1e-5) 
    {neg_haz_powerz=neg_haz_powerz+1
    skip<-1}
    else{
      if (CV==0) {temp<-m_bar}
      else {
        temp<-ceiling(rgamma(1,shape=CV^-2,rate=m_bar^-1*CV^-2))}
      nk[k,s] <- max(2,temp)
      Zk <- rep(trt[k],nk[k,s])
      Tk <- -log(runif(nk[k,s]))/c(Zk*delta+xi[k]+lambda0)
      Ck <- runif(nk[k,s])*tau
      data <- rbind(data,cbind(repmat(k,nk[k,s],1),Zk,pmin(Tk,Ck),Tk<=Ck))
    }
  }
  
  ID <- data[,1]
  Z <- data[,1+(1:ndelta),drop=FALSE]
  U <- data[,2+ndelta]
  Delta <- data[,ncol(data)]
  
  res=survMix(ID, Z, U, Delta)
  if (skip==0) {PVAL[s]=survPerm(ID, Z_cl=trt, U, Delta, pmt)} else {PVAL[s]=NA}
  DELTA[s,]=res$outest[,1]
  VARCZ[s,]=res$outest[,2]
  VARKC[s,]=res$outest[,3]
  VARMD[s,]=res$outest[,4]
  VARFG[s,]=res$outest[,5]
  VARMBN[s,]=res$outest[,6]
  n_cl_powerz[s]=res$outest[,7]
}

nsim<-nsim-sum(n_cl_powerz!=n)
DELTA<-DELTA[n_cl_powerz==n]
VARCZ_power<-VARCZ[n_cl_powerz==n]
VARKC_power<-VARKC[n_cl_powerz==n]
VARMD_power<-VARMD[n_cl_powerz==n]
VARFG_power<-VARFG[n_cl_powerz==n]
VARMBN_power<-VARMBN[n_cl_powerz==n]
################
#Power - z-test#
################
delta_est_powerz<-mean(DELTA)
sigma2_delta_est_powerz<-var(DELTA)
##Cai and Zeng
z_hat<-abs(DELTA-0)/sqrt(VARCZ_power)
powerz_CZ<-sum(z_hat>qnorm(1-alpha/2))/nsim

##Kauermann and Carroll
z_hat<-abs(DELTA-0)/sqrt(VARKC_power)
powerz_KC<-sum(z_hat>qnorm(1-alpha/2))/nsim

##Mancl and DeRouen
z_hat<-abs(DELTA-0)/sqrt(VARMD_power)
powerz_MD<-sum(z_hat>qnorm(1-alpha/2))/nsim

##Fay and Graubard
z_hat<-abs(DELTA-0)/sqrt(VARFG_power)
powerz_FG<-sum(z_hat>qnorm(1-alpha/2))/nsim

##Morel, Bokossa, and Neerchal
z_hat<-abs(DELTA-0)/sqrt(VARMBN_power)
powerz_MBN<-sum(z_hat>qnorm(1-alpha/2))/nsim

##Randomization test
power_permz<-mean(PVAL[n_cl_powerz==n]<0.05)
##############################
#Simulations part II - z-test#
##############################
delta<-rep(0,length(delta))
nsim <- nsim_init
cl <- 1:n

DELTA=matrix(NA,nsim,ndelta)
VARCZ=VARMD=VARKC=VARFG=VARMBN=matrix(NA,nsim,ndelta)
neg_haz_typeIz=0
n_cl_typeIz=rep(NA,nsim)
PVAL = rep(NA, nsim)

# Load permutation matrix
if(choose(n, n/2)<=50000){
  # Enumeration of all possible schemes
  emr <- t(combn(n, n / 2))
  R <- dim(emr)[1]
  pmt <- matrix(0, R, n)
  for (r in 1:R){
    pmt[r, emr[r,]] <- 1
  }
} else{
  S=50000
  # Simulation of permutations
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
  xi <- rnorm(n,mean=0,sd=sqrt(sigma2))
  trt <- rep(0, n)
  trt[sample(1:n,n/2,replace=F)] <- 1
  data <- matrix(nrow = 0, ncol = (3+ndelta))
  
  # Check if hazard is negative
  for (k in 1:n){
    if (trt[k]*delta+xi[k]+lambda0<1e-5) 
    {neg_haz_typeIz=neg_haz_typeIz+1
    skip<-1}
    else{
      if (CV==0) {temp<-m_bar}
      else {
        temp<-ceiling(rgamma(1,shape=CV^-2,rate=m_bar^-1*CV^-2))}
      nk[k,s] <- max(2,temp)
      Zk <- rep(trt[k],nk[k,s])
      Tk <- -log(runif(nk[k,s]))/c(Zk*delta+xi[k]+lambda0)
      Ck <- runif(nk[k,s])*tau
      data <- rbind(data,cbind(repmat(k,nk[k,s],1),Zk,pmin(Tk,Ck),Tk<=Ck))
    }
  }
  
  ID <- data[,1]
  Z <- data[,1+(1:ndelta),drop=FALSE]
  U <- data[,2+ndelta]
  Delta <- data[,ncol(data)]
  
  res=survMix(ID, Z, U, Delta)
  if (skip==0) {PVAL[s]=survPerm(ID, Z_cl=trt, U, Delta, pmt)} else {PVAL[s]=NA}
  DELTA[s,]=res$outest[,1]
  VARCZ[s,]=res$outest[,2]
  VARKC[s,]=res$outest[,3]
  VARMD[s,]=res$outest[,4]
  VARFG[s,]=res$outest[,5]
  VARMBN[s,]=res$outest[,6]
  n_cl_typeIz[s]=res$outest[,7]
}

nsim<-nsim-sum(n_cl_typeIz!=n)
DELTA<-DELTA[n_cl_typeIz==n]
VARCZ_typeI<-VARCZ[n_cl_typeIz==n]
VARKC_typeI<-VARKC[n_cl_typeIz==n]
VARMD_typeI<-VARMD[n_cl_typeIz==n]
VARFG_typeI<-VARFG[n_cl_typeIz==n]
VARMBN_typeI<-VARMBN[n_cl_typeIz==n]
#######################
#Type I-error - z-test#
#######################
delta_est_typeIz<-mean(DELTA)
sigma2_delta_est_typeIz<-var(DELTA)
##Cai and Zeng
z_hat<-abs(DELTA-0)/sqrt(VARCZ_typeI)
typeIz_CZ<-sum(z_hat>qnorm(1-alpha/2))/nsim

##Kauermann and Carroll
z_hat<-abs(DELTA-0)/sqrt(VARKC_typeI)
typeIz_KC<-sum(z_hat>qnorm(1-alpha/2))/nsim

##Mancl and DeRouen
z_hat<-abs(DELTA-0)/sqrt(VARMD_typeI)
typeIz_MD<-sum(z_hat>qnorm(1-alpha/2))/nsim

##Fay and Graubard
z_hat<-abs(DELTA-0)/sqrt(VARFG_typeI)
typeIz_FG<-sum(z_hat>qnorm(1-alpha/2))/nsim

##Morel, Bokossa, and Neerchal
z_hat<-abs(DELTA-0)/sqrt(VARMBN_typeI)
typeIz_MBN<-sum(z_hat>qnorm(1-alpha/2))/nsim

##Permutation Test
typeIz_perm<-mean(PVAL[n_cl_typeIz==n]<0.05)
#############################
#Simulations part I - t-test#
#############################
delta<-delta_init
nsim <- nsim_init
n <- N2
cl <- 1:n

DELTA=matrix(NA,nsim,ndelta)
VARCZ=VARMD=VARKC=VARFG=VARMBN=matrix(NA,nsim,ndelta)
neg_haz_powert=0
n_cl_powert=rep(NA,nsim)
PVAL = rep(NA, nsim)

# Load permutation matrix
if(choose(n, n/2)<=50000){
  # Enumeration of all possible schemes
  emr <- t(combn(n, n / 2))
  R <- dim(emr)[1]
  pmt <- matrix(0, R, n)
  for (r in 1:R){
    pmt[r, emr[r,]] <- 1
  }
} else{
  S=50000
  # Simulation of permutations
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
  xi <- rnorm(n,mean=0,sd=sqrt(sigma2))
  trt <- rep(0, n)
  trt[sample(1:n,n/2,replace=F)] <- 1
  data <- matrix(nrow = 0, ncol = (3+ndelta))
  
  # Check if hazard is negative
  for (k in 1:n){
    if (trt[k]*delta+xi[k]+lambda0<1e-5) 
    {neg_haz_powert=neg_haz_powert+1
    skip<-1}
    else{
      if (CV==0) {temp<-m_bar}
      else {
        temp<-ceiling(rgamma(1,shape=CV^-2,rate=m_bar^-1*CV^-2))}
      nk[k,s] <- max(2,temp)
      ind[k,s] <- as.numeric(2>temp)
      Zk <- rep(trt[k],nk[k,s])
      Tk <- -log(runif(nk[k,s]))/c(Zk*delta+xi[k]+lambda0)
      Ck <- runif(nk[k,s])*tau
      data <- rbind(data,cbind(repmat(k,nk[k,s],1),Zk,pmin(Tk,Ck),Tk<=Ck))
    }
  }
  
  ID <- data[,1]
  Z <- data[,1+(1:ndelta),drop=FALSE]
  U <- data[,2+ndelta]
  Delta <- data[,ncol(data)]
  
  res=survMix(ID, Z, U, Delta)
  if (skip==0) {PVAL[s]=survPerm(ID, Z_cl=trt, U, Delta, pmt)} else {PVAL[s]=NA}
  DELTA[s,]=res$outest[,1]
  VARCZ[s,]=res$outest[,2]
  VARKC[s,]=res$outest[,3]
  VARMD[s,]=res$outest[,4]
  VARFG[s,]=res$outest[,5]
  VARMBN[s,]=res$outest[,6]
  n_cl_powert[s]=res$outest[,7]
}

nsim<-nsim-sum(n_cl_powert!=n)
DELTA<-DELTA[n_cl_powert==n]
VARCZ_power<-VARCZ[n_cl_powert==n]
VARKC_power<-VARKC[n_cl_powert==n]
VARMD_power<-VARMD[n_cl_powert==n]
VARFG_power<-VARFG[n_cl_powert==n]
VARMBN_power<-VARMBN[n_cl_powert==n]
################
#Power - t-test#
################
delta_est_powert<-mean(DELTA)
sigma2_delta_est_powert<-var(DELTA)
##Cai and Zeng
t_hat<-abs(DELTA-0)/sqrt(VARCZ_power)
powert_CZ<-sum(t_hat>qt(1-alpha/2,n-1))/nsim

##Kauermann and Carroll
t_hat<-abs(DELTA-0)/sqrt(VARKC_power)
powert_KC<-sum(t_hat>qt(1-alpha/2,n-1))/nsim

##Mancl and DeRouen
t_hat<-abs(DELTA-0)/sqrt(VARMD_power)
powert_MD<-sum(t_hat>qt(1-alpha/2,n-1))/nsim

##Fay and Graubard
t_hat<-abs(DELTA-0)/sqrt(VARFG_power)
powert_FG<-sum(t_hat>qt(1-alpha/2,n-1))/nsim

##Morel, Bokossa, and Neerchal
t_hat<-abs(DELTA-0)/sqrt(VARMBN_power)
powert_MBN<-sum(t_hat>qt(1-alpha/2,n-1))/nsim

##Permutation Test
powert_perm<-mean(PVAL[n_cl_powert==n]<0.05)
##############################
#Simulations part II - t-test#
##############################
delta<-rep(0,length(delta))
nsim <- nsim_init
cl <- 1:n

DELTA=matrix(NA,nsim,ndelta)
VARCZ=VARMD=VARKC=VARFG=VARMBN=matrix(NA,nsim,ndelta)
neg_haz_typeIt=0
n_cl_typeIt=rep(NA,nsim)
PVAL = rep(NA, nsim)

# Load permutation matrix
if(choose(n, n/2)<=50000){
  # Enumeration of all possible schemes
  emr <- t(combn(n, n / 2))
  R <- dim(emr)[1]
  pmt <- matrix(0, R, n)
  for (r in 1:R){
    pmt[r, emr[r,]] <- 1
  }
} else{
  S=50000
  # Simulation of permutations
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
  xi <- rnorm(n,mean=0,sd=sqrt(sigma2))
  trt <- rep(0, n)
  trt[sample(1:n,n/2,replace=F)] <- 1
  data <- matrix(nrow = 0, ncol = (3+ndelta))
  
  # Check if hazard is negative
  for (k in 1:n){
    if (trt[k]*delta+xi[k]+lambda0<1e-5) 
    {neg_haz_typeIt=neg_haz_typeIt+1
    skip<-1}
    else{
      if (CV==0) {temp<-m_bar}
      else {
        temp<-ceiling(rgamma(1,shape=CV^-2,rate=m_bar^-1*CV^-2))}
      nk[k,s] <- max(2,temp)
      ind[k,s] <- as.numeric(2>temp)
      Zk <- rep(trt[k],nk[k,s])
      Tk <- -log(runif(nk[k,s]))/c(Zk*delta+xi[k]+lambda0)
      Ck <- runif(nk[k,s])*tau
      data <- rbind(data,cbind(repmat(k,nk[k,s],1),Zk,pmin(Tk,Ck),Tk<=Ck))
    }
  }
  
  ID <- data[,1]
  Z <- data[,1+(1:ndelta),drop=FALSE]
  U <- data[,2+ndelta]
  Delta <- data[,ncol(data)]
  
  res=survMix(ID, Z, U, Delta)
  if (skip==0) {PVAL[s]=survPerm(ID, Z_cl=trt, U, Delta, pmt)} else {PVAL[s]=0}
  DELTA[s,]=res$outest[,1]
  VARCZ[s,]=res$outest[,2]
  VARKC[s,]=res$outest[,3]
  VARMD[s,]=res$outest[,4]
  VARFG[s,]=res$outest[,5]
  VARMBN[s,]=res$outest[,6]
  n_cl_typeIt[s]=res$outest[,7]
}

nsim<-nsim-sum(n_cl_typeIt!=n)
DELTA<-DELTA[n_cl_typeIt==n]
VARCZ_typeI<-VARCZ[n_cl_typeIt==n]
VARKC_typeI<-VARKC[n_cl_typeIt==n]
VARMD_typeI<-VARMD[n_cl_typeIt==n]
VARFG_typeI<-VARFG[n_cl_typeIt==n]
VARMBN_typeI<-VARMBN[n_cl_typeIt==n]
######################
#TypeI-error - t-test#
######################
delta_est_typeIt<-mean(DELTA)
sigma2_delta_est_typeIt<-var(DELTA)
##Cai and Zeng
t_hat<-abs(DELTA-0)/sqrt(VARCZ_typeI)
typeIt_CZ<-sum(t_hat>qt(1-alpha/2,n-1))/nsim

##Kauermann and Carroll
t_hat<-abs(DELTA-0)/sqrt(VARKC_typeI)
typeIt_KC<-sum(t_hat>qt(1-alpha/2,n-1))/nsim

##Mancl and DeRouen
t_hat<-abs(DELTA-0)/sqrt(VARMD_typeI)
typeIt_MD<-sum(t_hat>qt(1-alpha/2,n-1))/nsim

##Fay and Graubard
t_hat<-abs(DELTA-0)/sqrt(VARFG_typeI)
typeIt_FG<-sum(t_hat>qt(1-alpha/2,n-1))/nsim

##Morel, Bokossa, and Neerchal
t_hat<-abs(DELTA-0)/sqrt(VARMBN_typeI)
typeIt_MBN<-sum(t_hat>qt(1-alpha/2,n-1))/nsim

##Permutation Test
typeIt_perm<-mean(PVAL[n_cl_typeIt==n]<0.05)


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
  "Treatment Effect (delta)  :     ",delta_init,"\n",
  "Clustering Effect (sigma2):     ",sigma2,"\n",
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
  "Negative Hazard Occurance:      ",neg_haz_powerz,"\n",
  "Number of Simulations Excluded: ",sum(n_cl_powerz!=N1),"\n",
  "Parameter Estimate:             ",delta_est_powerz,"\n",
  "Empirical Variance:             ",sigma2_delta_est_powerz,"\n",
  "\n",
  "Predicted:                      ",format(round(powerz, digits = 3), nsmall=3),"\n",
  "Cai and Zeng:                   ",format(round(powerz_CZ, digits = 3), nsmall=3),"\n",
  "Kauermann and Carroll:          ",format(round(powerz_KC, digits = 3), nsmall=3),"\n",
  "Mancl and DeRouen:              ",format(round(powerz_MD, digits = 3), nsmall=3),"\n",
  "Morel, Bokossa, and Neerchal:   ",format(round(powerz_MBN, digits = 3), nsmall=3),"\n",
  "#########################","\n",
  "###z-test type-I error###","\n",
  "#########################","\n",
  "Negative Hazard Occurance:      ",neg_haz_typeIz,"\n",
  "Number of Simulations Excluded: ",sum(n_cl_typeIz!=N1),"\n",
  "Parameter Estimate:             ",delta_est_typeIz,"\n",
  "Empirical Variance:             ",sigma2_delta_est_typeIz,"\n",
  "\n",
  "Cai and Zeng:                   ",format(round(typeIz_CZ, digits = 3), nsmall=3),"\n",
  "Kauermann and Carroll:          ",format(round(typeIz_KC, digits = 3), nsmall=3),"\n",
  "Mancl and DeRouen:              ",format(round(typeIz_MD, digits = 3), nsmall=3),"\n",
  "Morel, Bokossa, and Neerchal:   ",format(round(typeIz_MBN, digits = 3), nsmall=3),"\n",
  "##################","\n",
  "###t-test Power###","\n",
  "##################","\n",
  "Negative Hazard Occurance:      ",neg_haz_powert,"\n",
  "Number of Simulations Excluded: ",sum(n_cl_powert!=N2),"\n",
  "Parameter Estimate:             ",delta_est_powert,"\n",
  "Empirical Variance:             ",sigma2_delta_est_powert,"\n",
  "\n",
  "Predicted:                      ",format(round(powert, digits = 3), nsmall=3),"\n",
  "Cai and Zeng:                   ",format(round(powert_CZ, digits = 3), nsmall=3),"\n",
  "Kauermann and Carroll:          ",format(round(powert_KC, digits = 3), nsmall=3),"\n",
  "Mancl and DeRouen:              ",format(round(powert_MD, digits = 3), nsmall=3),"\n",
  "Morel, Bokossa, and Neerchal:   ",format(round(powert_MBN, digits = 3), nsmall=3),"\n",
  "Permutation:                    ",format(round(powert_perm, digits = 3), nsmall=3),"\n",
  "#########################","\n",
  "###t-test type-I error###","\n",
  "#########################","\n",
  "Negative Hazard Occurance:      ",neg_haz_typeIt,"\n",
  "Number of Simulations Excluded: ",sum(n_cl_typeIt!=N2),"\n",
  "Parameter Estimate:             ",delta_est_typeIt,"\n",
  "Empirical Variance:             ",sigma2_delta_est_typeIt,"\n",
  "\n",
  "Cai and Zeng:                   ",format(round(typeIt_CZ, digits = 3), nsmall=3),"\n",
  "Kauermann and Carroll:          ",format(round(typeIt_KC, digits = 3), nsmall=3),"\n",
  "Mancl and DeRouen:              ",format(round(typeIt_MD, digits = 3), nsmall=3),"\n",
  "Morel, Bokossa, and Neerchal:   ",format(round(typeIt_MBN, digits = 3), nsmall=3),"\n",
  "Permutation:                    ",format(round(typeIt_perm, digits = 3), nsmall=3),"\n",
  "\n",
  "Simulation Time: ",Elapsed,"(seconds)")
