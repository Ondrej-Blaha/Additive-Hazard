###########################################################################################################################
# Functions: one, two, three, four, five, six

# Author: Ondrej Blaha, Fan Li

# Creation Date: November 3, 2021 (R version 4.0.5)

# Purpose: This set of functions is used as a functional description of the six quantities denoted in our paper as
#          g1, g2, g3, g4, g5, g6. These functions reflect our choice treating the frailty parameter xi
#          as normally distributed.
#
#          These functions are used in the main code in order to calculate variance within the sample size calculation.
#          Technically, these functions return a single value f(t) for a selected time t and f(t,s) for selected
#          times t and s respectively. while sourcing information from the main code about the initial parameters
#          (lambda0,beta0,theta0,tau). However, these functions are used in the main code to integrate over using
#          pracma package functions 'integral' and 'dblquad'.

# Required Parameters: 
#          t,s: Time (scalars)

# Required Global Parameters:
#          lambda0: Initial baseline hazard
#          beta0: Initial treatment effect
#          theta0: Initial value of the frailty parameter
#          tau: Duration of the study

# Output:
#          f(t), f(t,s): values of the functions g1,...,g6 for specific time points t and s. 

# Example: Using global parameters: tau=1; lambda0=4; beta0=1.5; theta0=0.5
#          and supplying the functions 'one' and 'six' with the time value of t=0.5 
#          yields results of f(t): one(t); six(t)
#          [1] 0.05510155
#          [1] 0.01155467
#          Supplying functions 'two', 'three', 'four', and 'five' with t=0.5; s=0.25
#          yields results of f(t,s): two(t,s); three(t,s); four(t,s); five(t,s)
#          [1] 0.05673516
#          [1] -0.05840255
#          [1] -0.05686679
#          [1] 0.06002606

###########################################################################################################################
one<-function(t){
  # Censoring distribution
  Q_t <- 1-t/tau
  
  # X_bar
  mu_t <- 1/(1+exp(beta0*t))
  
  # Survival functions
  f_ij_1 <- exp(-lambda0*t-beta0*t+t^2*theta0/2)*(lambda0+beta0-t*theta0)
  f_ij_0 <- exp(-lambda0*t+t^2*theta0/2)*(lambda0-t*theta0)
  
  # Final Quantity
  one <- 1/2*Q_t*(1-mu_t)^2*f_ij_1 + 1/2*Q_t*(0-mu_t)^2*f_ij_0
  
  return(one)
}

two<-function(t,s){
  # Censoring distribution
  Q_ts <- (1-t/tau)*(1-s/tau)
  
  # X_bar
  mu_t <- 1/(1+exp(beta0*t))
  mu_s <- 1/(1+exp(beta0*s))
  
  # Survival functions
  f_ij_ik_1 <- exp(-lambda0*(t+s)-beta0*(t+s)+(t+s)^2*theta0/2)*(theta0+(-lambda0-beta0+(t+s)*theta0)^2)
  f_ij_ik_0 <- exp(-lambda0*(t+s)+(t+s)^2*theta0/2)*(theta0+(-lambda0+(t+s)*theta0)^2)
  
  # Final Quantity
  two <- 1/2*Q_ts*(1-mu_t)*(1-mu_s)*f_ij_ik_1 + 1/2*Q_ts*(0-mu_t)*(0-mu_s)*f_ij_ik_0
  
  return(two)
}

three<-function(t,s){
  # Censoring distribution
  Q_ts <- (1-t/tau)*(1-s/tau)
  
  # X_bar
  mu_t <- 1/(1+exp(beta0*t))
  mu_s <- 1/(1+exp(beta0*s))
  
  # Survival functions
  dt_F_ij_ik_1 <- exp(-lambda0*(t+s)-beta0*(t+s)+(t+s)^2*theta0/2)*(-lambda0-beta0+(t+s)*theta0)
  dt_F_ij_ik_0 <- exp(-lambda0*(t+s)+(t+s)^2*theta0/2)*(-lambda0+(t+s)*theta0)
  
  # Final Quantity
  three <- 1/2*Q_ts*(1-mu_t)*(1-mu_s)*(dt_F_ij_ik_1)*(lambda0-s*theta0+beta0) + 1/2*Q_ts*(0-mu_t)*(0-mu_s)*(dt_F_ij_ik_0)*(lambda0-s*theta0)
  
  return(three)
}

four<-function(t,s){
  # Censoring distribution
  Q_ts <- (1-t/tau)*(1-s/tau)
  
  # X_bar
  mu_t <- 1/(1+exp(beta0*t))
  mu_s <- 1/(1+exp(beta0*s))
  
  # Survival functions
  ds_F_ij_ik_1 <- exp(-lambda0*(t+s)-beta0*(t+s)+(t+s)^2*theta0/2)*(-lambda0-beta0+(t+s)*theta0)
  ds_F_ij_ik_0 <- exp(-lambda0*(t+s)+(t+s)^2*theta0/2)*(-lambda0+(t+s)*theta0)
  
  # Final Quantity
  four <- 1/2*Q_ts*(1-mu_t)*(1-mu_s)*(ds_F_ij_ik_1)*(lambda0-t*theta0+beta0) + 1/2*Q_ts*(0-mu_t)*(0-mu_s)*(ds_F_ij_ik_0)*(lambda0-t*theta0)
  
  return(four)
}

five<-function(t,s){
  # Censoring distribution
  Q_ts <- (1-t/tau)*(1-s/tau)
  
  # X_bar
  mu_t <- 1/(1+exp(beta0*t))
  mu_s <- 1/(1+exp(beta0*s))
  
  # Survival functions
  F_ij_ik_1 <- exp(-lambda0*(t+s)-beta0*(t+s)+(t+s)^2*theta0/2)
  F_ij_ik_0 <- exp(-lambda0*(t+s)+(t+s)^2*theta0/2)
  
  # Final Quantity
  five <- 1/2*Q_ts*(1-mu_t)*(1-mu_s)*F_ij_ik_1*(lambda0-t*theta0+beta0)*(lambda0-s*theta0+beta0) + 1/2*Q_ts*(0-mu_t)*(0-mu_s)*F_ij_ik_0*(lambda0-t*theta0)*(lambda0-s*theta0)
  
  return(five)
}

six<-function(t){
  # Censoring distribution
  Q_t <- 1-t/tau
  
  # X_bar
  mu_t <- 1/(1+exp(beta0*t))
  
  # Survival functions
  S_1 <- exp(-lambda0*t-beta0*t+t^2*theta0/2)
  S_0 <- exp(-lambda0*t+t^2*theta0/2)
  
  # Final Quantity
  six <- 1/2*Q_t*(0-mu_t)^2*S_0 + 1/2*Q_t*(1-mu_t)^2*S_1
  
  return(six)
}