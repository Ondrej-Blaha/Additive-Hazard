###########################################################################################################################
# Functions: one, two, three, four, five, six, h

# Author: Ondrej Blaha, Fan Li

# Creation Date: December 13, 2021 (R version 4.0.5)

# Purpose: This set of functions is used as a functional description of the seven quantities denoted in our paper as
#          g1, g2, g3, g4, g5, g6, h. These functions reflect our choice treating the frailty parameter xi
#          as normally distributed.
#
#          These functions are used in the main code in order to calculate variance within the sample size calculation.
#          Technically, these functions return a single value f(t) for a selected time t and f(t,s) for selected
#          times t and s respectively. while sourcing information from the main code about the initial parameters
#          (lambda0,delta,sigma2,tau). However, these functions are used in the main code to integrate over using
#          pracma package functions 'integral' and 'dblquad'.

# Required Parameters: 
#          t,s: Time (scalars)

# Required Global Parameters:
#          lambda0: Fixed baseline hazard (scalar)
#          delta:   Pre-specified treatment effect (scalar)
#          sigma2:  Pre-specified clustering effect (scalar)
#          tau:     Duration of the study (scalar)

# Output:
#          g1(t), g2(t,s), g3(t,s), g4(t,s), g5(t,s), g6(t,s), h(t): values of the functions g1,...,g6 and h for 
#          specific time points t and s. 

# Example: Using global parameters: tau=1; lambda0=4; delta=1.5; sigma2=0.5
#          and supplying the functions 'one' and 'h' with the time value of t=0.5 
#          yields results of g1(t), h(t): one(t); h(t)
#          [1] 0.05510155
#          [1] 0.01155467
#          Supplying functions 'two', 'three', 'four', 'five', and 'six' with t=0.5; s=0.25
#          yields results of g2(t,s),...,g6(t,s): two(t,s); three(t,s); four(t,s); five(t,s); six(t,s)
#          [1] 0.05673516
#          [1] -0.03881071
#          [1] -0.07645863
#          [1] 0.02040359
#          [1] 0.03962246

###########################################################################################################################
one<-function(t){
  # Censoring distribution
  G_t <- (1-t/tau)
  
  # Z_bar
  mu_t <- 1/(1+exp(delta*t))
  
  # Density functions
  f0_t <- exp(-lambda0*t+sigma2*t^2/2)*(lambda0-sigma2*t)
  f1_t <- exp(-lambda0*t+sigma2*t^2/2-delta*t)*(lambda0-sigma2*t+delta)

  # Final quantity
  one <- (1-pi)*G_t*mu_t^2*f0_t + pi*G_t*(1-mu_t)^2*f1_t
  
  return(one)
}

two<-function(t,s){
  # Censoring distribution
  G_ts <- (1-t/tau)*(1-s/tau)
  
  # Z_bar
  mu_t <- 1/(1+exp(delta*t))
  mu_s <- 1/(1+exp(delta*s))
  
  # Density functions
  f0_ts <- exp(-lambda0*(t+s)+sigma2*(t+s)^2/2)*(sigma2+(-lambda0+sigma2*(t+s))^2)
  f1_ts <- exp(-lambda0*(t+s)+sigma2*(t+s)^2/2-delta*(t+s))*(sigma2+(-lambda0+sigma2*(t+s)-delta)^2)

  # Final Quantity
  two <- (1-pi)*G_ts*mu_t*mu_s*f0_ts + pi*G_ts*(1-mu_t)*(1-mu_s)*f1_ts
  
  return(two)
}

three<-function(t,s){
  # Censoring distribution
  G_ts <- (1-t/tau)*(1-s/tau)
  
  # Z_bar
  mu_t <- 1/(1+exp(delta*t))
  mu_s <- 1/(1+exp(delta*s))
  
  # Survival functions
  dt_S0_ts <- exp(-lambda0*(t+s)+sigma2*(t+s)^2/2)*(-lambda0+sigma2*(t+s))
  ds_S0_ts <- exp(-lambda0*(t+s)+sigma2*(t+s)^2/2)*(-lambda0+sigma2*(t+s))

  # Final quantity
  three <- (1-pi)*G_ts*mu_t*mu_s*(dt_S0_ts*(lambda0-s*sigma2) + ds_S0_ts*(lambda0-t*sigma2))
  
  return(three)
}

four<-function(t,s){
  # Censoring distribution
  G_ts <- (1-t/tau)*(1-s/tau)
  
  # Z_bar
  mu_t <- 1/(1+exp(delta*t))
  mu_s <- 1/(1+exp(delta*s))
  
  # Survival functions
  dt_S1_ts <- exp(-lambda0*(t+s)+sigma2*(t+s)^2/2-delta*(t+s))*(-lambda0+sigma2*(t+s)-delta)
  ds_S1_ts <- exp(-lambda0*(t+s)+sigma2*(t+s)^2/2-delta*(t+s))*(-lambda0+sigma2*(t+s)-delta)
  
  # Final quantity
  four <- pi*G_ts*(1-mu_t)*(1-mu_s)*(dt_S1_ts*(lambda0-s*sigma2+delta) + ds_S1_ts*(lambda0-t*sigma2+delta))
  
  return(four)
}

five<-function(t,s){
  # Censoring distribution
  G_ts <- (1-t/tau)*(1-s/tau)
  
  # Z_bar
  mu_t <- 1/(1+exp(delta*t))
  mu_s <- 1/(1+exp(delta*s))
  
  # Survival function
  S0_ts <- exp(-lambda0*(t+s)+sigma2*(t+s)^2/2)

  # Final quantity
  five <- (1-pi)*G_ts*mu_t*mu_s*S0_ts*(lambda0-s*sigma2)*(lambda0-t*sigma2)
  
  return(five)
}

six<-function(t,s){
  # Censoring distribution
  G_ts <- (1-t/tau)*(1-s/tau)
  
  # Z_bar
  mu_t <- 1/(1+exp(delta*t))
  mu_s <- 1/(1+exp(delta*s))
  
  # Survival function
  S1_ts <- exp(-lambda0*(t+s)+sigma2*(t+s)^2/2-delta*(t+s))

  # Final quantity
  six <- pi*G_ts*(1-mu_t)*(1-mu_s)*S1_ts*(lambda0-s*sigma2+delta)*(lambda0-t*sigma2+delta)
  
  return(six)
}

h<-function(t){
  # Censoring distribution
  G_t <- (1-t/tau)

  # Survival functions
  S0_t <- exp(-lambda0*t+sigma2*t^2/2)
  S1_t <- exp(-lambda0*t+sigma2*t^2/2-delta*t)
  
  # Final quantity
  h <- G_t/(1/((1-pi)*S0_t) + 1/(pi*S1_t))
  
  return(h)
}