###########################################################################################################################
# Procedure: VIF

# Author: Fan Li, Ondrej Blaha

# Creation Date: November 3, 2021 (R version 4.0.5)

# Purpose: This procedure is the source of values displayed in Figure 2 of the paper.
#          It sources global (initial) parameters: lambda0, beta0, theta0, CV, m_bar, and tau described below
#          and returns values of ICC, VIF, and estimated sample size (number of clusters)
#          as well as estimated sandwich variance of the treatment effect.
#
#          This procedure is not used within the main code and serves exclusively as a source for
#          data example in Figure 2.

# Dependencies: pracma, Rfast - libraries
#               Integrals.R - File containing supporting functions

# Required Global Parameters:
#          lambda0: Initial baseline hazard
#          beta0: Initial treatment effect
#          theta0: Initial clustering effect
#          CV: Desired coefficient of variation
#          m_bar: Desired average cluster size
#          tau: Duration of the study (maximum follow-up time)

# Output:
#          ICC: Estimated intraclass correlation
#          varbeta: Estimated sandwich variance for treatment effect
#          VIF: Variance inflation factor
#          nhat: Estimated sample size (number of clusters)

# Example: Supplying the procedure with the following global parameters (inspired by the STRIDE trial):
#          lambda0=0.2; beta0=-0.1; theta0=0.03; CV=0.5; m_bar=64; tau=3
#          yields the following results: VIF()
#          $rho
#          [1] 0.3376434
#          $VIF
#          [1] 27.67383
#          $varbeta
#          [1] 0.1622737
#          $nhat
#          [1] 130

###########################################################################################################################
library(pracma)
library(Rfast)
source("Integrals.R")

VIF<-function(){
  # ICC
  rho = (dblquad(two,0,tau,0,tau,dim=2) + 
           dblquad(three,0,tau,0,tau,dim=2) +
           dblquad(four,0,tau,0,tau,dim=2) +
           dblquad(five,0,tau,0,tau,dim=2))/integral(one,0,tau)
  
  # Sandwich variance
  B=integral(one,0,tau) + 
    ((1+CV^2)*m_bar-1)*(dblquad(two,0,tau,0,tau,dim=2) + 
                         dblquad(three,0,tau,0,tau,dim=2) +
                         dblquad(four,0,tau,0,tau,dim=2) +
                         dblquad(five,0,tau,0,tau,dim=2))
  A=m_bar*(integral(six,0,tau))^2
  varbeta <- B/A
  
  # VIF 
  VIF=(1+((1+CV^2)*m_bar-1)*rho)
  
  # Sample Size
  # Get initial value using z-test
  nhat=ceiling((qnorm(1-0.05/2)+qnorm(1-0.2))^2*(varbeta)/(beta0)^2)
  # Iterative computation for t-test
  while(nhat < (qt(1-0.05/2,df=nhat-1)+qt(1-0.2,df=nhat-1))^2*(varbeta)/(beta0)^2){nhat <- nhat+1}
  
  return(list(rho=rho,VIF=VIF,varbeta=varbeta,nhat=nhat))
}