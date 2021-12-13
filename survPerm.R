###########################################################################################################################
# Macro: survPerm

# Author: Ondrej Blaha, Fan Li

# Creation Date: December 13, 2021 (R version 4.0.5)

# Purpose: This macro is used to calculate the permutation p-value for the randomization test which is
#          inspired by the Braun and Feng 2001 JASA paper.
#
#          This macro is used within the main code a delivers one p-value per one out of the nsim simulations
#          while using the data generated within the main code as input.          

# Required Parameters: 
#          ID: Identification vector
#          Z_cl: Treatment allocation (on cluster level)
#          U: Vector of outcomes (observed survival time)
#          Delta: Censoring indicator vector
#          pmt: Permutation matrix of treatment allocation (capped at 50,000 permutations)

# Output:
#          Permutation p-value (scalar)

# Example: Supplying the macro with the following data: ID_ex=rep(1:6,3); trt_ex=c(rep(0,3),rep(1,3)); set.seed(6747583)
#          U_ex=runif(length(ID_ex),min=0,max=3); Delta_ex=rbinom(length(ID_ex),1,0.3)
#          pmt_ex=matrix(0, 20, 6); for (r in 1:20){pmt_ex[r, t(combn(6,3))[r,]]<-1}
#          yields the following results: survPerm(ID=ID_ex,Z_cl=trt_ex,U=U_ex,Delta=Delta_ex,pmt=pmt_ex)
#          [1] 0.95
#
#          Disclaimer: the presented example is a minimal working example and the data generation process and
#          the data itself supplied to the macro do not represent the actual data used within the main program.

###########################################################################################################################
survPerm=function(ID, Z_cl, U, Delta, pmt){
  
  # Sort observations by time
  b <- order(U)
  U <- sort(U)
  ID <- ID[b]
  Delta <- Delta[b]
  
  # Increment of time
  dU <- U-c(0,U[1:(length(U)-1)])
  # Each row is an individual, each column is a specific time point
  IndUU <- (t(repmat(U,length(U),1))>=repmat(U,length(U),1))
  UID <- sort(unique(ID))
  IDind <- zeros(length(UID), length(U))
  for (i in 1:length(UID)){
    IDind[i, ID==UID[i]] <- 1
  }
  
  # General set up
  n <- length(UID)
  nu <- length(U)
  # The rate of counting process of event, or dN(t)
  NN <- diag(Delta)
  # Sum of U(t) at each t
  temp0 <- t(IndUU)%*%repmat(1,nu,1)
  
  # Estimation of marginal baseline hazard
  dHU <- colSums(IndUU*NN)/c(temp0)
  HU <- cumsum(dHU)
  
  # Permutation distribution of score
  ######################################################################################
  dH <- dHU
  epsilon <- NN-repmat(t(dH),nu,1)
  score <- as.numeric(IDind%*%(as.matrix(IndUU*epsilon)%*%repmat(1,nu,1)))
  ######################################################################################
  
  obs = sum(Z_cl*score)
  Pdist = c(pmt%*%score)
  
  # Permutation p-value
  return(mean(abs(Pdist) >= abs(obs)))
}