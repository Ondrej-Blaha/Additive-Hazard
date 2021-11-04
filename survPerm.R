###########################################################################################################################
# Macro: survPerm

# Author: Ondrej Blaha, Fan Li

# Creation Date: November 3, 2021 (R version 4.0.5)

# Purpose: This macro is used to calculate the permutation p-value for the randomization test which is
#          inspired by the Braun and Feng 2001 JASA paper.
#
#          This macro is used within the main code a delivers one p-value per one out of the nsim_init simulations
#          while using the data generated within the main code as input.          

# Required Parameters: 
#          ID: Identification vector
#          X_cl: Treatment allocation (on cluster level)
#          Y: Vector of outcomes (survival time)
#          Delta: Censoring indicator vector
#          pmt: Permutation matrix of treatment allocation (capped at 50,000 permutations)

# Output:
#          Permutation p-value (scalar)

# Example: Supplying the macro the following data: ID_ex=rep(1:6,3); trt_ex=c(rep(0,3),rep(1,3)); set.seed(6747583)
#          Y_ex=runif(length(ID_ex),min=0,max=3); Delta_ex=rbinom(length(ID_ex),1,0.3)
#          pmt_ex=matrix(0, 20, 6); for (r in 1:20){pmt_ex[r, t(combn(6,3))[r,]]<-1}
#          yields the following results: survPerm(ID=ID_ex,X_cl=trt_ex,Y=Y_ex,Delta=Delta_ex,pmt=pmt_ex)
#          [1,] 0.95
#
#          Disclaimer: the presented example is a minimal working example and the data generation process and
#          the data itself supplied to the macro do not represent the actual data used within the main program.

###########################################################################################################################
survPerm=function(ID, X_cl, Y, Delta, pmt){
  
  # Sort observations by time
  b <- order(Y)
  Y <- sort(Y)
  ID <- ID[b]
  Delta <- Delta[b]
  
  # Increment of time
  dY <- Y-c(0,Y[1:(length(Y)-1)])
  # Each row is an individual, each column is a specific time point
  IndYY <- (t(repmat(Y,length(Y),1))>=repmat(Y,length(Y),1))
  UID <- sort(unique(ID))
  IDind <- zeros(length(UID), length(Y))
  for (i in 1:length(UID)){
    IDind[i, ID==UID[i]] <- 1
  }
  
  # General set up
  n <- length(UID)
  ny <- length(Y)
  # The rate of counting process of event, or dN(t)
  NN <- diag(Delta)
  # Sum of Y(t) at each t
  temp0 <- t(IndYY)%*%repmat(1,ny,1)
  
  # Estimation of marginal baseline hazard
  dHY <- colSums(IndYY*NN)/c(temp0)
  HY <- cumsum(dHY)
  
  # Permutation distribution of score
  ######################################################################################
  dH <- dHY
  epsilon <- NN-repmat(t(dH),ny,1)
  score <- as.numeric(IDind%*%(as.matrix(IndYY*epsilon)%*%repmat(1,ny,1)))
  ######################################################################################
  
  obs = sum(X_cl*score)
  Pdist = c(pmt%*%score)
  
  # Permutation p-value
  return(mean(abs(Pdist) >= abs(obs)))
}