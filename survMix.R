###########################################################################################################################
# Macro: survMix (contains functions: GEEN, RESULTS)

# Author: Ondrej Blaha, Fan Li

# Creation Date: December 13, 2021 (R version 4.0.5)

# Purpose: This macro is used to calculate parameter estimates and bias-corrected and uncorrected
#          sandwich variance estimates. It takes observed survival times on entry in form of 
#          (ID vector, vector of treatment allocation, vector of survival times, censoring indicator vector) and
#          outputs a vector of estimates (treatment effect, CZ uncorrected sandwich variance,
#          KC, MD, FG, MBN types of corrected sandwich variances, number of clusters).
#
#          This macro is connected to the main code which supplies survMix macro with generated data
#          and uses output values for each simulation to process and summarize results for the entire
#          set of nsim simulations.

# Required Parameters: 
#          ID: Identification vector
#          Z: Design matrix (vector of treatment allocation)
#          U: Vector of outcomes (observed survival time)
#          Delta: Censoring indicator vector

# Output:
#          delta_est: treatment effect estimate
#          robust: uncorrected sandwich variance estimate
#          varKC: Kauermann and Carroll type bias-corrected sandwich variance estimate
#          varMD: Mancl and DeRouen type bias-corrected sandwich variance estimate
#          varFG: Fay and Graubard type bias-corrected sandwich variance estimate
#          varMBN: Morel, Bokossa, and Neerchal type bias-corrected sandwich variance estimate
#          n_cl: number of clusters

# Example: Supplying the macro the following data: ID_ex=rep(1:6,3); Z_ex=rep(c(rep(0,3),rep(1,3)),3); set.seed(67937163)
#          U_ex=runif(length(ID_ex),min=0,max=3); Delta_ex=rbinom(length(ID_ex),1,0.3)
#          yields the following results: survMix(ID=ID_ex,Z=Z_ex,U=U_ex,Delta=Delta_ex)
#                 Estimate  BC0 (CZ)   BC1 (KC)   BC2 (MD)   BC3 (FG)  BC4 (MBN) # Clusters
#          [1,] -0.1061628 0.0242168 0.03481818 0.05254111 0.03481818 0.06239724          6
#
#          Disclaimer: the presented example is a minimal working example and the data generation process and
#          the data itself supplied to the macro do not represent the actual data used within the main program.

###########################################################################################################################
survMix=function(ID, Z, U, Delta){
  
  GEEN <- function(ID, Z, U, Delta, dU, IndUU){
    # General set up
    UID <- sort(unique(ID))
    n <- length(UID)
    nu <- length(U)
    ndelta <- dim(as.matrix(Z))[2]
    # the rate of counting process of event, or dN(t)
    NN <- diag(Delta)
    
    # Sum of U(t) at each t
    temp0 <- t(IndUU)%*%repmat(1,nu,1)
    # bar{Z}(t) at each t (row) and for each dimension (column)
    barZ <- (t(IndUU)%*%Z)/repmat(temp0,1,ndelta)
    nom <- zeros(ndelta,1)
    denom <- zeros(ndelta,ndelta)
    
    # Least-squares estimation of delta
    for (k in 1:ndelta){
      tempk <- repmat(as.matrix(as.matrix(Z)[,k]),1,nu)-repmat(barZ[,k],nu,1)
      nom[k] <- sum(IndUU*tempk*NN)
      for (s in 1:ndelta){
        temps <- repmat(as.matrix(as.matrix(Z)[,s]),1,nu)-repmat(barZ[,s],nu,1)
        denom[k,s] <- sum((IndUU*tempk*temps)%*%dU)
      }
    }
    delta_est <- solve(denom)%*%nom
    
    # Estimation of marginal and conditional baseline hazard
    Zdelta_est <- Z%*%delta_est
    dHU <- colSums(IndUU*(NN-(Zdelta_est)%*%t(dU)))/c(temp0)
    HU <- cumsum(dHU)
    
    # Obtain martingale increment
    dH <- dHU
    epsilon <- NN-repmat(t(dH),nu,1)-(Z%*%delta_est)%*%t(dU)
    nom2 <- zeros(n,ndelta)
    denom2 <- zeros(ndelta,ndelta)
    
    # Influence function and variance of delta est. allowing for four types of bias-corrections
    ###########################################################################################
    Omega <- Omega_star <- array(0,c(ndelta,ndelta,n))
    for (k in 1:ndelta){
      tempk <- repmat(as.matrix(as.matrix(Z)[,k]),1,nu)-repmat(barZ[,k],nu,1)
      nom2[,k] <- IDind%*%(as.matrix(IndUU*tempk*epsilon)%*%repmat(1,nu,1))
      for (s in 1:ndelta){
        temps <- repmat(as.matrix(as.matrix(Z)[,s]),1,nu)-repmat(barZ[,s],nu,1)
        temps_star <- repmat(as.matrix(as.matrix(Z)[,s]),1,nu)
        Omega[k,s,] <- IDind%*%((IndUU*tempk*temps)%*%dU)
        Omega_star[k,s,] <- IDind%*%((IndUU*tempk*temps_star)%*%dU)
        denom2[k,s] <- sum(Omega[k,s,])
      }
    }
    ########################################################################################
    # Model-based variance
    inustar <- 1/denom2
    inustartr <- inustar # Transpose is itself given scalar
    # Sandwich variance
    U_i_array <- U_i_bc1_array <- U_i_bc2_array <- 
      U_i_bc3_array <- array(0,c(ndelta,1,n))
    UUtran <- UUbc1 <- UUbc2 <- UUbc3 <- 0
    for(i in 1:n){ 
      ######################################################################################
      Bcm <- 1/c(1-Omega_star[,,i]*inustar)
      Hi <- 1/sqrt(1-min(0.75,c(Omega_star[,,i]*inustar)))
      ######################################################################################
      U_i <- nom2[i,]
      U_c <- Bcm*U_i
      U_d <- Hi*U_i
      UUtran_c <- (U_i)^2
      UUtran <- UUtran+UUtran_c
      UUbc1_c <- (U_c)*(U_i)
      UUbc1 <- UUbc1+UUbc1_c
      UUbc2_c <- (U_c)^2
      UUbc2 <- UUbc2+UUbc2_c
      UUbc3_c <- (U_d)^2
      UUbc3 <- UUbc3+UUbc3_c
      U_i_array[,,i] <- U_i
      U_i_bc1_array[,,i] <- sqrt(Bcm)*U_i
      U_i_bc2_array[,,i] <- U_c
      U_i_bc3_array[,,i] <- U_d
    }
    
    # BC0 or usual Sandwich estimator of Cai and Zeng (2011);     
    robust=c(inustar*UUtran*inustartr)
    
    # BC1 or Variance estimator that extends Kauermann and Carroll (2001);
    varKC=c(inustar*(UUbc1+(UUbc1))*inustartr/2)
    
    # BC2 or Variance estimator that extends Mancl and DeRouen (2001);
    varMD=c(inustar*UUbc2*inustartr)
    
    # BC3 or Variance estimator that extends Fay and Graubard (2001);
    varFG=c(inustar*UUbc3*inustartr)
    
    # BC4 or Variance estimator that extends Morel, Bokossa, and Neerchal (2003);
    deltan=min(0.5, ndelta/(n-ndelta))
    psin=max(1, ((nu-1)*n)/((nu-ndelta)*(n-1))*c(UUtran*inustar)/ndelta)
    varMBN=((nu-1)*n)/((nu-ndelta)*(n-1))*robust + c(deltan*psin*inustar)
    
    return(list(delta_est=delta_est,Zdelta_est=Zdelta_est,dH=dH,HU=HU,naive=inustar,
                robust=robust,varKC=varKC,varMD=varMD,varFG=varFG,varMBN=varMBN,
                U_i_array=U_i_array,U_i_bc1_array=U_i_bc1_array,
                U_i_bc2_array=U_i_bc2_array,U_i_bc3_array=U_i_bc3_array,
                UUtran=UUtran,UUbc1=UUbc1,UUbc2=UUbc2,UUbc3=UUbc3,n_cl=n))
  }
  
  RESULTS=function(delta_est, robust, varKC, varMD, varFG, varMBN, n_cl){
    BC0=(robust)
    BC1=(varKC)
    BC2=(varMD)
    BC3=(varFG)
    BC4=(varMBN)
    outest=cbind(delta_est,BC0,BC1,BC2,BC3,BC4,n_cl)
    colnames(outest)<-c("Estimate","BC0 (CZ)","BC1 (KC)","BC2 (MD)","BC3 (FG)","BC4 (MBN)","# Clusters")
    return(list(outest=outest))
  }
  
  # Sort observations by time
  b <- order(U)
  U <- sort(U)
  Z <- Z[b]
  ID <- ID[b]
  Delta <- Delta[b]
  
  # Increment of time
  dU <- U-c(0,U[1:(length(U)-1)])
  IndUU <- (t(repmat(U,length(U),1))>=repmat(U,length(U),1))
  UID <- sort(unique(ID))
  IDind <- zeros(length(UID), length(U))
  for (i in 1:length(UID)){
    IDind[i, ID==UID[i]] <- 1
  }
  
  # Estimate delta
  BETA_RES <- GEEN(ID, Z, U, Delta, dU, IndUU)
  delta_est <- BETA_RES$delta_est
  Zdelta_est <- BETA_RES$Zdelta_est
  naive <- BETA_RES$naive
  dH <- BETA_RES$dH
  HU <- BETA_RES$HU
  robust <- BETA_RES$robust
  varKC <- BETA_RES$varKC
  varMD <- BETA_RES$varMD
  varFG <- BETA_RES$varFG
  varMBN <- BETA_RES$varMBN
  n_cl <- BETA_RES$n_cl
  
  # Final Results
  RESULTS(delta_est, robust, varKC, varMD, varFG, varMBN, n_cl)
}