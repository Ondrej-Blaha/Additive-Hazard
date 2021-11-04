###########################################################################################################################
# Macro: survMix (contains functions: GEEN, RESULTS)

# Author: Ondrej Blaha, Fan Li

# Creation Date: November 3, 2021 (R version 4.0.5)

# Purpose: This macro is used to calculate parameter estimates and bias-corrected and uncorrected
#          sandwich variance estimates. It takes survival data on entry in form of 
#          (ID vector, matrix of covariates, vector of survival times, censoring indicator vector) and
#          outputs a vector of estimates (treatment effect, CZ uncorrected sandwich variance,
#          KC, MD, FG, MBN types of corrected sandwich variances, number of clusters).
#
#          This macro is connected to the main code which supplies survMix macro with generated data
#          and uses output values for each simulation to process and summarize results for the entire
#          set of nsim_init simulations.

# Required Parameters: 
#          ID: Identification vector
#          X: Design matrix (model covariates)
#          Y: Vector of outcomes (survival time)
#          Delta: Censoring indicator vector

# Output:
#          beta: treatment effect estimate
#          robust: uncorrected sandwich variance estimate
#          varKC: Kauermann and Carroll type bias-corrected sandwich variance estimate
#          varMD: Mancl and DeRouen type bias-corrected sandwich variance estimate
#          varFG: Fay and Graubard type bias-corrected sandwich variance estimate
#          varMBN: Morel, Bokossa, and Neerchal type bias-corrected sandwich variance estimate
#          n_cl: number of clusters

# Example: Supplying the macro the following data: ID_ex=rep(1:6,3); X_ex=rep(c(rep(0,3),rep(1,3)),3); set.seed(67937163)
#          Y_ex=runif(length(ID_ex),min=0,max=3); Delta_ex=rbinom(length(ID_ex),1,0.3)
#          yields the following results: survMix(ID=ID_ex,X=X_ex,Y=Y_ex,Delta=Delta_ex)
#                 Estimate  BC0 (CZ)   BC1 (KC)   BC2 (MD)   BC3 (FG)  BC4 (MBN) # Clusters
#          [1,] -0.1061628 0.0242168 0.03481818 0.05254111 0.03481818 0.06239724          6
#
#          Disclaimer: the presented example is a minimal working example and the data generation process and
#          the data itself supplied to the macro do not represent the actual data used within the main program.

###########################################################################################################################
survMix=function(ID, X, Y, Delta){
  
  GEEN <- function(ID, X, Y, Delta, dY, IndYY){
    # General set up
    UID <- sort(unique(ID))
    n <- length(UID)
    ny <- length(Y)
    nbeta <- dim(as.matrix(X))[2]
    # the rate of counting process of event, or dN(t)
    NN <- diag(Delta)
    
    # Sum of Y(t) at each t
    temp0 <- t(IndYY)%*%repmat(1,ny,1)
    # bar{X}(t) at each t (row) and for each dimension (column)
    barX <- (t(IndYY)%*%X)/repmat(temp0,1,nbeta)
    nom <- zeros(nbeta,1)
    denom <- zeros(nbeta,nbeta)
    
    # Least-squares estimation of beta
    for (k in 1:nbeta){
      tempk <- repmat(as.matrix(as.matrix(X)[,k]),1,ny)-repmat(barX[,k],ny,1)
      nom[k] <- sum(IndYY*tempk*NN)
      for (s in 1:nbeta){
        temps <- repmat(as.matrix(as.matrix(X)[,s]),1,ny)-repmat(barX[,s],ny,1)
        denom[k,s] <- sum((IndYY*tempk*temps)%*%dY)
      }
    }
    beta <- solve(denom)%*%nom
    
    # Estimation of marginal and conditional baseline hazard
    Xbeta <- X%*%beta
    dHY <- colSums(IndYY*(NN-(Xbeta)%*%t(dY)))/c(temp0)
    HY <- cumsum(dHY)
    
    # Obtain martingale increment
    dH <- dHY
    epsilon <- NN-repmat(t(dH),ny,1)-(X%*%beta)%*%t(dY)
    nom2 <- zeros(n,nbeta)
    denom2 <- zeros(nbeta,nbeta)
    
    # Influence function and variance of beta allowing for four types of bias-correction
    ########################################################################################
    Omega <- Omega_star <- array(0,c(nbeta,nbeta,n))
    for (k in 1:nbeta){
      tempk <- repmat(as.matrix(as.matrix(X)[,k]),1,ny)-repmat(barX[,k],ny,1)
      nom2[,k] <- IDind%*%(as.matrix(IndYY*tempk*epsilon)%*%repmat(1,ny,1))
      for (s in 1:nbeta){
        temps <- repmat(as.matrix(as.matrix(X)[,s]),1,ny)-repmat(barX[,s],ny,1)
        temps_star <- repmat(as.matrix(as.matrix(X)[,s]),1,ny)
        Omega[k,s,] <- IDind%*%((IndYY*tempk*temps)%*%dY)
        Omega_star[k,s,] <- IDind%*%((IndYY*tempk*temps_star)%*%dY)
        denom2[k,s] <- sum(Omega[k,s,])
      }
    }
    ########################################################################################
    # Model-based variance
    inustar <- 1/denom2
    inustartr <- inustar # Transpose is itself given scalar
    # Sandwich variance
    U_i_array <- U_i_bc1_array <- U_i_bc2_array <- 
      U_i_bc3_array <- array(0,c(nbeta,1,n))
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
    deltan=min(0.5, nbeta/(n-nbeta))
    psin=max(1, ((ny-1)*n)/((ny-nbeta)*(n-1))*c(UUtran*inustar)/nbeta)
    varMBN=((ny-1)*n)/((ny-nbeta)*(n-1))*robust + c(deltan*psin*inustar)
    
    return(list(beta=beta,Xbeta=Xbeta,dH=dH,HY=HY,naive=inustar,
                robust=robust,varKC=varKC,varMD=varMD,varFG=varFG,varMBN=varMBN,
                U_i_array=U_i_array,U_i_bc1_array=U_i_bc1_array,
                U_i_bc2_array=U_i_bc2_array,U_i_bc3_array=U_i_bc3_array,
                UUtran=UUtran,UUbc1=UUbc1,UUbc2=UUbc2,UUbc3=UUbc3,n_cl=n))
  }
  
  RESULTS=function(beta, robust, varKC, varMD, varFG, varMBN, n_cl){
    BC0=(robust)
    BC1=(varKC)
    BC2=(varMD)
    BC3=(varFG)
    BC4=(varMBN)
    outbeta=cbind(beta,BC0,BC1,BC2,BC3,BC4,n_cl)
    colnames(outbeta)<-c("Estimate","BC0 (CZ)","BC1 (KC)","BC2 (MD)","BC3 (FG)","BC4 (MBN)","# Clusters")
    return(list(outbeta=outbeta))
  }
  
  # Sort observations by time
  b <- order(Y)
  Y <- sort(Y)
  X <- X[b]
  ID <- ID[b]
  Delta <- Delta[b]
  
  # Increment of time
  dY <- Y-c(0,Y[1:(length(Y)-1)])
  IndYY <- (t(repmat(Y,length(Y),1))>=repmat(Y,length(Y),1))
  UID <- sort(unique(ID))
  IDind <- zeros(length(UID), length(Y))
  for (i in 1:length(UID)){
    IDind[i, ID==UID[i]] <- 1
  }
  
  # Estimate beta
  BETA_RES <- GEEN(ID, X, Y, Delta, dY, IndYY)
  beta <- BETA_RES$beta
  Xbeta <- BETA_RES$Xbeta
  naive <- BETA_RES$naive
  dH <- BETA_RES$dH
  HY <- BETA_RES$HY
  robust <- BETA_RES$robust
  varKC <- BETA_RES$varKC
  varMD <- BETA_RES$varMD
  varFG <- BETA_RES$varFG
  varMBN <- BETA_RES$varMBN
  n_cl <- BETA_RES$n_cl
  
  # Final Results
  RESULTS(beta, robust, varKC, varMD, varFG, varMBN, n_cl)
}