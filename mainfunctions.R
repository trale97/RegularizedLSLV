#### Title: Regularized LS-LV method
#### Author: Tra Le 
#### Supervisor: Dr. Katrijn Van Deun
#### Created: February 1, 2023
#### Last modified: April 7, 2023

#########################################################################################
#####                  Preliminaries (only the measurement part)                    #####
#########################################################################################

##1. USLPCA function 
USLPCA <- function(DATA, R, CARD, MaxIter, eps) {
  J <- dim(DATA)[2]
  I <- dim(DATA)[1]
  convAO <- 0
  iter <- 1
  Lossc <- 1
  Lossvec <- Lossc        #Used to compute convergence criterium
  loadings <- INITLOADINGS(DATA, R, CARD)
  A <- loadings
  scores <- matrix(nrow = I, ncol = R)
  if (sum(length(CARD)) > 1){
    CARDc <- J-CARD #Number of loadings -per PC- that have to be zero
    while (convAO == 0) {
      #1. Update component scores
      XP <- DATA %*% loadings/ (sqrt(I))
      scores <- sqrt(I) * svd(XP, R, R)$u %*% t(svd(XP, R, R)$v)
      A <- t(DATA) %*% scores/ I
      #Calculate loss
      Loss <- RESIDUAL(DATA,scores,loadings)
      Lossvec <- c(Lossvec,Loss)
      #2. Update loadings
      loadings <- A
      for (r in 1:R){
        ind <- sort(abs(A[,r]), index.return = TRUE)
        loadings[ind$ix[1:CARDc[r]],r] <- 0
      } 
      #Calculate loss
      Lossu <- RESIDUAL(DATA,scores,loadings)
      Lossvec <- c(Lossvec,Lossu)
      if (iter > MaxIter) {
        convAO <- 1
      }
      if (abs(Lossc-Lossu) < eps) {
        convAO <- 1
      }
      iter <- iter + 1
      Lossc <- Lossu
    }} else {
      CARDc <- J*R-CARD #Number of loadings -over all PCs- that have to be zero
      while (convAO == 0) {
        #1. Update component scores
        XP <- DATA %*% loadings / (sqrt(I))
        scores <- sqrt(I) * svd(XP, R, R)$u %*% t(svd(XP, R, R)$v)
        A <- t(DATA) %*% scores / I
        #Calculate loss
        Loss <- RESIDUAL(DATA,scores,loadings)
        Lossvec <- c(Lossvec,Loss)
        #2. Update loadings
        loadings <- A
        ind <- sort(abs(A), index.return = TRUE)
        loadings[ind$ix[1:CARDc]] <- 0
        #Calculate loss
        Lossu <- RESIDUAL(DATA,scores,loadings)
        Lossvec <- c(Lossvec,Lossu)
        if (iter > MaxIter ) {
          convAO <- 1
        }
        if (abs(Lossc-Lossu) < eps) {
          convAO <- 1
        }
        iter <- iter + 1
        Lossc <- Lossu
      }
    }
  uslpca <- list('scores' = scores, 'loadings' = loadings, 'Lossvec' = Lossvec, 'Residual' = Loss)
  return(uslpca)
}

##2. Regularized LSLV with lasso 
RLSLV <- function(DATA, R, lambda, MaxIter, eps){
  convAO <- 0
  iter <- 1
  Lossc <- 1
  Lossvec <- Lossc   #Used to compute convergence criterium
  svd1 <- svd(DATA, R, R)
  loadings <- matrix(rnorm(J*R), nrow = J, ncol = R) #initialize P 
  scores <- svd1$u 
  diffT <- 0
  diffP <- 0
  while (convAO == 0) {
    iter0 <- 1
    Losst <- 1
    Lossvec0 <- Losst
    convT0 <- 0
    Lossvec1 <- 1
    #1. Update component scores 
    while(convT0 == 0){
      Lossu1old <- LOSS(DATA,scores,loadings,lambda)
      E <- DATA - scores%*%t(loadings) 
      for (r in 1:R){
        Er <- E + scores[,r]%*%t(loadings[,r])
        num <- Er%*%loadings[,r]
        scores[,r] <- sqrt(I)*num/sqrt(sum(num^2))
        Lossu1 <- LOSS(DATA,scores,loadings,lambda)
        diffT <- c(diffT,Lossu1old-Lossu1)
        Lossu1old <- Lossu1
        Lossvec1 <- c(Lossvec1, Lossu1)
      }
      #t(scores)%*%scores
      #Calculate loss
      Lossu0 <- LOSS(DATA,scores,loadings,lambda)
      Lossvec0 <- c(Lossvec0,Lossu0)
      # check convergence
      if (iter0 > MaxIter) {
        convT0 <- 1
      }
      if (abs(Losst-Lossu0) < eps){
        convT0 <- 1
      }
      iter0 <- iter0 + 1
      Losst <- Lossu0
    }
    Loss <- LOSS(DATA, scores, loadings, lambda)
    
    #2. Update loadings
    Lossu1old <- LOSS(DATA,scores,loadings,lambda)
    E <- DATA - scores%*%t(loadings) 
    for (r in 1:R){
      Er <- E+scores[,r]%*%t(loadings[,r])
      crosstEr <- t(Er)%*%scores[,r]
      loadings[,r]<-sign(crosstEr)*apply(cbind(abs(crosstEr)-lambda/2,0),1,max)/I
    }
  
    #Calculate loss
    Lossu <- LOSS(DATA,scores,loadings,lambda)
    Lossvec <- c(Lossvec,Lossu)
    if (iter > MaxIter) {
      convAO <- 1
    }
    if (abs(Lossc-Lossu) < eps) {
      convAO <- 1
    }
    iter <- iter + 1
    Lossc <- Lossu
  }
  return_rlslv <- list()
  return_rlslv$loadings <- loadings
  return_rlslv$scores <- scores
  return_rlslv$loss <- Loss
  return_rlslv$lossvec <- Lossvec
  
  return(return_rlslv)
}

###############################################################################
############################SUBFUNCTIONS#######################################
###############################################################################

##1. LOSS
RESIDUAL <- function(DATA, SCORES, LOADINGS){
  XHAT <- SCORES%*%t(LOADINGS)
  res <- sum(rowSums((XHAT-DATA)^2))
  
  return(res)
}

##2. INTIAL LOADINGS
INITLOADINGS <- function(DATA, R, CARD){
  J <- dim(DATA)[2]
  I <- dim(DATA)[1]
  P <- matrix(rnorm(J*R), ncol = R, nrow = J)
  if (length(CARD)>1){
    CARDc <- J-CARD
    for (r in 1:R){
      ind <- sort(abs(P[,r]), index.return = TRUE)
      P[ind$ix[1:CARDc[r]],r] <- 0
    } 
  } else {
    CARDc <- J*R-CARD
    ind <- sort(abs(P), index.return = TRUE)
    P[ind$ix[1:CARDc]] <- 0
  }
  loadings <- P
  return(loadings)
}


##3. LOSS for lasso
LOSS <- function(DATA, SCORES, LOADINGS, LAMBDA){
  XHAT <- SCORES%*%t(LOADINGS)
  res <- sum(rowSums((XHAT-DATA)^2))
  penalty <- sum(abs(LOADINGS))
  loss <- res+LAMBDA*penalty
  return(loss)
}

#########################################################################################
#####                             IS for MODEL SELECTION                            #####
#########################################################################################
IS_USLPCA <- function(X, R, INIT, card, MaxIter, eps){
  #1. Obtain USLPCA solution
  USLPCA <- USLPCA(X, R, INIT, CARD = card, MaxIter, eps)
  #2. Calculate BIC
  ssrespca <- 0
  J <- dim(X)[2]
  va <- 0
  nrzero <- 0
  vzero <- 0
  #calculation SSRES PCA
  vzero <- sum(rowSums(X^2))
  d <- svd(X)$d
  ssrespca <- ssrespca+sum((d[(R+1):length(d)])^2)
  loading <- c(round(USLPCA$loadings,3))
  va <- va+sum((d[1:R])^2)
  nrzero <- nrzero+sum(loading==0)
  dfjspca <- 0
  for (jr in 1:(J*R)){
    dfjspca <- dfjspca+length(unique(loading[jr][loading[jr]!=0]))
  }
  bicvalue <- USLPCA$Residual/ssrespca+log(dfjspca)
  vs <- vzero-USLPCA$Residual
  nrcoef <- J*R
  nrzeqcoef <- nrcoef-dfjspca
  IS <- list()
  IS$value <- va*vs/vzero^2*nrzeqcoef/nrcoef
  IS$vaf <- 1-(USLPCA$Residual/vzero)
  IS$propunique <- dfjspca/(nrcoef-nrzero)
  IS$propzero <- nrzero/nrcoef
  #return(bicvalue)
  return(IS)
}

#########################################################################################
#####                             Multistart procedure                              #####
#########################################################################################
MULTISTART <- function(DATA, R, CARD, MaxIter, eps, nstarts, lambda){
  if(missing(nstarts)){
    nstarts <- 20
  } 
  
  if(missing(method)){
    method <- "correlated"
  }
  
  if(R == 1){
    stop("Parameter R = 1 is not allowed.")
  }
  
  Pout3d <- list()
  Tout3d <- list()
  LOSS <- array()
  LOSSvec <- list()
  
  for (n in 1:nstarts){
    if(method == "orthogonal"){
      result <- USLPCA(DATA, R, CARD, MaxIter, eps)
    } else if (method == "correlated"){
      result <- RLSLV(DATA, R, lambda, MaxIter, eps)
    }
    Pout3d[[n]] <- result$loadings
    Tout3d[[n]] <- result$scores
    LOSS[n] <- result$Residual
    LOSSvec[[n]] <- result$Lossvec
  }
  # choose solution with lowest loss value
  k <- which(LOSS == min(LOSS))
  if (length(k)>1){
    pos <- sample(1:length(k), 1)
    k <- k[pos]
  }
  
  return_varselect <- list()
  return_varselect$loadings <- Pout3d[[k]]
  return_varselect$scores <- Tout3d[[k]]
  return_varselect$Lossvec <- LOSSvec
  
  return(return_varselect)
}
