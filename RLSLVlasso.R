#### Title: Regularized LS-LV method
#### Author: Tra Le 
#### Supervisor: Dr. Katrijn Van Deun
#### Created: February 1, 2023
#### Last modified: May 29, 2023

#########################################################################################
#####                  RLSLV with LASSO for correlated components                   #####
#########################################################################################

##1. RLSLV with LASSO
RLSLV <- function(DATA, R, lambda, MaxIter, eps){
  convAO <- 0
  iter <- 1
  Lossc <- 1
  Lossvec <- Lossc   #Used to compute convergence criterium
  svd1 <- svd(DATA, R, R)
  P1 <- matrix(rnorm(J*R), ncol = R, nrow = J)
  P2 <- svd1$v %*% diag(svd1$d[1:R])/sqrt(I)
  loadings <- P1*.3 + P2*.7 #initialize P 
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

##2. LOSS for LASSO
LOSS <- function(DATA, SCORES, LOADINGS, LAMBDA){
  XHAT <- SCORES%*%t(LOADINGS)
  res <- sum(rowSums((XHAT-DATA)^2))
  penalty <- sum(abs(LOADINGS))
  loss <- res+LAMBDA*penalty
  return(loss)
}

