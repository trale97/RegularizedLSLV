### Model selection for USLPCA using Index of Sparseness

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