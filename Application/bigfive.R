#### Title: Regularized LS-LV
#### Author: Tra Le

#########################################################################################
######                         Empirical Data Application                           #####
#########################################################################################

#### Big Five Data
#install.packages("qgraph")
library(qgraph)
data("big5")
DATA <- scale(big5)
source('CCLSLV')
source('LSLVLASSO.R')

nzero <- 0:239
PS <- nzero/240

### CCLSLV
CARDvec <- 240*5-PS*240*5
CARDvec2 <- CARDvec/5
IS2 <- matrix(ncol = 3, nrow = length(CARDvec2)) # index of spareseness for model selection
for (j in 1:length(CARDvec2)){
  a <- IS_CCLSLV(DATA, R = 5, card = rep(CARDvec2[j],5), MaxIter = 20, eps = 10^-6, nstarts = 5)
  IS2[j,] <- c(CARDvec2[j],a$value, a$vaf) 
}

CCLSLV_model <- MULTISTART_CCLSLV(DATA, R = 5, CARD = rep(CARDvec2[which.max(IS2[,2])], 5), MaxIter = 200, eps = 10^-6, nstarts = 10)
CCLSLV_loadings <- CCLSLV_model$loadings

### LSLV-LASSO (need to try several lambda values to get the same level of sparseness with CCLSLV)
LSLVLASSO_model <- MULTISTART_LSLVLASSO(DATA, R = 5, MaxIter = 200, eps = 10^-6, nstarts = 10, lambda = 208)
LSLVLASSO_loadings <- LSLVLASSO_model$loadings
LSLVLASSO_scores <- LSLVLASSO_model$scores
V_a <- sum((LSLVLASSO_scores %*% t(LSLVLASSO_loadings))^2)
V_oo <- sum(DATA^2)
PEV <- V_a/V_oo #formula from Zhengguo
# double check with another formula
PEV <- 1-(LSLVLASSO_model$Loss/V_oo)  

### SGCCA (try several sparsity values)
blocks <- list(bl1 = DATA)
fit_sgcca <- rgcca(blocks = blocks, method = "sgcca", sparsity = .37, superblock = FALSE, ncomp = 5,
                   scheme = "factorial",
                   comp_orth = TRUE,
                   verbose = TRUE)
nonzero1 <- sum(fit_sgcca$a$bl1[,1]!=0)
nonzero2 <- sum(fit_sgcca$a$bl1[,2]!=0)
nonzero3 <- sum(fit_sgcca$a$bl1[,3]!=0)
nonzero4 <- sum(fit_sgcca$a$bl1[,4]!=0)
nonzero5 <- sum(fit_sgcca$a$bl1[,5]!=0)

