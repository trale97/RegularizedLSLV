#### Title: Regularized LS-LV
#### Author: Tra Le
#### Created: April 10, 2023
#### Last modified: April 11, 2023

#########################################################################################
######                         USLPCA simulation analysis                           #####
#########################################################################################

# Clear workspace
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
rm(list = ls())

# load library
library(doParallel)

# load simulation setup
load("Info_simulation.RData")
Info_matrix = Infor_simulation$design_matrix_replication
Ndatasets = Infor_simulation$n_data_sets
MatrixInfo = Infor_simulation$design_matrix_replication

betas <- c(-0.02, 0.15, .5)

# register cores
no_cores <- detectCores()
c1 <- makePSOCKcluster(floor( no_cores*.4))
registerDoParallel(c1)

# run simulation
Simulation = foreach(i = 1:Ndatasets,
                     .packages = c("lm.beta"),
                     .combine = rbind)%dopar%{
                       source('USLPCA.R', local = TRUE)
                       
                       # Load data
                       load(paste0("CORRELATED/composite",i,".RData")) # composite-based
                       Out1 = list(data = out$data, W = out$W, Ptrue = out$Ptrue, 
                                   n_items = out$n_items, s_size = out$s_size)
                       
                       
                       load(paste0("CORRELATED/factor",i,".RData")) #factor-based
                       Out2 = list(data = out$data, Y = out$Y, Ptrue = out$Ptrue, 
                                   n_items = out$n_items, s_size = out$s_size)
                       
                       I <- Out2$s_size 
                       K <- Out2$n_items
                       data1 <- scale(Out1$data)
                       data2 <- scale(Out2$data)
                       ## Measurement model
                       # population values
                       Ptrue1 <- Out1$Ptrue
                       Ptrue2 <- Out2$Ptrue
                       
                       result1 <- MULTISTART_USLPCA(DATA = data1, R = 3, CARD = rep(K,3), MaxIter = 200, eps = 10^-6, nstarts = 10)
                       result2 <- MULTISTART_USLPCA(DATA = data2, R = 3, CARD = rep(K,3), MaxIter = 200, eps = 10^-6, nstarts = 10)
                       
                       Pmatrix1 <- result1$loadings
                       Pmatrix2 <- result2$loadings
                       
                       # zero/non-zero structure recovery rate + Tucker's congruence
                       perm <- gtools::permutations(3, 3)
                       corrate1 <- vector(length = nrow(perm))
                       corrate2 <- vector(length = nrow(perm))
                       spcr_tucongrT1 <- vector(length = nrow(perm))
                       spcr_tucongrT2 <- vector(length = nrow(perm))
                       for (p in 1:nrow(perm)) {
                         corrate1[p] <-num_correct(Ptrue1, round(Pmatrix1[,perm[p,]],3))
                         corrate2[p] <-num_correct(Ptrue2, round(Pmatrix2[,perm[p,]],3))
                         spcr_tucongrT1[p] <- sum(diag(abs(psych::factor.congruence(Ptrue1,Pmatrix1[,perm[p,]]))))/3
                         spcr_tucongrT2[p] <- sum(diag(abs(psych::factor.congruence(Ptrue2,Pmatrix2[,perm[p,]]))))/3
                       }
                       corrate_composite <- max(corrate1)
                       corrate_factor <- max(corrate2)
                       tucongr_composite <- max(spcr_tucongrT1)
                       tucongr_factor <- max(spcr_tucongrT2)
                       
                       index1 <- which.max(corrate1)
                       index2 <- which.max(corrate2)
                       
                       Pmatrix1 <- Pmatrix1[,perm[index1,]]
                       Pmatrix2 <- Pmatrix2[,perm[index2,]]
                       corsign1 <- sign(diag(cor(Ptrue1, Pmatrix1)))
                       corsign2 <- sign(diag(cor(Ptrue2, Pmatrix2)))
                       
                       Pmatrix1 <- Pmatrix1%*%diag(corsign1)
                       Pmatrix2 <- Pmatrix2%*%diag(corsign2)
                       
                       ## Structural model
                       Tmatrix1 <- result1$scores
                       Tmatrix2 <- result2$scores
                       Tmatrix1 <- Tmatrix1[,perm[index1,]]%*%diag(corsign1) # account for permutation + reflection
                       Tmatrix2 <- Tmatrix2[,perm[index2,]]%*%diag(corsign2)
                       
                       #beta1 <- lm.beta(lm(data1[,(K*3+1)] ~ Tmatrix1))
                       #beta2 <- lm.beta(lm(data2[,(K*3+1)] ~ Tmatrix2))
                       
                       ## measure summary
                       measure_summary <- matrix(c('Factor', Out2$s_size, Out2$n_items, corrate_factor, tucongr_factor,
                                                   'Composite', Out1$s_size, Out1$n_items, corrate_composite, tucongr_composite), nrow = 2, byrow = T)
                       colnames(measure_summary) <- c("Model", "I", "K", "Zero_Rec", "Tucker")
                       # Saving performance measures
                       output = list(Pmatrix_factor = Pmatrix2, Tmatrix_factor = Tmatrix2,
                                     Pmatrix_composite = Pmatrix1, Tmatrix_composite = Tmatrix1,
                                     measures = measure_summary, Pfactor = Ptrue2, Pcomposite = Ptrue1)
                       
                       save(output, file = paste0("CORRELATED/USLPCA", i, ".RData"))
                       
                     }
# stop Cluster
stopCluster(c1)

########################### Model selection
# Clear workspace
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
rm(list = ls())

# load library
library(doParallel)

# load simulation setup
load("Info_simulation.RData")
Info_matrix = Infor_simulation$design_matrix_replication
Ndatasets = Infor_simulation$n_data_sets
MatrixInfo = Infor_simulation$design_matrix_replication

betas <- c(-0.02, 0.15, .5)

# register cores
no_cores <- detectCores()
c1 <- makePSOCKcluster(floor( no_cores*.4))
registerDoParallel(c1)

# run simulation
Simulation = foreach(i = 1:Ndatasets,
                     .packages = c("lm.beta"),
                     .combine = rbind)%dopar%{
                       source('USLPCA.R', local = TRUE)
                       source('IS_USLPCA.R', local = TRUE)
                       
                       # Load data
                       load(paste0("CORRELATED/composite",i,".RData")) # composite-based
                       Out1 = list(data = out$data, W = out$W, Ptrue = out$Ptrue, 
                                   n_items = out$n_items, s_size = out$s_size)
                       
                       
                       load(paste0("CORRELATED/factor",i,".RData")) #factor-based
                       Out2 = list(data = out$data, Y = out$Y, Ptrue = out$Ptrue, 
                                   n_items = out$n_items, s_size = out$s_size)
                       
                       I <- Out2$s_size 
                       K <- Out2$n_items
                       data1 <- scale(Out1$data)
                       data2 <- scale(Out2$data)
                       ## Measurement model
                       # population values
                       Ptrue1 <- Out1$Ptrue
                       Ptrue2 <- Out2$Ptrue
                       
                       cardmat <- matrix(rep(c(3, 5, 7, 9), 3),nrow = 4)
                       avec <- matrix(nrow = 4, ncol = 2)
                       bvec <- matrix(nrow = 4, ncol = 2)
                       for (j in 1:4){
                         a <- IS_USLPCA(X = data1, R = 3, card = cardmat[j,], MaxIter = 20, eps = 10^-6, nstarts = 5)
                         b <- IS_USLPCA(X = data2, R = 3, card = cardmat[j,], MaxIter = 20, eps = 10^-6, nstarts = 5)
                         avec[j,] <- c(cardmat[j,1],a$value)
                         bvec[j,] <- c(cardmat[j,1],b$value)
                       }
                       selmodelindex1 <- which.max(avec[,2]) 
                       selmodelindex2 <- which.max(bvec[,2])
                       
                       result1 <- MULTISTART_USLPCA(DATA = data1, R = 3, CARD = rep(avec[selmodelindex1,1],3), MaxIter = 200, eps = 10^-6, nstarts = 10)
                       result2 <- MULTISTART_USLPCA(DATA = data2, R = 3, CARD = rep(avec[selmodelindex2,1],3), MaxIter = 200, eps = 10^-6, nstarts = 10)
                       
                       Pmatrix1 <- result1$loadings
                       Pmatrix2 <- result2$loadings
                       
                       # zero/non-zero structure recovery rate + Tucker's congruence
                       perm <- gtools::permutations(3, 3)
                       corrate1 <- vector(length = nrow(perm))
                       corrate2 <- vector(length = nrow(perm))
                       spcr_tucongrT1 <- vector(length = nrow(perm))
                       spcr_tucongrT2 <- vector(length = nrow(perm))
                       for (p in 1:nrow(perm)) {
                         corrate1[p] <-num_correct(Ptrue1, round(Pmatrix1[,perm[p,]],3))
                         corrate2[p] <-num_correct(Ptrue2, round(Pmatrix2[,perm[p,]],3))
                         spcr_tucongrT1[p] <- sum(diag(abs(psych::factor.congruence(Ptrue1,Pmatrix1[,perm[p,]]))))/3
                         spcr_tucongrT2[p] <- sum(diag(abs(psych::factor.congruence(Ptrue2,Pmatrix2[,perm[p,]]))))/3
                       }
                       corrate_composite <- max(corrate1)
                       corrate_factor <- max(corrate2)
                       tucongr_composite <- max(spcr_tucongrT1)
                       tucongr_factor <- max(spcr_tucongrT2)
                       
                       index1 <- which.max(corrate1)
                       index2 <- which.max(corrate2)
                       
                       Pmatrix1 <- Pmatrix1[,perm[index1,]]
                       Pmatrix2 <- Pmatrix2[,perm[index2,]]
                       corsign1 <- sign(diag(cor(Ptrue1, Pmatrix1)))
                       corsign2 <- sign(diag(cor(Ptrue2, Pmatrix2)))
                       
                       Pmatrix1 <- Pmatrix1%*%diag(corsign1)
                       Pmatrix2 <- Pmatrix2%*%diag(corsign2)
                       
                       ## Structural model
                       Tmatrix1 <- result1$scores
                       Tmatrix2 <- result2$scores
                       Tmatrix1 <- Tmatrix1[,perm[index1,]]%*%diag(corsign1) # account for permutation + reflection
                       Tmatrix2 <- Tmatrix2[,perm[index2,]]%*%diag(corsign2)
                       
                       #beta1 <- lm.beta(lm(data1[,(K*3+1)] ~ Tmatrix1))
                       #beta2 <- lm.beta(lm(data2[,(K*3+1)] ~ Tmatrix2))
                       
                       ## measure summary
                       measure_summary <- matrix(c('Factor', Out2$s_size, Out2$n_items, corrate_factor, tucongr_factor,
                                                   'Composite', Out1$s_size, Out1$n_items, corrate_composite, tucongr_composite), nrow = 2, byrow = T)
                       colnames(measure_summary) <- c("Model", "I", "K", "Zero_Rec", "Tucker")
                       # Saving performance measures
                       output = list(Pmatrix_factor = Pmatrix2, Tmatrix_factor = Tmatrix2,
                                     Pmatrix_composite = Pmatrix1, Tmatrix_composite = Tmatrix1,
                                     measures = measure_summary, Pfactor = Ptrue2, Pcomposite = Ptrue1)
                       
                       save(output, file = paste0("CORRELATED/IS/USLPCA", i, ".RData"))
                       
                     }
# stop Cluster
stopCluster(c1)