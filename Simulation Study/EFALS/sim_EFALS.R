#### Title: Regularized LS-LV
#### Author: Tra Le
#### Created: April 10, 2023
#### Last modified: April 11, 2023

#########################################################################################
######                         CCLSLV simulation analysis                           #####
#########################################################################################

# Clear workspace
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
rm(list = ls())

# load simulation setup
load("ORTHOGONAL/Info_simulation.RData")
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
                     .packages = c("lavaan", "lm.beta"),
                     .combine = rbind)%dopar%{
                       
                       # Load data
                       load(paste0("ORTHOGONAL2/factor",i,".RData")) #factor-based
                       Out = list(data = out$data, Y = out$Y, Ptrue = out$Ptrue, 
                                   n_items = out$n_items, s_size = out$s_size)
                       
                       # population values
                       Ptrue <- Out$Ptrue
                       I <- Out$s_size 
                       K <- Out$n_items
                       data <- scale(Out$data)
                       ## Measurement model
                       model1 <- MULTISTART_USLPCA(DATA = data[,1:(K*3)], R = 3, CARD = rep(K,3), 
                                                   MaxIter = 200, eps = 10^-6, nstarts = 10)
                       model2 <- MULTISTART_EFALS(DATA = data[,1:(K*3)], R = 3, CARD = rep(K,3), 
                                                  MaxIter = 200, eps = 10^-6, nstarts = 10)
                       
                       
                       Pmatrix1 <- model1$loadings
                       Pmatrix2 <- model2$loadings
                       
                       # zero/non-zero structure recovery rate + Tucker's congruence
                       perm <- gtools::permutations(3, 3)
                       corrate1 <- vector(length = nrow(perm))
                       corrate2 <- vector(length = nrow(perm))
                       spcr_tucongrT1 <- vector(length = nrow(perm))
                       spcr_tucongrT2 <- vector(length = nrow(perm))
                       for (p in 1:nrow(perm)) {
                         corrate1[p] <-num_correct(Ptrue, round(Pmatrix1[,perm[p,]]))
                         corrate2[p] <-num_correct(Ptrue, round(Pmatrix2[,perm[p,]]))
                         spcr_tucongrT1[p] <- sum(diag(abs(psych::factor.congruence(Ptrue,Pmatrix1[,perm[p,]]))))/3
                         spcr_tucongrT2[p] <- sum(diag(abs(psych::factor.congruence(Ptrue,Pmatrix2[,perm[p,]]))))/3
                       }
                       corrateUSLPCA <- max(corrate1)
                       corrateEFALS <- max(corrate2)
                       tucongrUSLPCA <- max(spcr_tucongrT1)
                       tucongrEFALS <- max(spcr_tucongrT2)
                       
                       index1 <- which.max(corrate1)
                       index2 <- which.max(corrate2)
                       
                       Pmatrix1 <- Pmatrix1[,perm[index1,]]
                       Pmatrix2 <- Pmatrix2[,perm[index2,]]
                       corsign1 <- sign(diag(cor(Ptrue1, Pmatrix1)))
                       corsign2 <- sign(diag(cor(Ptrue2, Pmatrix2)))
                       
                       Pmatrix1 <- Pmatrix1%*%diag(corsign1)
                       Pmatrix2 <- Pmatrix2%*%diag(corsign2)
                       
                       ## Structural model
                       Tmatrix1 <- model1$scores
                       Tmatrix2 <- model2$scores
                       Tmatrix1 <- Tmatrix1[,perm[index1,]]%*%diag(corsign1) # account for permutation + reflection
                       Tmatrix2 <- Tmatrix2[,perm[index2,]]%*%diag(corsign2)
                       
                       beta1 <- lm(data1[,(K*3+1)] ~ 0 + Tmatrix1)
                       beta2 <- lm(data1[,(K*3+1)] ~ 0 + Tmatrix2)
                       
                       ## measure summary
                       measure_summary <- matrix(c('EFALS', Out2$s_size, Out2$n_items, corrateEFALS, tucongrEFALS,
                                                   'USLPCA', Out1$s_size, Out1$n_items, corrateUSLPCA, tucongrUSLPCA), nrow = 2, byrow = T)
                       colnames(measure_summary) <- c("Method", "I", "K", "Zero_Rec", "Tucker")
                       # Saving performance measures
                       output = list(Pmat_EFALS = Pmatrix2, Tmat_EFALS = Tmatrix2,
                                     Pmat_USLPCA = Pmatrix1, Tmat_USLPCA = Tmatrix1,
                                     betaUSLPCA = beta1, betaEFALS = beta2,
                                     measures = measure_summary, Ptrue = Ptrue)
                       
                       save(output, file = paste0("ORTHOGONAL2/USLPCA", i, ".RData"))
                       
                     }
# stop Cluster
stopCluster(c1)
