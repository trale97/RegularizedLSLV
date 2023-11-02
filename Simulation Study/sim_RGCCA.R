#### Title: Simulation study
#### Author: Tra Le
#### Created: May 31, 2023
#### Last modified: 

#########################################################################################
######                          RGCCA simulation analysis                           #####
#########################################################################################

# Clear workspace
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
rm(list = ls())
#library(dplyr)
library(doParallel)

# load simulation setup
load("CORRELATED/Info_simulation.RData")
Info_matrix = Infor_simulation$design_matrix_replication
Ndatasets = Infor_simulation$n_data_sets
MatrixInfo = Infor_simulation$design_matrix_replication

# register cores
no_cores <- detectCores()
c1 <- makePSOCKcluster(floor( no_cores*.4))
registerDoParallel(c1)

# run simulation
Simulation3 = foreach(i = 1:Ndatasets,
                      .packages = c("RGCCA", "lm.beta", "dplyr"),
                      .combine = rbind)%dopar%{
                        source('num_correct.R', local = T)
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
                        Wtrue <- Out1$W
                        Ptrue2 <- Out2$Ptrue
                        
                        blocks1 <- list(bl1 = data1, bl2 = data1, bl3 = data1)
                        blocks2 <- list(bl1 = data2, bl2 = data2, bl3 = data2)
                        
                        sparsity1 <- as.data.frame(matrix(c(.5, .45, .38, .6, .5, .4, .6, .5, .4), ncol = 3, byrow = T))
                        sparsity2 <- as.data.frame(matrix(c(.57, .5, .38, .6, .5, .4, .6, .5, .4), ncol = 3, byrow = T))
                        sparsity1 <- sparsity1 %>% slice(rep(1:n(), each = 300))
                        sparsity2 <- sparsity2 %>% slice(rep(1:n(), each = 300))
                        
                        fit_sgcca1 <- rgcca(blocks = blocks1, connection =  1-diag(3), method = "sgcca", sparsity = sparsity1[i,], superblock = FALSE, ncomp = c(1,1,1),
                                           scheme = "factorial",
                                           comp_orth = TRUE,
                                           verbose = TRUE)
                        fit_sgcca2 <- rgcca(blocks = blocks2, connection =  1-diag(3), method = "sgcca", sparsity = sparsity2[i,], superblock = FALSE, ncomp = c(1,1,1),
                                            scheme = "factorial",
                                            comp_orth = TRUE,
                                            verbose = TRUE)
                        
                        w1 <- cbind(fit_sgcca1$astar$bl1, fit_sgcca1$astar$bl2, fit_sgcca1$astar$bl3)
                        w2 <- cbind(fit_sgcca2$astar$bl1, fit_sgcca2$astar$bl2, fit_sgcca2$astar$bl3)
                        Tmatrix1 <- data1%*%w1
                        Tmatrix2 <- data2%*%w2
                        Pmatrix1 <- t(solve(t(Tmatrix1)%*%Tmatrix1)%*%t(Tmatrix1)%*%data1)
                        Pmatrix2 <- t(solve(t(Tmatrix2)%*%Tmatrix2)%*%t(Tmatrix2)%*%data2)
                      
                        
                        
                        #LV1_1 <- plot(fit_sgcca1, type = "loadings", block = 1, comp = 1, display_order = FALSE)
                        #LV2_1 <- plot(fit_sgcca1, type = "loadings", block = 1, comp = 2, display_order = FALSE)
                        #LV3_1 <- plot(fit_sgcca1, type = "loadings", block = 1, comp = 3, display_order = FALSE)
                        
                        #LV1_2 <- plot(fit_sgcca2, type = "loadings", block = 1, comp = 1, display_order = FALSE)
                        #LV2_2 <- plot(fit_sgcca2, type = "loadings", block = 1, comp = 2, display_order = FALSE)
                        #LV3_2 <- plot(fit_sgcca2, type = "loadings", block = 1, comp = 3, display_order = FALSE)
                        
                        #Pmatrix1 <- cbind(LV1_1$data$x, LV2_1$data$x, LV3_1$data$x)
                        #Pmatrix2 <- cbind(LV1_2$data$x, LV2_2$data$x, LV3_2$data$x)
                        #dim(Pmatrix1)
                        #dim(Pmatrix2)
                        
                        # zero/non-zero structure recovery rate + Tucker's congruence
                        perm <- gtools::permutations(3, 3)
                        corrate1 <- vector(length = nrow(perm))
                        corrate2 <- vector(length = nrow(perm))
                        spcr_tucongrT1 <- vector(length = nrow(perm))
                        spcr_tucongrT2 <- vector(length = nrow(perm))
                        for (p in 1:nrow(perm)) {
                          corrate1[p] <- num_correct(Wtrue, round(w1[,perm[p,]],1))
                          corrate2[p] <- num_correct(Ptrue2, round(Pmatrix2[,perm[p,]],1))
                          spcr_tucongrT1[p] <- sum(diag(abs(psych::factor.congruence(Wtrue,w1[,perm[p,]]))))/3
                          spcr_tucongrT2[p] <- sum(diag(abs(psych::factor.congruence(Ptrue2,Pmatrix2[,perm[p,]]))))/3
                        }
                        corrate_composite <- max(corrate1)
                        corrate_factor <- max(corrate2)
                        tucongr_composite <- max(spcr_tucongrT1)
                        tucongr_factor <- max(spcr_tucongrT2)
                        
                        index1 <- which.max(spcr_tucongrT1)
                        index2 <- which.max(spcr_tucongrT2)
                        
                        w1 <- w1[,perm[index1,]]
                        Pmatrix2 <- Pmatrix2[,perm[index2,]]
                        corsign1 <- sign(diag(cor(Wtrue, w1)))
                        corsign2 <- sign(diag(cor(Ptrue2, Pmatrix2)))
                        
                        w1 <- w1%*%diag(corsign1)
                        Pmatrix2 <- Pmatrix2%*%diag(corsign2)
                        
                        ## Structural model
                        Tmatrix1 <- Tmatrix1[,perm[index1,]]%*%diag(corsign1) # account for permutation + reflection
                        Tmatrix2 <- Tmatrix2[,perm[index2,]]%*%diag(corsign2)
                        
                        beta1 <- lm.beta(lm(Tmatrix1[,3] ~ Tmatrix1[,1] + Tmatrix1[,2]))
                        beta2 <- lm.beta(lm(Tmatrix2[,3] ~ Tmatrix2[,1] + Tmatrix2[,2]))
                        #beta1 <- c(beta1$standardized.coefficients[[2]], beta1$standardized.coefficients[[3]])
                        #beta2 <- c(beta2$standardized.coefficients[[2]], beta2$standardized.coefficients[[3]])
                        
                        ## measure summary
                        measure_summary <- matrix(c('Factor', Out2$s_size, Out2$n_items, corrate_factor, tucongr_factor,
                                                    'Composite', Out1$s_size, Out1$n_items, corrate_composite, tucongr_composite), nrow = 2, byrow = T)
                        colnames(measure_summary) <- c("Model", "I", "K", "Zero_Rec", "Tucker")
                        # Saving performance measures
                        output = list(Pmatrix_factor = Pmatrix2, Tmatrix_factor = Tmatrix2,
                                      Wmatrix_composite = w1, Tmatrix_composite = Tmatrix1,
                                      beta_factor = beta2, beta_composite = beta1,
                                      measures = measure_summary, Pfactor = Ptrue2, Wcomposite = Wtrue, Wfactor = w2)
                        
                        save(output, file = paste0("RGCCA_new/CORRELATED/RGCCA", i, ".RData"))
                        
                      }
# stop Cluster
stopCluster(c1)

