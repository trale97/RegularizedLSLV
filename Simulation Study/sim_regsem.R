# models for regsem
model3_regsem<- '
f1 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
f2 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
f3 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9

f1 ~~ 1*f1
f2 ~~ 1*f2
f3 ~~ 1*f3

f1 ~~ 0*f2
f1 ~~ 0*f3
f2 ~~ 0*f3

V10 ~ f1 + f2 + f3'

model5_regsem <- '
f1 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15
f2 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15
f3 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15

f1 ~~ 1*f1
f2 ~~ 1*f2
f3 ~~ 1*f3

f1 ~~ 0*f2
f1 ~~ 0*f3
f2 ~~ 0*f3

V16 ~ f1 + f2 + f3'

model7_regsem <- '
f1 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + V21
f2 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + V21
f3 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + V21

f1 ~~ 1*f1
f2 ~~ 1*f2
f3 ~~ 1*f3

f1 ~~ 0*f2
f1 ~~ 0*f3
f2 ~~ 0*f3

V22 ~ f1 + f2 + f3'

# register cores
`no_cores <- detectCores()
c1 <- makePSOCKcluster(floor( no_cores*.4))
registerDoParallel(c1)

# run simulation
Simulation_regSEM3 = foreach(i = 1:300,
                             .packages = c("regsem"),
                             .combine = rbind)%dopar%{
                               
                               # Load data
                               load(paste0("ORTHOGONAL2/composite",i,".RData")) # composite-based
                               Out1 = list(data = out$data, W = out$W, Ptrue = out$Ptrue, 
                                           n_items = out$n_items, s_size = out$s_size)
                               
                               
                               load(paste0("ORTHOGONAL2/factor",i,".RData")) #factor-based
                               Out2 = list(data = out$data, Ptrue = out$Ptrue, 
                                           n_items = out$n_items, s_size = out$s_size)
                               
                               I <- Out1$s_size 
                               K <- Out1$s_size
                               Ptrue1 <- Out1$Ptrue
                               Ptrue2 <- Out2$Ptrue
                               data1 <- scale(Out1$data)
                               data2 <- scale(Out2$data)
                               
                               # regsem
                               fit1 <- sem(model3_regsem, data = data1, se = "none")
                               a <- extractMatrices(fit1)$A
                               reg.out1 <- cv_regsem(model = fit1, n.lambda = 50, type = "lasso", pars_pen = a[1:9, 11:13], jump = .01, lambda.start = .01, mult.start = T)
                               loadings1 <- unname(reg.out1$final_pars[1:27])
                               Pmatrix1 <- matrix(c(loadings1), ncol = 3)
                               
                               fit2 <- sem(model3_regsem, data = data2, se = "none")
                               reg.out2 <- cv_regsem(model = fit2, n.lambda = 50, type = "lasso", pars_pen = a[1:9, 11:13], jump = .001, lambda.start = .01, mult.start = T)
                               loadings2 <- unname(reg.out2$final_pars[1:27])
                               Pmatrix2 <- matrix(c(loadings2), ncol = 3)
                               
                               # account for permutation + reflection
                               perm <- gtools::permutations(3, 3)
                               corrate1 <- vector(length = nrow(perm))
                               corrate2 <- vector(length = nrow(perm))
                               spcr_tucongrT1 <- vector(length = nrow(perm))
                               spcr_tucongrT2 <- vector(length = nrow(perm))
                               for (p in 1:nrow(perm)) {
                                 corrate1[p] <-num_correct(Ptrue1, Pmatrix1[,perm[p,]])
                                 corrate2[p] <-num_correct(Ptrue2, Pmatrix2[,perm[p,]])
                                 spcr_tucongrT1[p] <- sum(diag(abs(psych::factor.congruence(Ptrue1,Pmatrix1[,perm[p,]]))))/3
                                 spcr_tucongrT2[p] <- sum(diag(abs(psych::factor.congruence(Ptrue2,Pmatrix2[,perm[p,]]))))/3
                               }
                               
                               Pmatrix1 <- Pmatrix1[,perm[which.max(spcr_tucongrT1),]]
                               Pmatrix2 <- Pmatrix2[,perm[which.max(spcr_tucongrT2),]]
                               corsign1 <- sign(diag(cor(Ptrue1, Pmatrix1)))
                               corsign2 <- sign(diag(cor(Ptrue2, Pmatrix2)))
                               
                               Pmatrix1 <- Pmatrix1%*%diag(corsign1)
                               Pmatrix2 <- Pmatrix2%*%diag(corsign2)
                               
                               ## Structural model
                               beta1 <- unname(reg.out1$final_pars[28:30])[perm[which.max(spcr_tucongrT1),]]
                               beta2 <- unname(reg.out2$final_pars[28:30])[perm[which.max(spcr_tucongrT2),]]
                               beta1 <- beta1*corsign1
                               beta2 <- beta2*corsign2
                               
                               # zero/non-zero structure recovery rate + Tucker's congruen
                               corrate_composite <- max(corrate1)
                               corrate_factor <- max(corrate2)
                               tucongr_composite <- max(spcr_tucongrT1)
                               tucongr_factor <- max(spcr_tucongrT2)
                               
                               ## measure summary
                               measure_summary <- matrix(c('Factor', Out2$s_size, Out2$n_items, corrate_factor, tucongr_factor,
                                                           'Composite', Out1$s_size, Out1$n_items, corrate_composite, tucongr_composite), nrow = 2, byrow = T)
                               colnames(measure_summary) <- c("Model", "I", "K", "Zero_Rec", "Tucker")
                               
                               # Saving performance measures
                               output = list(Pmatrix_factor = Pmatrix2, Ptrue_factor = Ptrue2,
                                             Pmatrix_composite = Pmatrix1, Ptrue_composite = Ptrue1,
                                             beta_composite = beta1, beta_factor = beta2,
                                             measures = measure_summary)
                               
                               save(output, file = paste0("ORTHOGONAL2/regSEM", i, ".RData"))
                               
                             }
# stop Cluster
stopCluster(c1)        

#### 5 items
rm(list = ls())
no_cores <- detectCores()
c1 <- makePSOCKcluster(floor( no_cores*.4))
registerDoParallel(c1)
Simulation_regSEM5 = foreach(i = 301:600,
                             .packages = c("regsem"),
                             .combine = rbind)%dopar%{
                               
                               # Load data
                               load(paste0("ORTHOGONAL2/composite",i,".RData")) # composite-based
                               Out1 = list(data = out$data, W = out$W, Ptrue = out$Ptrue, 
                                           n_items = out$n_items, s_size = out$s_size)
                               
                               
                               load(paste0("ORTHOGONAL2/factor",i,".RData")) #factor-based
                               Out2 = list(data = out$data, Ptrue = out$Ptrue, 
                                           n_items = out$n_items, s_size = out$s_size)
                               
                               I <- Out1$s_size 
                               K <- Out1$s_size
                               Ptrue1 <- Out1$Ptrue
                               Ptrue2 <- Out2$Ptrue
                               data1 <- scale(Out1$data)
                               data2 <- scale(Out2$data)
                               
                               # regsem
                               fit1 <- sem(model5_regsem, data = data1, se = "none")
                               a <- extractMatrices(fit1)$A
                               reg.out1 <- cv_regsem(model = fit1, n.lambda = 50, type = "lasso", pars_pen = a[1:15, 17:19], jump = .01, lambda.start = .01, mult.start = T)
                               loadings1 <- unname(reg.out1$final_pars[1:45])
                               Pmatrix1 <- matrix(c(loadings1), ncol = 3)
                               
                               fit2 <- sem(model5_regsem, data = data2, se = "none")
                               reg.out2 <- cv_regsem(model = fit2, n.lambda = 50, type = "lasso", pars_pen = a[1:15, 17:19], jump = .001, lambda.start = .01, mult.start = T)
                               loadings2 <- unname(reg.out2$final_pars[1:45])
                               Pmatrix2 <- matrix(c(loadings2), ncol = 3)
                               
                               # account for permutation + reflection
                               perm <- gtools::permutations(3, 3)
                               corrate1 <- vector(length = nrow(perm))
                               corrate2 <- vector(length = nrow(perm))
                               spcr_tucongrT1 <- vector(length = nrow(perm))
                               spcr_tucongrT2 <- vector(length = nrow(perm))
                               for (p in 1:nrow(perm)) {
                                 corrate1[p] <-num_correct(Ptrue1, Pmatrix1[,perm[p,]])
                                 corrate2[p] <-num_correct(Ptrue2, Pmatrix2[,perm[p,]])
                                 spcr_tucongrT1[p] <- sum(diag(abs(psych::factor.congruence(Ptrue1,Pmatrix1[,perm[p,]]))))/3
                                 spcr_tucongrT2[p] <- sum(diag(abs(psych::factor.congruence(Ptrue2,Pmatrix2[,perm[p,]]))))/3
                               }
                               
                               Pmatrix1 <- Pmatrix1[,perm[which.max(spcr_tucongrT1),]]
                               Pmatrix2 <- Pmatrix2[,perm[which.max(spcr_tucongrT2),]]
                               corsign1 <- sign(diag(cor(Ptrue1, Pmatrix1)))
                               corsign2 <- sign(diag(cor(Ptrue2, Pmatrix2)))
                               
                               Pmatrix1 <- Pmatrix1%*%diag(corsign1)
                               Pmatrix2 <- Pmatrix2%*%diag(corsign2)
                               
                               ## Structural model
                               beta1 <- unname(reg.out1$final_pars[46:48])[perm[which.max(spcr_tucongrT1),]]
                               beta2 <- unname(reg.out2$final_pars[46:48])[perm[which.max(spcr_tucongrT2),]]
                               beta1 <- beta1*corsign1
                               beta2 <- beta2*corsign2
                               
                               # zero/non-zero structure recovery rate + Tucker's congruen
                               corrate_composite <- max(corrate1)
                               corrate_factor <- max(corrate2)
                               tucongr_composite <- max(spcr_tucongrT1)
                               tucongr_factor <- max(spcr_tucongrT2)
                              
                               ## measure summary
                               measure_summary <- matrix(c('Factor', Out2$s_size, Out2$n_items, corrate_factor, tucongr_factor,
                                                           'Composite', Out1$s_size, Out1$n_items, corrate_composite, tucongr_composite), nrow = 2, byrow = T)
                               colnames(measure_summary) <- c("Model", "I", "K", "Zero_Rec", "Tucker")
                               
                               # Saving performance measures
                               output = list(Pmatrix_factor = Pmatrix2, Ptrue_factor = Ptrue2,
                                             Pmatrix_composite = Pmatrix1, Ptrue_composite = Ptrue1,
                                             beta_composite = beta1, beta_factor = beta2,
                                             measures = measure_summary)
                               
                               save(output, file = paste0("ORTHOGONAL2/regSEM", i, ".RData"))
                               
                             }
# stop Cluster
stopCluster(c1)

#### 7 items
rm(list = ls())
