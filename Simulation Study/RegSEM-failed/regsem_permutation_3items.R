rm(list = ls())

model3_regsem <- '
f1 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
f2 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
f3 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9

f1 ~~ 1*f1
f2 ~~ 1*f2
f3 ~~ 1*f3

f1 ~~ f2
f3 ~ f1 + f2'

model3_regsem2 <- '
f1 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
f2 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
f3 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9

f1 ~~ 1*f1
f2 ~~ 1*f2
f3 ~~ 1*f3

f1 ~~ f2
f3 ~ f2 + f1'

model3_regsem3 <- '
f1 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
f2 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
f3 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9

f1 ~~ 1*f1
f2 ~~ 1*f2
f3 ~~ 1*f3

f1 ~~ f3
f2 ~ f1 + f3'

model3_regsem4 <- '
f1 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
f2 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
f3 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9

f1 ~~ 1*f1
f2 ~~ 1*f2
f3 ~~ 1*f3

f1 ~~ f3
f2 ~ f3 + f1'

model3_regsem5 <- '
f1 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
f2 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
f3 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9

f1 ~~ 1*f1
f2 ~~ 1*f2
f3 ~~ 1*f3

f2 ~~ f3
f1 ~ f2 + f3'

model3_regsem6 <- '
f1 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
f2 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
f3 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9

f1 ~~ 1*f1
f2 ~~ 1*f2
f3 ~~ 1*f3

f2 ~~ f3
f1 ~ f3 + f2'

model3_list <- list(model3_regsem, model3_regsem2, model3_regsem3, model3_regsem4, model3_regsem5, model3_regsem6)

# register cores
library(doParallel)
no_cores <- detectCores()
c1 <- makePSOCKcluster(floor( no_cores*.4))
registerDoParallel(c1)

Simulation_regSEM7 = foreach(i = 701:900,
                             .packages = c("regsem"),
                             .combine = rbind)%dopar%{
                               source('findLasso_regSEM.R', local = T)
                               
                               # Load data
                               load(paste0("CORRELATED/composite",i,".RData")) # composite-based
                               Out1 = list(data = out$data, W = out$W, Ptrue = out$Ptrue, 
                                           n_items = out$n_items, s_size = out$s_size)
                               
                               
                               load(paste0("CORRELATED/factor",i,".RData")) #factor-based
                               Out2 = list(data = out$data, Ptrue = out$Ptrue, 
                                           n_items = out$n_items, s_size = out$s_size)
                               
                               I <- Out1$s_size 
                               K <- Out1$n_items
                               Ptrue1 <- Out1$Ptrue
                               Ptrue2 <- Out2$Ptrue
                               data1 <- scale(Out1$data)
                               data2 <- scale(Out2$data)
                               
                               # regsem
                               Pmatrix1_list <- list()
                               Pmatrix2_list <- list()
                               beta1_list <- list()
                               beta2_list <- list()
                               
                               corrate1 <- vector(length = length(model3_list))
                               corrate2 <- vector(length = length(model3_list))
                               spcr_tucongrT1 <- vector(length = length(model3_list))
                               spcr_tucongrT2 <- vector(length = length(model3_list))
                               
                               
                               for (j in 1:length(model3_list)){
                                 fit1 <- sem(model3_list[[j]], data = data1, se = "none")
                                 a <- extractMatrices(fit1)$A
                                 #lasso_para1 <- findLasso(fit = fit1, Ptrue = Ptrue1, maxItr = 1000, lassou = 0.1)
                                 reg.out1 <- cv_regsem(model = fit1, n.lambda = 50, type = "lasso", pars_pen = a[1:9, 10:12], jump = .01, lambda.start = .1, mult.start = T)
                                 loadings1 <- unname(reg.out1$final_pars[1:27])
                                 Pmatrix1 <- matrix(c(loadings1), ncol = 3)
                                 Pmatrix1_list[[j]] <- Pmatrix1
                                 beta1 <- unname(reg.out1$final_pars[28:30])
                                 beta1_list[[j]] <- beta1
                                 
                                 fit2 <- sem(model3_list[[j]], data = data2, se = "none")
                                 #lasso_para2 <- findLasso(fit = fit2, Ptrue = Ptrue2, maxItr = 1000, lassou = 0.1)
                                 reg.out2 <- cv_regsem(model = fit2, n.lambda = 50, type = "lasso", pars_pen = a[1:9, 10:12], jump = .01, lambda.start = .1, mult.start = T)
                                 loadings2 <- unname(reg.out2$final_pars[1:27])
                                 Pmatrix2 <- matrix(c(loadings2), ncol = 3)
                                 Pmatrix2_list[[j]] <- Pmatrix2
                                 beta2 <- unname(reg.out2$final_pars[28:30])
                                 beta2_list[[j]] <- beta2
                                 
                                 corrate1[j] <-num_correct(Ptrue1, Pmatrix1)
                                 corrate2[j] <-num_correct(Ptrue2, Pmatrix2)
                                 spcr_tucongrT1[j] <- sum(diag(abs(psych::factor.congruence(Ptrue1,Pmatrix1))))/3
                                 spcr_tucongrT2[j] <- sum(diag(abs(psych::factor.congruence(Ptrue2,Pmatrix2))))/3
                               }
                               
                               # zero/non-zero structure recovery rate + Tucker's congruen
                               corrate_composite <- max(corrate1, na.rm = T)
                               corrate_factor <- max(corrate2, na.rm = T)
                               tucongr_composite <- max(spcr_tucongrT1, na.rm = T)
                               tucongr_factor <- max(spcr_tucongrT2, na.rm = T)
                               
                               Pmatrix1 <- Pmatrix1_list[[which.max(spcr_tucongrT1)]]
                               Pmatrix2 <- Pmatrix1_list[[which.max(spcr_tucongrT2)]]
                               
                               # structural model
                               beta1 <- beta1_list[[which.max(spcr_tucongrT1)]]
                               beta2 <- beta2_list[[which.max(spcr_tucongrT1)]]
                               
                               
                               ## measure summary
                               measure_summary <- matrix(c('Factor', Out2$s_size, Out2$n_items, corrate_factor, tucongr_factor,
                                                           'Composite', Out1$s_size, Out1$n_items, corrate_composite, tucongr_composite), nrow = 2, byrow = T)
                               colnames(measure_summary) <- c("Model", "I", "K", "Zero_Rec", "Tucker")
                               
                               # Saving performance measures
                               output = list(Pmatrix_factor = Pmatrix2, Ptrue_factor = Ptrue2,
                                             Pmatrix_composite = Pmatrix1, Ptrue_composite = Ptrue1,
                                             beta_composite = beta1, beta_factor = beta2,
                                             measures = measure_summary)
                               
                               save(output, file = paste0("CORRELATED/IS/Nov24/regSEM", i, ".RData"))
                               
                             }

# stop Cluster
stopCluster(c1)  