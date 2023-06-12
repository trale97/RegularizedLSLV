setwd("C:/Users/onno1234/Documents/Thesis")
dir.create("CORRELATED")

library(doParallel)
library(dplyr)
n <- c(20, 100, 500) # sample size
k <- c(3, 5, 7)

design_matrix <- expand.grid(s_size = n, n_items = k)
design_matrix_replication <- design_matrix %>% slice(rep(1:n(), each = 100))
Infor_simulation = list(n_data_sets = nrow(design_matrix_replication), n_replications = 100,
                        design_matrix_replication = design_matrix_replication)
save(Infor_simulation, file = "CORRELATED/Info_simulation.RData")

# setup
list1 <- replicate(300,corr3, simplify = F)
list2 <- replicate(300, corr5, simplify = F)
list3 <- replicate(300, corr7, simplify = F)
corrmat <- append(append(list1, list2), list3)

Ptrue3 <- replicate(300, Ptrue3, simplify = F)
Ptrue5 <- replicate(300, Ptrue5, simplify = F)
Ptrue7 <- replicate(300, Ptrue7, simplify = F)
Ptruelist <- append(append(Ptrue3,Ptrue5), Ptrue7)

W3 <- replicate(300, w3, simplify = F)
W5 <- replicate(300, w5, simplify = F)
W7 <- replicate(300, w7, simplify = F)
Wlist <- append(append(W3, W5), W7)

betas <- c(-0.02, 0.15, .5)

# register cores
no_cores <- detectCores() -1
c1 <- makePSOCKcluster(no_cores)
registerDoParallel(c1)

# start simulating the data
results_sim1_data1 <- foreach(i = 1:nrow(design_matrix_replication),
                              .options.RNG = 2197,
                              .packages = c("MASS"),
                              .combine = rbind) %dopar%{
                                # set the specific values of the parameters 
                                I = design_matrix_replication$s_size[i]
                                K = design_matrix_replication$n_items[i]
                                corr = corrmat[[i]]
                                Ptrue = Ptruelist[[i]]
                                W = Wlist[[i]]
                                
                                data <- mvrnorm(n = I, mu = rep(0, K*3), Sigma = corr, empirical = F)
                                
                                out <- list(data = data, Ptrue = Ptrue, W = W, corr = corr,
                                            s_size = I, n_items = K, teller = i)
                                save(out, file = paste0("CORRELATED/composite",i,".RData"))  
                              }

stopCluster(c1)

## Factor-based ----------------------------------------------------------------
rm(list = ls())
design_matrix_factor <- expand.grid(s_size = n, n_items = k)
design_matrix_replication_factor <- design_matrix_factor %>% slice(rep(1:n(), each = 100))
Infor_simulation_factor = list(n_data_sets = nrow(design_matrix_replication_factor), n_replications = 100,
                               design_matrix_replication_factor = design_matrix_replication_factor)
save(Infor_simulation_factor, file = "ORTHOGONAL/Info_simulation_factor.RData")

# Setup
corr3f <- replicate(300,corr3_f, simplify = F)
corr5f <- replicate(300, corr5_f, simplify = F)
corr7f <- replicate(300, corr7_f, simplify = F)
corrmat <- append(append(corr3f, corr5f), corr7f)

Ptrue3f <- matrix(c(.7, .8, .9, rep(0,9), .7, .8, .9, rep(0,9), .7, .8, .9), ncol = 3)
Ptrue5f <- matrix(c(.6, .6, .7, .8, .8, rep(0, 15), .6, .6, .7, .8, .8, rep(0, 15), .6, .6, .7, .8, .8), ncol = 3)
Ptrue7f <- matrix(c(.6, .6, .7, .7, .7, .8, .8, rep(0,21), 
                    .6, .6, .7, .7, .7, .8, .8, rep(0,21), .6, .6, .7, .7, .7, .8, .8), ncol = 3)
Ptrue3f <- replicate(300, Ptrue3_f, simplify = F)
Ptrue5f <- replicate(300, Ptrue5_f, simplify = F)
Ptrue7f <- replicate(300, Ptrue7_f, simplify = F)
PtrueList <- append(append(Ptrue3f, Ptrue5f), Ptrue7f)


# register cores
no_cores <- detectCores() -1
c1 <- makePSOCKcluster(no_cores)
registerDoParallel(c1)

results_sim2_data2 <- foreach(i = 1:nrow(design_matrix_replication),
                              .options.RNG = 2197,
                              .packages = c("MASS"),
                              .combine = rbind) %dopar%{
                                # set the specific values of the parameters 
                                I = design_matrix_replication$s_size[i]
                                K = design_matrix_replication$n_items[i]
                                corr = corrmat[[i]]
                                Ptrue = PtrueList[[i]]
                                
                                data <- mvrnorm(n = I, mu = rep(0, K*3), Sigma = corr, empirical = F)
                                
                                out <- list(data = data, Ptrue = Ptrue, corr = corr,
                                            s_size = I, n_items = K, teller = i)
                                save(out, file = paste0("CORRELATED/factor",i,".RData"))  
                              }

stopCluster(c1)


