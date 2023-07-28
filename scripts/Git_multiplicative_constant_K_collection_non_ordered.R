## script to launch simulations
## author : Perrine Lacroix
## date : May, 23 2023

## This code creates the model collections where the ordered on variables in not respected and characteristics of each model are computed. 
## The used data sets are the same as for the deterministic ordered model collections. 


rm(list=objects())

# Packages loading
library(bolasso)
library(SLOPE)
library(knockoff)
library(caret)
library(doMC)
library(xtable)

prefix = "~/Documents/GitHub/Trade_off_FDR_PR"
setwd(paste(prefix,"/data",sep =""))
# Parameters of simulation
nbr_it = 1000   # Data sets number for the empirical estimations of FDR and PR. 
nbr_it_bound = 100   # Data sets number to compute the FDR bounds and to study the variability of the FDR bounds.
n = 50  # number of observations
p = 50   # number of variables
n_temps = 2*n*nbr_it   # Required total iterations number to get nbr_it train sets for the empirical estimation of FDR, with nbr_it validation sets for the empirical estimation of PR. 
n_temps_bounds = n*nbr_it_bound   # Required total iterations number to get nbr_it_bound data sets for the FDR bounds and the variability studies of the FDR bounds.
sigma_2 <- 1   # Variance of the noise in the Gaussian linear regression
m=50   # number of resamples for random model collections
alpha_FDR <- 0.05    # threshold for the knockoff method

## parameters of the created data for the saving
config_data = paste0(paste("sigma2",sigma_2,"p",p,"n",n,"nbr_it",nbr_it,sep="_"),".RData")


#####################################################################################################################
###### Data are the same as for deterministic ordered model collections
#####################################################################################################################
setwd(paste(prefix,"/data",sep =""))
config_data = paste0(paste("sigma2",sigma_2,
                           "p",p,"n",n,"nbr_it",nbr_it,sep="_"),".RData")
load(paste("matrix_X",config_data,sep="_"))
config_data = paste0(paste("beta_true_0_20","sigma2",sigma_2,"p",p,"n",n,"nbr_it",nbr_it,sep="_"),".RData")
load(paste("vector_Y",config_data,sep="_"))
# load(paste("smaller_noise_vector_Y",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("larger_noise_closed_values_vector_Y",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article
load(paste("vector_Y_noise_bounds",config_data,sep="_"))
# load(paste("vector_Y_smaller_noise_bounds",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("vector_Y_larger_noise_closed_values_bounds",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article



#####################################################################################################################
###### Creation of the deterministic model collections where order on variables is not respected 
#####################################################################################################################
## Here, the dimension of the true model is fixed to 10
## According to the construction of beta_0, X_1 is the most relevant variable, then X_2, then X_3, ...
## We change the order of variables by creating different deterministic model collections:
## 1) Random permutation of active variables
## 2) Exchanging variables 9-11 and 10-12      (Exchanging the least important variables)
## 3) Exchanging variables 1-16 and 2-20      (the most important variables are pushed far into the collection)
## 4) Random permutations of active variables that vary between each iteration
## 5) Random permutations of variables 1:12 that vary between each iteration
## 6) Random permutations of variables 1:15 that vary between each iteration


###### Creation of the deterministic model collections used for empirical estimations where order on variables is not respected 
set.seed(1234)    # for the randomness of permutations
dimension <- 0:p
path_collection <- list()
support <- list()   # list for the supports of variable in the collection
support_num <- list()    # list for the supports of variable index in the collection
exchang1 <- c(sample(1:10, replace = FALSE),11:p)    # Random permutation of active variables
support[[1]] <- list()
support_num[[1]] <- list()
for (k in 2:(p+1)){
  support[[1]][[k]] <- colnames(data_metric)[exchang1[1:(k-1)]]
  support_num[[1]][[k]] <- exchang1[1:(k-1)]
}
exchang2 <- c(1:8,11,10,9,12:p)    # Exchanging variables 9-11 and 10-12      (Exchanging the least important variables)
support[[2]] <- list()
support_num[[2]] <- list()
for (k in 2:(p+1)){
  support[[2]][[k]] <- colnames(data_metric)[exchang2[1:(k-1)]]
  support_num[[2]][[k]] <- exchang2[1:(k-1)]
}
exchang3 <- c(16,20,3:15,1,17:19,2,21:p)     # Exchanging variables 1-16 and 2-20      (the most important variables are pushed far into the collection)
support[[3]] <- list()
support_num[[3]] <- list()
for (k in 2:(p+1)){
  support[[3]][[k]] <- colnames(data_metric)[exchang3[1:(k-1)]]
  support_num[[3]][[k]] <- exchang3[1:(k-1)]
}
exchang4 <- lapply(1:nbr_it, function(it) c(sample(1:10, replace = FALSE),11:p))     # Random permutations of active variables that vary between each iteration
support[[4]] <- list()
support_num[[4]] <- list()
for (it in 1:nbr_it){
  support[[4]][[it]] <- list()
  support_num[[4]][[it]]  <- list()
  for (k in 2:(p+1)){
    support[[4]][[it]][[k]] <- colnames(data_metric)[exchang4[[it]][1:(k-1)]]
    support_num[[4]][[it]][[k]] <- exchang4[[it]][1:(k-1)]
  }
}
exchang5 <- lapply(1:nbr_it, function(it) c(sample(1:12, replace = FALSE),13:p))      # Random permutations of variables 1:12 that vary between each iteration
support[[5]] <- list()
support_num[[5]] <- list()
for (it in 1:nbr_it){
  support[[5]][[it]] <- list()
  support_num[[5]][[it]]  <- list()
  for (k in 2:(p+1)){
    support[[5]][[it]][[k]] <- colnames(data_metric)[exchang5[[it]][1:(k-1)]]
    support_num[[5]][[it]][[k]] <- exchang5[[it]][1:(k-1)]
  }
}
exchang6 <- lapply(1:nbr_it, function(it) c(sample(1:15, replace = FALSE),16:p))      # Random permutations of variables 1:15 that vary between each iteration
support[[6]] <- list()
support_num[[6]] <- list()
for (it in 1:nbr_it){
  support[[6]][[it]] <- list()
  support_num[[6]][[it]]  <- list()
  for (k in 2:(p+1)){
    support[[6]][[it]][[k]] <- colnames(data_metric)[exchang6[[it]][1:(k-1)]]
    support_num[[6]][[it]][[k]] <- exchang6[[it]][1:(k-1)]
  }
}
pb = txtProgressBar(min = 1, max = 6, initial = 1)
l = 11    # the dimension of the true model is fixed to 10
Y_select <- Y_select_matrix[,l]
for (x in 1:6){
  path_collection[[x]] <- list()
  for (k in seq(1,n*nbr_it,by=n)){
    j = which(seq(1,n*nbr_it,by=n) ==k)
    if (x<=3){
      support_num_bis <- support_num[[x]]
      support_bis <- support[[x]]
    }else{
      support_num_bis <- support_num[[x]][[j]]
      support_bis <- support[[x]][[j]]
    }
    set.seed(j+1234)
    index <- 2:(p+1)
    re_estimation <- list()
    for (i in index){
      re_estimation[[i]] <- lm(Y_select[k:(k+(n-1))] ~ data_select[,support_num_bis[[i]]] -1)  # No intercept, mean square estimators onto each model
      # The beta coefficients are re-estimated onto each model of the collection by the MSE minization
    }
    beta_new <- sapply(re_estimation, function(a) as.vector(a$coefficients))
    for (i in index){
      names(beta_new[[i]]) <- support_bis[[i]]
    }
    beta_estimator <- matrix(0, nrow = length(support_bis), ncol = p)
    for (i in 2:length(support_bis)){
      beta_estimator[i,support_num_bis[[i]]] <- as.numeric(beta_new[[i]])
    }
    LS = sapply(1:length(dimension), function(i)
      (1/nrow(data_select))*sum((Y_select[k:(k+(n-1))] - data_select %*% beta_estimator[i,])^2)) # The least squared values are computed for each beta of the collection
    data_frame_temp_1 <- data.frame(LS=LS,dim=dimension,complexite=lchoose(p,dimension))      # characteristics of the model collection
    data_frame_temp_2 <- data.frame(beta=beta_estimator)
    path_collection[[x]][[j]]=list(val1 = data_frame_temp_1, val2 = data_frame_temp_2)
  }
  setTxtProgressBar(pb,x)
}
save(path_collection,file=paste("path_collection_deterministic",config_data,sep="_"))
# save(path_collection,file=paste("path_collection_smaller_noise_deterministic",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# save(path_collection,file=paste("path_collection_larger_noise_closed_values_deterministic",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article


###### Creation of the deterministic model collections used for FDR bound evaluations where order on variables is not respected 
path_collection_bounds <- list()
pb = txtProgressBar(min = 1, max = 6, initial = 1)
l=11      # the dimension of the true model is fixed to 10
Y_bounds <- Y_bounds_matrix[,l]
for (x in 1:6){
  path_collection_bounds[[x]] <- list()
  for (k in seq(1,n*nbr_it_bound,by=n)){
    j = which(seq(1,n*nbr_it_bound,by=n) ==k)
    if (x<=3){
      support_num_bis <- support_num[[x]]
      support_bis <- support[[x]]
    }else{
      support_num_bis <- support_num[[x]][[j]]
      support_bis <- support[[x]][[j]]
    }
    set.seed(j+1234)
    index <- 2:(p+1)
    re_estimation <- list()
    for (i in index){
      re_estimation[[i]] <- lm(Y_bounds[k:(k+(n-1))] ~ data_select[,support_num_bis[[i]]] -1)  # No intercept
      # The beta coefficients are re-estimated onto each model of the collection by the MSE minization
    }
    beta_new <- sapply(re_estimation, function(a) as.vector(a$coefficients))
    for (i in index){
      names(beta_new[[i]]) <- support_bis[[i]]
    }
    beta_estimator <- matrix(0, nrow = length(support_bis), ncol = p)
    for (i in 2:length(support_bis)){
      beta_estimator[i,support_num_bis[[i]]] <- as.numeric(beta_new[[i]])
    }
    LS = sapply(1:length(dimension), function(i)
      (1/nrow(data_select))*sum((Y_bounds[k:(k+(n-1))] - data_select %*% beta_estimator[i,])^2))  # The least squared values are computed for each beta of the collection
    data_frame_temp_1 <- data.frame(LS=LS,dim=dimension,complexite=lchoose(p,dimension))
    data_frame_temp_2 <- data.frame(beta=beta_estimator)
    path_collection_bounds[[x]][[j]]=list(val1 = data_frame_temp_1, val2 = data_frame_temp_2)
  }
  setTxtProgressBar(pb,x)
}
save(path_collection_bounds,file=paste("path_collection_bounds_deterministic",config_data,sep="_"))
# save(path_collection_bounds,file=paste("path_collection_bounds_smaller_noise_deterministic",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# save(path_collection_bounds,file=paste("path_collection_bounds_larger_noise_closed_values_deterministic",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article

#####################################################################################################################



#####################################################################################################################
###### Creation of the random model collections (order on variables is not respected)
#####################################################################################################################
## According to the construction of beta_0, X_1 is the most relevant variable, then X_2, then X_3, ...
## We generate the order of variables by creating different random model collection:
## 1) Bolasso with R function bolasso of the R package bolasso with m=50 resamples
## 2) Slope with R function SLOPE of the R package SLOPE with m=50 resamples
## 3) The Knockoff method with R function knockoff.filter of the R package knockoff with alpha=0.05
## 4) Random Forests with R function rfe of the R package caret with m=50 resamples.

###### Creation of the random model collections used for empirical estimations (order on variables is not respected)
set.seed(1234)
l = 11    # the dimension of the true model is fixed to 10
## Bolasso
Y_select <- Y_select_matrix[,l]
path_collection_bolasso <- list()
pb = txtProgressBar(min = 1, max = length(seq(1,n*nbr_it,by=n)), initial = 1)
for (k in seq(1,n*nbr_it,by=n)){
  j = which(seq(1,n*nbr_it,by=n) ==k)
  set.seed(j+1234)
  bolasso_result <- bolasso(x = data_select, y = Y_select[k:(k+(n-1))], progress = FALSE, n.boot = m, implement = 'glmnet')
  support <- lapply(bolasso_result, function(x) support <- 1*(as.matrix(x$glmnet.fit$beta)!=0))
  support_num <- lapply(1:m, function(i) apply(support[[i]], 1, sum)/dim(support[[i]])[2])
  frequency <- list()
  frequence <- matrix(rep(0,n*p),nrow=p,ncol=m)
  for (it in 1:m){
    frequence[,it] <- as.numeric(support_num[[it]])     # variable frequencies on each regularization path
  }
  frequence_vect <- apply(frequence,1,sum)/m      # total variable frequencies averaged on the resamples and the regularization paths
  order_index <- sort(frequence_vect,index.return= TRUE,decreasing = TRUE)$ix
  support_num <- list()    # list for the supports of variable index in the collection
  for(x in 2:(p+1)){
    support_num[[x]] <- order_index[1:(x-1)]     # construction of the model support
  }
  dimension <- unlist(lapply(support_num, function(x) length(x)))
  index <- 2:(p+1)
  re_estimation <- list()
  for (i in index){
    re_estimation[[i]] <- lm(Y_select[k:(k+(n-1))] ~ data_select[,support_num[[i]]] -1)  # No intercept, mean square estimators onto each model
    # The beta coefficients are re-estimated onto each model of the collection by the MSE minization
  }
  beta_new <- sapply(re_estimation, function(a) as.vector(a$coefficients))
  beta_estimator <- matrix(0, nrow = length(support_num), ncol = p)
  for (i in 2:length(support_num)){
    beta_estimator[i,support_num[[i]]] <- as.numeric(beta_new[[i]])
  }
  LS = sapply(1:length(dimension), function(i)
    (1/nrow(data_select))*sum((Y_select[k:(k+(n-1))] - data_select %*% beta_estimator[i,])^2))  # The least squared values are computed for each beta of the collection
  data_frame_temp_1 <- data.frame(LS=LS,dim=dimension,complexite=lchoose(p,dimension))      # characteristics of the model collection
  data_frame_temp_2 <- data.frame(beta=beta_estimator)
  path_collection_bolasso[[j]] <- list(val1 = data_frame_temp_1, val2 = data_frame_temp_2, val3 = order_index, val4 = frequence_vect, val5 = support_num)
  setTxtProgressBar(pb,j)
}

## Slope
Y_select <- Y_select_matrix[,l]
path_collection_slope <- list()
pb = txtProgressBar(min = 1, max = length(seq(1,n*nbr_it,by=n)), initial = 1)
for (k in seq(1,n*nbr_it,by=n)){
  j = which(seq(1,n*nbr_it,by=n) ==k)
  set.seed(j+1234)
  index_resampling <- lapply(1:m, function(x) x <- unique(sort(sample(1:n, n, replace = TRUE))))     # creation of the resamples
  sampling_regressor <- lapply(index_resampling, function(x) sampling <- data_select[x,])
  sampling_regressed <- lapply(index_resampling, function(x) sampling <- Y_select[k:(k+(n-1))][x])
  slope_results <- lapply(1:m, function(i) SLOPE(x = sampling_regressor[[i]], y = sampling_regressed[[i]], family = "gaussian", intercept = TRUE))
  matrix_freq <- matrix(0, nrow = m, ncol = p)
  for (l in 1:m){
    matrix_temp <- matrix(0, nrow = (length(slope_results[[l]]$nonzeros==TRUE)/p), ncol = p)
    for (x in seq(1,length(slope_results[[l]]$nonzeros==TRUE),by=p)){
      z = which(seq(1,length(slope_results[[l]]$nonzeros==TRUE),by=p) ==x)
      matrix_temp[z,] <- 1*slope_results[[l]]$nonzeros[x:(x+(p-1))]
    }
    matrix_freq[l,] <- apply(matrix_temp, 2, sum)/nrow(matrix_temp)     # variable frequencies on each regularization path
  }
  frequence_vect <- apply(matrix_freq,2,sum)/m      # total variable frequencies averaged on the resamples and the regularization paths
  order_index <- sort(frequence_vect,index.return= TRUE,decreasing = TRUE)$ix
  support_num <- list()    # list for the supports of variable index in the collection
  for(x in 2:(p+1)){
    support_num[[x]] <- order_index[1:(x-1)]    # construction of the model support
  }
  dimension <- unlist(lapply(support_num, function(x) length(x)))
  index <- 2:(p+1)
  re_estimation <- list()
  for (i in index){
    re_estimation[[i]] <- lm(Y_select[k:(k+(n-1))] ~ data_select[,support_num[[i]]] -1)  # No intercept, mean square estimators onto each model
    # The beta coefficients are re-estimated onto each model of the collection by the MSE minization
  }
  beta_new <- sapply(re_estimation, function(a) as.vector(a$coefficients))
  beta_estimator <- matrix(0, nrow = length(support_num), ncol = p)
  for (i in 2:length(support_num)){
    beta_estimator[i,support_num[[i]]] <- as.numeric(beta_new[[i]])
  }
  LS = sapply(1:length(dimension), function(i)
    (1/nrow(data_select))*sum((Y_select[k:(k+(n-1))] - data_select %*% beta_estimator[i,])^2))  # The least squared values are computed for each beta of the collection
  data_frame_temp_1 <- data.frame(LS=LS,dim=dimension,complexite=lchoose(p,dimension))      # characteristics of the model collection
  data_frame_temp_2 <- data.frame(beta=beta_estimator)
  path_collection_slope[[j]] <- list(val1 = data_frame_temp_1, val2 = data_frame_temp_2, val3 = order_index, val4 = frequence_vect, val5 = support_num)
  setTxtProgressBar(pb,j)
}

# Knockoffs
Y_select <- Y_select_matrix[,l]
path_collection_knockoffs <- list()
pb = txtProgressBar(min = 1, max = length(seq(1,n*nbr_it,by=n)), initial = 1)
for (k in seq(1,n*nbr_it,by=n)){
  j = which(seq(1,n*nbr_it,by=n) ==k)
  set.seed(j+1234)
  knockoffs_result <- knockoff.filter(X = data_select, y = Y_select[k:(k+(n-1))], knockoffs = create.second_order,statistic = stat.lasso_lambdasmax, fdr = alpha_FDR)
  order_index <- sort(knockoffs_result$statistic, index.return = TRUE, decreasing = TRUE)$ix
  support_num <- list()    # list for the supports of variable index in the collection
  for(x in 2:(p+1)){
    support_num[[x]] <- order_index[1:(x-1)]     # construction of the model support
  }
  dimension <- unlist(lapply(support_num, function(x) length(x)))
  index <- 2:(p+1)
  re_estimation <- list()
  for (i in index){
    re_estimation[[i]] <- lm(Y_select[k:(k+(n-1))] ~ data_select[,support_num[[i]]] -1)  # No intercept, mean square estimators onto each model
    # The beta coefficients are re-estimated onto each model of the collection by the MSE minization
  }
  beta_new <- sapply(re_estimation, function(a) as.vector(a$coefficients))
  beta_estimator <- matrix(0, nrow = length(support_num), ncol = p)
  for (i in 2:length(support_num)){
    beta_estimator[i,support_num[[i]]] <- as.numeric(beta_new[[i]])
  }
  LS = sapply(1:length(dimension), function(i)
    (1/nrow(data_select))*sum((Y_select[k:(k+(n-1))] - data_select %*% beta_estimator[i,])^2))  # The least squared values are computed for each beta of the collection
  data_frame_temp_1 <- data.frame(LS=LS,dim=dimension,complexite=lchoose(p,dimension))      # characteristics of the model collection
  data_frame_temp_2 <- data.frame(beta=beta_estimator)
  path_collection_knockoffs[[j]] <- list(val1 = data_frame_temp_1, val2 = data_frame_temp_2, val3 = order_index,  val5 = support_num)
  setTxtProgressBar(pb,j)
}

## Random forests
Y_select <- Y_select_matrix[,l]
path_collection_RF <- list()
pb = txtProgressBar(min = 1, max = length(seq(1,n*nbr_it,by=n)), initial = 1)
for (k in seq(1,n*nbr_it,by=n)){
  j = which(seq(1,n*nbr_it,by=n) ==k)
  set.seed(j+1234)
  RF_results <- rfe(x = data_select, y = Y_select[k:(k+(n-1))], sizes = p, rfeControl = rfeControl(functions = rfFuncs, number = m, method = "boot", saveDetails = TRUE), rerank = TRUE)   # RandomForests
  # the sizes parameter is the number of features that should be retained.
  # rfFuncs is the function for random forests
  # the rerank = TRUE oarameter is for recomputing the order on variables after elimination of least important variables (backward strategy)
  # the method = boot paramater is for the bootstrap strategy
  # the number parameter is the number of resamples
  order_index <- as.numeric(sapply(RF_results$optVariables, function(x) which(x==colnames(data_select))))
  support_num <- list()    # list for the supports of variable index in the collection
  for(x in 2:(p+1)){
    support_num[[x]] <- order_index[1:(x-1)]    # construction of the model support
  }
  dimension <- unlist(lapply(support_num, function(x) length(x)))
  index <- 2:(p+1)
  re_estimation <- list()
  for (i in index){
    re_estimation[[i]] <- lm(Y_select[k:(k+(n-1))] ~ data_select[,support_num[[i]]] -1)  # No intercept, mean square estimators onto each model
    # The beta coefficients are re-estimated onto each model of the collection by the MSE minization
  }
  beta_new <- sapply(re_estimation, function(a) as.vector(a$coefficients))
  beta_estimator <- matrix(0, nrow = length(support_num), ncol = p)
  for (i in 2:length(support_num)){
    beta_estimator[i,support_num[[i]]] <- as.numeric(beta_new[[i]])
  }
  LS = sapply(1:length(dimension), function(i)
    (1/nrow(data_select))*sum((Y_select[k:(k+(n-1))] - data_select %*% beta_estimator[i,])^2)) # The least squared values are computed for each beta of the collection
  data_frame_temp_1 <- data.frame(LS=LS,dim=dimension,complexite=lchoose(p,dimension))      # characteristics of the model collection
  data_frame_temp_2 <- data.frame(beta=beta_estimator)
  path_collection_RF[[j]] <- list(val1 = data_frame_temp_1, val2 = data_frame_temp_2, val3 = order_index, val5 = support_num)
  setTxtProgressBar(pb,j)
}


###### Creation of the random model collections used for FDR bound evaluations (order on variables is not respected)
set.seed(1234)
l=11   # the dimension of the true model is fixed to 10
## Bolasso
path_collection_bounds_bolasso <- list()
Y_bounds <- Y_bounds_matrix[,l]
pb = txtProgressBar(min = 1, max = length(seq(1,n*nbr_it_bound,by=n)), initial = 1)
for (k in seq(1,n*nbr_it_bound,by=n)){
  j = which(seq(1,n*nbr_it_bound,by=n) ==k)
  set.seed(j+1234)
  bolasso_result <- bolasso(x = data_select, y = Y_bounds[k:(k+(n-1))], progress = FALSE, n.boot = m, implement = 'glmnet')
  support <- lapply(bolasso_result, function(x) support <- 1*(as.matrix(x$glmnet.fit$beta)!=0))
  support_num <- lapply(1:m, function(i) apply(support[[i]], 1, sum)/dim(support[[i]])[2])
  frequency <- list()
  frequence <- matrix(rep(0,n*p),nrow=p,ncol=m)
  for (it in 1:m){
    frequence[,it] <- as.numeric(support_num[[it]])     # variable frequencies on each regularization path
  }
  frequence_vect <- apply(frequence,1,sum)/m      # total variable frequencies averaged on the resamples and the regularization paths
  order_index <- sort(frequence_vect,index.return= TRUE,decreasing = TRUE)$ix
  support_num <- list()    # list for the supports of variable index in the collection
  for(x in 2:(p+1)){
    support_num[[x]] <- order_index[1:(x-1)]     # construction of the model support
  }
  index <- 2:(p+1)
  dimension <- unlist(lapply(support_num, function(x) length(x)))
  re_estimation <- list()
  for (i in index){
    re_estimation[[i]] <- lm(Y_bounds[k:(k+(n-1))] ~ data_select[,support_num[[i]]] -1)  # No intercept, mean square estimators onto each model
    # The beta coefficients are re-estimated onto each model of the collection by the MSE minization
  }
  beta_new <- sapply(re_estimation, function(a) as.vector(a$coefficients))
  beta_estimator <- matrix(0, nrow = length(support_num), ncol = p)
  for (i in 2:length(support_num)){
    beta_estimator[i,support_num[[i]]] <- as.numeric(beta_new[[i]])
  }
  LS = sapply(1:length(dimension), function(i)
    (1/nrow(data_select))*sum((Y_bounds[k:(k+(n-1))] - data_select %*% beta_estimator[i,])^2)) # The least squared values are computed for each beta of the collection
  data_frame_temp_1 <- data.frame(LS=LS,dim=dimension,complexite=lchoose(p,dimension))      # characteristics of the model collection
  data_frame_temp_2 <- data.frame(beta=beta_estimator)
  path_collection_bounds_bolasso[[j]] <- list(val1 = data_frame_temp_1, val2 = data_frame_temp_2, val3 = order_index, val4 = frequence_vect, val5 = support_num)
  setTxtProgressBar(pb,j)
}

## Slope
path_collection_bounds_slope <- list()
Y_bounds <- Y_bounds_matrix[,l]
pb = txtProgressBar(min = 1, max = length(seq(1,n*nbr_it_bound,by=n)), initial = 1)
for (k in seq(1,n*nbr_it_bound,by=n)){
  j = which(seq(1,n*nbr_it_bound,by=n) ==k)
  set.seed(j+1234)
  index_resampling <- lapply(1:m, function(x) x <- unique(sort(sample(1:n, n, replace = TRUE))))     # creation of the resamples
  sampling_regressor <- lapply(index_resampling, function(x) sampling <- data_select[x,])
  sampling_regressed <- lapply(index_resampling, function(x) sampling <- Y_bounds[k:(k+(n-1))][x])
  slope_results <- lapply(1:m, function(i) SLOPE(x = sampling_regressor[[i]], y = sampling_regressed[[i]], family = "gaussian", intercept = TRUE))
  matrix_freq <- matrix(0, nrow = m, ncol = p)
  for (l in 1:m){
    matrix_temp <- matrix(0, nrow = (length(slope_results[[l]]$nonzeros==TRUE)/p), ncol = p)
    for (x in seq(1,length(slope_results[[l]]$nonzeros==TRUE),by=p)){
      z = which(seq(1,length(slope_results[[l]]$nonzeros==TRUE),by=p) ==x)
      matrix_temp[z,] <- 1*slope_results[[l]]$nonzeros[x:(x+(p-1))]
    }
    matrix_freq[l,] <- apply(matrix_temp, 2, sum)/nrow(matrix_temp)     # variable frequencies on each regularization path
  }
  frequence_vect <- apply(matrix_freq,2,sum)/m      # total variable frequencies averaged on the resamples and the regularization paths
  order_index <- sort(frequence_vect,index.return= TRUE,decreasing = TRUE)$ix
  support_num <- list()    # list for the supports of variable index in the collection
  for(x in 2:(p+1)){
    support_num[[x]] <- order_index[1:(x-1)]    # construction of the model support
  }
  index <- 2:(p+1)
  dimension <- unlist(lapply(support_num, function(x) length(x)))
  re_estimation <- list()
  for (i in index){
    re_estimation[[i]] <- lm(Y_bounds[k:(k+(n-1))] ~ data_select[,support_num[[i]]] -1)  # No intercept, mean square estimators onto each model
    # The beta coefficients are re-estimated onto each model of the collection by the MSE minization
  }
  beta_new <- sapply(re_estimation, function(a) as.vector(a$coefficients))
  beta_estimator <- matrix(0, nrow = length(support_num), ncol = p)
  for (i in 2:length(support_num)){
    beta_estimator[i,support_num[[i]]] <- as.numeric(beta_new[[i]])
  }
  LS = sapply(1:length(dimension), function(i)
    (1/nrow(data_select))*sum((Y_bounds[k:(k+(n-1))] - data_select %*% beta_estimator[i,])^2))   # The least squared values are computed for each beta of the collection
  data_frame_temp_1 <- data.frame(LS=LS,dim=dimension,complexite=lchoose(p,dimension))      # characteristics of the model collection
  data_frame_temp_2 <- data.frame(beta=beta_estimator)
  path_collection_bounds_slope[[j]] <- list(val1 = data_frame_temp_1, val2 = data_frame_temp_2, val3 = order_index, val4 = frequence_vect, val5 = support_num)
  setTxtProgressBar(pb,j)
}

## Knockoffs
path_collection_bounds_knockoffs <- list()
Y_bounds <- Y_bounds_matrix[,l]
pb = txtProgressBar(min = 1, max = length(seq(1,n*nbr_it_bound,by=n)), initial = 1)
for (k in seq(1,n*nbr_it_bound,by=n)){
  j = which(seq(1,n*nbr_it_bound,by=n) ==k)
  set.seed(j+1234)
  knockoffs_result <- knockoff.filter(X = data_select, y = Y_bounds[k:(k+(n-1))], knockoffs = create.second_order,statistic = stat.lasso_lambdasmax, fdr = alpha_FDR)
  order_index <- sort(knockoffs_result$statistic, index.return = TRUE, decreasing = TRUE)$ix
  support_num <- list()    # list for the supports of variable index in the collection
  for(x in 2:(p+1)){
    support_num[[x]] <- order_index[1:(x-1)]     # construction of the model support
  }
  index <- 2:(p+1)
  dimension <- unlist(lapply(support_num, function(x) length(x)))
  re_estimation <- list()
  for (i in index){
    re_estimation[[i]] <- lm(Y_bounds[k:(k+(n-1))] ~ data_select[,support_num[[i]]] -1)  # No intercept, mean square estimators onto each model
    # The beta coefficients are re-estimated onto each model of the collection by the MSE minization
  }
  beta_new <- sapply(re_estimation, function(a) as.vector(a$coefficients))
  beta_estimator <- matrix(0, nrow = length(support_num), ncol = p)
  for (i in 2:length(support_num)){
    beta_estimator[i,support_num[[i]]] <- as.numeric(beta_new[[i]])
  }
  LS = sapply(1:length(dimension), function(i)
    (1/nrow(data_select))*sum((Y_bounds[k:(k+(n-1))] - data_select %*% beta_estimator[i,])^2)) # The least squared values are computed for each beta of the collection
  data_frame_temp_1 <- data.frame(LS=LS,dim=dimension,complexite=lchoose(p,dimension))      # characteristics of the model collection
  data_frame_temp_2 <- data.frame(beta=beta_estimator)
  path_collection_bounds_knockoffs[[j]] <- list(val1 = data_frame_temp_1, val2 = data_frame_temp_2, val3 = order_index, val5 = support_num)
  setTxtProgressBar(pb,j)
}

## Random forests
path_collection_bounds_RF <- list()
Y_bounds <- Y_bounds_matrix[,l]
pb = txtProgressBar(min = 1, max = length(seq(1,n*nbr_it_bound,by=n)), initial = 1)
for (k in seq(1,n*nbr_it_bound,by=n)){
  j = which(seq(1,n*nbr_it_bound,by=n) ==k)
  set.seed(j+1234)
  RF_results <- rfe(x = data_select, y = Y_bounds[k:(k+(n-1))], sizes = p, rfeControl = rfeControl(functions = rfFuncs, number = m, method = "boot", saveDetails = TRUE), rerank = TRUE)   # RandomForests
  # the sizes parameter is the number of features that should be retained.
  # rfFuncs is the function for random forests
  # the rerank = TRUE oarameter is for recomputing the order on variables after elimination of least important variables (backward strategy)
  # the method = boot paramater is for the bootstrap strategy
  # the number parameter is the number of resamples
  order_index <- as.numeric(sapply(RF_results$optVariables, function(x) which(x==colnames(data_select))))
  support_num <- list()    # list for the supports of variable index in the collection
  for(x in 2:(p+1)){
    support_num[[x]] <- order_index[1:(x-1)]    # construction of the model support
  }
  dimension <- unlist(lapply(support_num, function(x) length(x)))
  index <- 2:(p+1)
  re_estimation <- list()
  for (i in index){
    re_estimation[[i]] <- lm(Y_bounds[k:(k+(n-1))] ~ data_select[,support_num[[i]]] -1)  # No intercept, mean square estimators onto each model
    # The beta coefficients are re-estimated onto each model of the collection by the MSE minization
  }
  beta_new <- sapply(re_estimation, function(a) as.vector(a$coefficients))
  beta_estimator <- matrix(0, nrow = length(support_num), ncol = p)
  for (i in 2:length(support_num)){
    beta_estimator[i,support_num[[i]]] <- as.numeric(beta_new[[i]])
  }
  LS = sapply(1:length(dimension), function(i)
    (1/nrow(data_select))*sum((Y_bounds[k:(k+(n-1))] - data_select %*% beta_estimator[i,])^2))  # The least squared values are computed for each beta of the collection
  data_frame_temp_1 <- data.frame(LS=LS,dim=dimension,complexite=lchoose(p,dimension))      # characteristics of the model collection
  data_frame_temp_2 <- data.frame(beta=beta_estimator)
  path_collection_bounds_RF[[j]] <- list(val1 = data_frame_temp_1, val2 = data_frame_temp_2, val3 = order_index, val5 = support_num)
  setTxtProgressBar(pb,j)
}

path_collection_random <- list(path_collection_bolasso,path_collection_slope,path_collection_RF,path_collection_knockoffs)
save(path_collection_random,file=paste("path_collection_random",config_data,sep="_"))
path_collection_bounds_random <- list(path_collection_bounds_bolasso,path_collection_bounds_slope,path_collection_bounds_RF,path_collection_bounds_knockoffs)
save(path_collection_bounds_random,file=paste("path_collection_bounds_random",config_data,sep="_"))

# # the configuration 2 of the scenario (ii)  of the article
# path_collection_random <- list(path_collection_bolasso,path_collection_slope,path_collection_RF,path_collection_knockoffs)    
# save(path_collection_random,file=paste("path_collection_smaller_noise_random",config_data,sep="_"))
# path_collection_bounds_random <- list(path_collection_bounds_bolasso,path_collection_bounds_slope,path_collection_bounds_RF,path_collection_bounds_knockoffs)
# save(path_collection_bounds_random,file=paste("path_collection_bounds_smaller_noise_random",config_data,sep="_"))

# # the configuration 3 of the scenario (ii)  of the article
# path_collection_random <- list(path_collection_bolasso,path_collection_slope,path_collection_RF,path_collection_knockoffs)    
# save(path_collection_random,file=paste("path_collection_larger_noise_closed_values_random",config_data,sep="_"))
# path_collection_bounds_random <- list(path_collection_bounds_bolasso,path_collection_bounds_slope,path_collection_bounds_RF,path_collection_bounds_knockoffs)
# save(path_collection_bounds_random,file=paste("path_collection_bounds_larger_noise_closed_values_random",config_data,sep="_"))
#####################################################################################################################





#####################################################################################################################
###### Comparaison of the generated model collections
#####################################################################################################################
# To compare the model collections, we calculate the proportion of active variables in models of size 5, 10, 15 and 20 of each of model collections.
# It allows to quantify the ability to discriminate between active and non active variables

load(paste("path_collection_deterministic",config_data,sep="_"))
load(paste("path_collection_random",config_data,sep="_"))
# load(paste("path_collection_smaller_noise_deterministic",config_data,sep="_"))
# load(paste("path_collection_smaller_noise_random",config_data,sep="_"))
# load(paste("path_collection_larger_noise_closed_values_deterministic",config_data,sep="_"))
# load(paste("path_collection_larger_noise_closed_values_random",config_data,sep="_"))

matrix_TP_P <- matrix(0, nrow = 4, ncol = length(path_collection)+length(path_collection_random))
colnames(matrix_TP_P) <- c("permutation_1_10","exchanging9_10_11_12","exchanging1_2_16_20","permutation_1_10","permutation_1_12","permutation_1_15","bolasso","slope","RandomForests","Knockoff")
rownames(matrix_TP_P) <- c("model_size_5","model_size_10","model_size_15","model_size_20")
index_temp <- c(5,10,15,20)
TP_possible <- c(5,10,10,10)
for (l in 1:length(index_temp)){
  for (y in 1:length(path_collection)){
    matrix_TP_P[l,y] <- sum(sapply(1:nbr_it, function(it) sum(1*sapply(which(path_collection[[y]][[it]]$val2[(index_temp[l]+1),]!=0), function(x) is.element(x,1:10)))/TP_possible[l]))/nbr_it
  }
  for (y in 1:length(path_collection_random)){
    matrix_TP_P[l,(y+length(path_collection))] <- sum(sapply(1:nbr_it, function(it) sum(1*sapply(which(path_collection_random[[y]][[it]]$val5[[(index_temp[l]+1)]]!=0), function(x) is.element(x,1:10)))/TP_possible[l]))/nbr_it
  }
}
xtable(matrix_TP_P, type = "latex", file = "matrice.tex",digits=3)
