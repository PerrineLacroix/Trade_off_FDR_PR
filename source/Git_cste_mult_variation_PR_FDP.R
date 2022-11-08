# This function source computes the least-squared, the predictive risk and the false discovery proportion values for each model of the given collections. 
# The least-squared, the predictive risk and the false discovery proportion values of the selected models from each value of considered multiplicative constant K are also extracted, 
# both with known and unknown variance. 
# In this last case, the variance is estimated by the slope heuristics algorithm implemented in the R package Capushe. 
cste_mult_variation_LS_PR_FDR <- function(cste_mult_vect,nbr_it,path,data_select,data_metric,Y_metric,size_beta_true,n,p,sigma_2){
  LS_path <- list()   # least-squared values along all the paths
  PR_path <- list()   # predictive risk values along all the paths
  FDP_path <- list()   # false discovery proportion along all the paths
  PR_hat <- matrix(nrow=nbr_it, ncol=length(cste_mult_vect))   # predictive risk of the selected model for each path and each K
  FDP_hat <- matrix(nrow=nbr_it, ncol=length(cste_mult_vect))   # false discovery proportion of the selected model for each path and each K
  model_hat <- matrix(nrow=nbr_it, ncol=length(cste_mult_vect))   # selected model for each path and each K
  ########## Computation of PR, FDP and the selected model for each nbr_it iterations (for empirical estimations)
  for (k in seq(1,n*nbr_it,by=n)){
    j = which(seq(1,n*nbr_it,by=n) ==k)
    LS_path[[j]] <- path[[j]]$val1$LS
    PR_path[[j]] <- sapply(1:length(path[[j]]$val1$dim), function(i) (1/nrow(data_metric))*sum((Y_metric[k:(k+(n-1))] - data_metric %*% as.numeric(path[[j]]$val2[i,]))^2))
            # predictive risk computation by using the validation data sets
    FDP_path[[j]] <- rep(0,length(path[[j]]$val1$dim))
    FDP_path[[j]][which(path[[j]]$val1$dim > size_beta_true)] <- sapply(which(path[[j]]$val1$dim > size_beta_true), function(i) (length(which(path[[j]]$val2[i,] !=0)) - size_beta_true)/length(which(path[[j]]$val2[i,] !=0)))
            # as soon as D_m > D_m^*, the number of FP in the model is (D_m - D_m^*) and the number of P is D_m
    model_temp <- c()
    for (i in cste_mult_vect){
      model_temp <- c(model_temp,which.min(i*sigma_2*(path[[j]]$val1$dim/n) + LS_path[[j]]))    
            # Model selected with respect to the multiplicative constante K
    }
    PR_hat[j,] <- as.numeric(sapply(1:length(cste_mult_vect), function(i) PR_path[[j]][model_temp[i]]))  
    FDP_hat[j,] <- as.numeric(sapply(1:length(cste_mult_vect), function(i) FDP_path[[j]][model_temp[i]]))   
    model_hat[j,] <- model_temp
  }
  ########## Computation of PR, FDP and the selected model for only one iteration (for data-driven evaluations)
  #### with known variance
  model_temp <- c()
  for (i in cste_mult_vect){
    model_temp <- c(model_temp,which.min(i*sigma_2*(path[[1]]$val1$dim/n) + path[[1]]$val1$LS))    
            # Model selected with respect to the multiplicative constante K
  }
  model_2 <- which.min(2*sigma_2*(path[[1]]$val1$dim/n) + path[[1]]$val1$LS)       
            # Model selected for the multiplicative constant 2
  LS_one_dataset <- as.numeric(sapply(1:length(cste_mult_vect), function(i) path[[1]]$val1$LS[model_temp[i]]))    
              # least-sqared values of the selected model for each K
  PR_one_dataset <- as.numeric(sapply(1:length(cste_mult_vect), function(i) PR_path[[1]][model_temp[i]]))    
              # predictive risk of the selected model for each K      
  diff_PR_one_dataset <- as.numeric(sapply(1:length(cste_mult_vect), function(i) (1/nrow(data_select))*sum((data_select %*% as.numeric(path[[1]]$val2[model_2,]) - data_select %*% as.numeric(path[[1]]$val2[model_temp[i],]))^2)))   
              # difference between the prediction risk obtained with K=2 and thoses obtained with the others multiplicative constant K
  #### with unknown variance and by using the slope heuristics principle in 1D
  forcap <- cbind(path[[1]]$val1$dim,path[[1]]$val1$dim/n,path[[1]]$val1$dim,path[[1]]$val1$LS)
  ResCapushe <- capushe(data = forcap)      # Application of the slope heuristics principle in 1D.
  model_temp_slope <- c()
  for (i in cste_mult_vect){
    model_temp_slope <- c(model_temp_slope,which.min(i*as.numeric(ResCapushe@DDSE@graph$reg$coefficients[2])*(path[[1]]$val1$dim/n) + path[[1]]$val1$LS))   
              # Model selected with respect to the multiplicative constante
  }
  model_temp_2 <- which.min(2*as.numeric(ResCapushe@DDSE@graph$reg$coefficients[2])*(path[[1]]$val1$dim/n) + path[[1]]$val1$LS)       # Model selected for the constant 2
              # Model selected for the multiplicative constant 2
  LS_one_dataset_slope <- as.numeric(sapply(1:length(cste_mult_vect), function(i) path[[1]]$val1$LS[model_temp_slope[i]]))
              # least-sqared values of the selected model for each K
  PR_one_dataset_slope <- as.numeric(sapply(1:length(cste_mult_vect), function(i) PR_path[[1]][model_temp_slope[i]]))
              # predictive risk of the selected model for each K      
  diff_PR_one_dataset_slope <- as.numeric(sapply(1:length(cste_mult_vect), function(i) (1/nrow(data_select))*sum((data_select %*% as.numeric(path[[1]]$val2[model_temp_2,]) - data_select %*% as.numeric(path[[1]]$val2[model_temp_slope[i],]))^2)))   
              # difference between the prediction risk obtained with K=2 and thoses obtained with the others multiplicative constant K
  res1 <- data.frame(PR_path = I(PR_path), FDP_path = I(FDP_path))      
              # characteristics of the paths
  res2 <- data.frame(model_hat = I(model_hat), PR_hat = I(PR_hat), FDP_hat = I(FDP_hat))     
              # characteristics of the selected models for nbr_it iterations
  res3 <- data.frame(model_one_dataset = I(model_temp), LS_one_dataset = I(LS_one_dataset), 
                     PR_one_dataset = I(PR_one_dataset), diff_PR_one_dataset = I(diff_PR_one_dataset))     
              # characteristics of the selected model for one iteration with a known variance 
  res4 <- data.frame(model_one_dataset_slope = I(model_temp_slope), LS_one_dataset_slope = I(LS_one_dataset_slope), 
                     PR_one_dataset_slope = I(PR_one_dataset_slope), diff_PR_one_dataset_slope = I(diff_PR_one_dataset_slope))     
              # characteristics of the selected model for one iteration with an unknown variance
  res5 <- data.frame(model_2 = I(model_2),model_temp_2 = I(model_temp_2),kappa_slope = I(as.numeric(ResCapushe@DDSE@graph$reg$coefficients[2])))       
              # characteristics of the applied slope heuristics
  return(list(path = res1, empirical = res2, known_variance = res3, unknown_variance = res4, parameter = res5))
}