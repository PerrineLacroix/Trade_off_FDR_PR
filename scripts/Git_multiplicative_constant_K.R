## script to launch simulations
## author : perrine lacroix
## date : April, 10 2024

## This code is the main R script. 
## From the created data sets and the generated model collections (from Git_multiplicative_constant_K_creation_data_and_collection.R and Git_multiplicative_constant_K_collection_non_ordered.R R functions),
## the PR and the FDR are firstly computed for each K and for each iteration. 
## Then, the P_{2r} term in the FDR expression is estimated by an empirical process.
## Lower bounds and upper bounds are computed, for each K and for all the \Tilde{K} considered in the article.
## The rest of the R script stands for plottings and for comparaison between our algorithm and some existing variable selection procedures.


rm(list=objects())

prefix = "~/Documents/GitHub/Trade_off_FDR_PR"

# Packages loading
library(xtable)
library(capushe)
library(LINselect)
library(knockoff)
library(caret)

## Loading of coded function
setwd(paste(prefix,"/source",sep =""))
source("Git_cste_mult_variation_PR_FDP.R")
source("Git_empirical_estimation.R")
source("Git_P_2r_estimation.R")
source("Git_lower_bounds.R")
source("Git_upper_bounds.R")
source("Git_lower_bound_estimated_one_dataset.R")
source("Git_upper_bound_estimated_one_dataset.R")
source("Git_empirical_estimation_collection_non_ordered.R")
source("Git_algorithm_trade_off.R")
source("Git_LinSelect.R")

## Loading of the dataset
setwd(paste(prefix,"/data",sep =""))
# Parameters of simulation
nbr_it = 1000   # Data sets number for the empirical estimations of FDR and PR. 
nbr_it_bound = 100   # Data sets number to compute the FDR bounds and to study the variability of the FDR bounds.
n = 50  # number of observations
p = 50   # number of variables
sigma_2 <- 1  # Variance of the noise in the Gaussian linear regression
cste_mult_vect <- seq(0.1,10,length.out = 20)    # The considered K constants
cste_mult_vect <- c(cste_mult_vect[which(cste_mult_vect<2)],2,cste_mult_vect[which(cste_mult_vect>2)])    # To add the constant 2 in the studied K
asympt_number <- 5000    # Number for the empirical estimation of the P_{2r} term in the FDR expression. 
## for the existing variable selection procedure: 
K_fold = 50      # for K-fold CV
max.steps_LIN = NULL      # by default for LinSelect
dmax_LIN = NULL      # by default for LinSelect
coeff_mult_pen_LIN = 1.1     # > 1      # by default for LinSelect
coeff_mult_proj_LIN = 0.5     # > 0      # by default for LinSelect
alpha_FDR <- 0.05     # for Knockoffs and our algorithm
## for the chosen model collections
ordered_variable <- TRUE    # TRUE for the deterministic ordered model collections; FALSE for the deterministic non ordered model collections and for random collections 
random_collection <- FALSE    # TRUE for random collections; FALSE for the deterministic ordered model collections and for the deterministic non ordered model collections

# Loading the generated data sets
config_data = paste0(paste("sigma2",sigma_2,"p",p,"n",n,"nbr_it",nbr_it,sep="_"),".RData")
load(paste("matrix_X",config_data,sep="_"))
config_data = paste0(paste("beta_true_0_20","sigma2",sigma_2,"p",p,"n",n,"nbr_it",nbr_it,sep="_"),".RData")
load(paste("vector_Y",config_data,sep="_"))
# load(paste("smaller_noise_vector_Y",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("larger_noise_closed_values_vector_Y",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article
load(paste("vector_Y_noise_bounds",config_data,sep="_"))
# load(paste("vector_Y_smaller_noise_bounds",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("vector_Y_larger_noise_closed_values_bounds",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article
if (ordered_variable == TRUE){
  index_plot <- 1:ncol(vrai_beta_matrix)
}else{
  if (random_collection == FALSE){
    l = 11       # the dimension of the true model is fixed to 10
    index_plot <- 1:6    # 6 deterministic non ordered model collections
    config_data = paste0(paste("deterministic_beta_true_0_20","sigma2",sigma_2,"p",p,"n",n,"nbr_it",nbr_it,sep="_"),".RData")
  }else{
    l = 11       # the dimension of the true model is fixed to 10
    index_plot <- 1:4      # 4 random model collections
    config_data = paste0(paste("random_beta_true_0_20","sigma2",sigma_2,"p",p,"n",n,"nbr_it",nbr_it,sep="_"),".RData")
  }
}
load(paste("path_collection",config_data,sep="_"))
# load(paste("path_collection_smaller_noise",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("path_collection_larger_noise_closed_values",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article
load(paste("path_collection_bounds",config_data,sep="_"))
# load(paste("path_collection_bounds_smaller_noise",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("path_collection_bounds_larger_noise_closed_values",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article
if (random_collection == TRUE){
  path_collection <- path_collection_random
  path_collection_bounds <- path_collection_bounds_random
}


#####################################################################################################################
###### Computation of the PR and FDP of the selected model for each K and for each iteration
#####################################################################################################################
#### The following lines give the PR and the FDP of the selected model for each considered K and on all the nbr_it data sets
#### PR_FDP_Hat_K is a matrix of size nbr_it x length(cste_mult_vect)
if (ordered_variable == TRUE){
  PR_FDP_hat_K <- empirical_estimation(cste_mult_vect,nbr_it,n,p,sigma_2,path_collection,Y_metric_matrix,data_select,data_metric)
}else{
  PR_FDP_hat_K <- empirical_estimation_non_ordered(cste_mult_vect,nbr_it,n,p,sigma_2,path_collection,Y_metric_matrix,data_select,data_metric)
}
setwd(paste(prefix,"/data",sep =""))
save(PR_FDP_hat_K, file = paste("PR_FDP_hat_K",config_data,sep="_"))
# save(PR_FDP_hat_K, file = paste("smaller_noise_PR_FDP_hat_K",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# save(PR_FDP_hat_K, file = paste("larger_noise_closed_values_PR_FDP_hat_K",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article





#####################################################################################################################
###### Empirical estimation of the P_{2r} term in the FDR expression
#####################################################################################################################
setwd(paste(prefix,"/data",sep =""))
load(paste("PR_FDP_hat_K",config_data,sep="_"))
# load(paste("smaller_noise_PR_FDP_hat_K",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("larger_noise_closed_values_PR_FDP_hat_K",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article

# P_{2r} is deterministic if r is known. r varies between D_m_star+1 and q.
# We generate P_{2r} for r in [D_m_star +1,q] for all possible D_m_star (in [0,q])
# Then, P_{2r} will be used for appropriated r (according to D_m_star or D_m_hat)
list_P_2r <- P_2r_estimation(n,p,asympt_number,cste_mult_vect)
setwd(paste(prefix,"/data",sep =""))
save(list_P_2r,file=paste0(paste("estimation_of_P_2r_for_all_D_m_in_0_q","p",p,"n",n,"asympt_number",asympt_number,sep="_"),".RData"))




#####################################################################################################################
###### LOWER BOUNDS and UPPER BOUNDS
#####################################################################################################################

############# LOWER BOUNDS
setwd(paste(prefix,"/data",sep =""))
load(paste0(paste("estimation_of_P_2r_for_all_D_m_in_0_q","p",p,"n",n,"asympt_number",asympt_number,sep="_"),".RData"))
detail <- FALSE
if (ordered_variable == FALSE){
  l = 11     # the dimension of the true model is fixed to 10
  true_beta <- matrix(NA,nrow = nrow(vrai_beta_matrix), ncol = length(path_collection))
  for (i in 1:length(path_collection)){
    true_beta[,i] <- vrai_beta_matrix[,l]
  }
  vrai_beta_matrix <- true_beta
}
list_lower_bounds <- lower_bounds(detail,vrai_beta_matrix,n,p,cste_mult_vect,sigma_2,path_collection,path_collection_bounds,list_P_2r)
save(list_lower_bounds,file=paste("lower_bounds_all",config_data,sep="_"))
# save(list_lower_bounds,file=paste("smaller_noise_lower_bounds_all",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# save(list_lower_bounds,file=paste("larger_noise_closed_values_lower_bounds_all",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article

############# UPPER BOUNDS
load(paste0(paste("estimation_of_P_2r_for_all_D_m_in_0_q","p",p,"n",n,"asympt_number",asympt_number,sep="_"),".RData"))
load(paste("lower_bounds_all",config_data,sep="_"))
# load(paste("smaller_noise_lower_bounds_all",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("larger_noise_closed_values_lower_bounds_all",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article
detail <- FALSE
comparison <- TRUE
with_hat <- TRUE
list_upper_bounds <- upper_bounds(detail,vrai_beta_matrix,n,p,cste_mult_vect,sigma_2,path_collection,path_collection_bounds,list_P_2r)
save(list_upper_bounds,file=paste("upper_bounds_all",config_data,sep="_"))
# save(list_upper_bounds,file=paste("smaller_noise_upper_bounds_all",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# save(list_upper_bounds,file=paste("larger_noise_closed_values_upper_bounds_all",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article





#####################################################################################################################
###### Plotting
#####################################################################################################################
setwd(paste(prefix,"/data",sep =""))
load(paste("PR_FDP_hat_K",config_data,sep="_"))
load(paste0(paste("estimation_of_P_2r_for_all_D_m_in_0_q","p",p,"n",n,"asympt_number",asympt_number,sep="_"),".RData"))
load(paste("lower_bounds_all",config_data,sep="_"))
load(paste("upper_bounds_all",config_data,sep="_"))
# load(paste("smaller_noise_lower_bounds_all",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("smaller_noise_upper_bounds_all",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("larger_noise_closed_values_lower_bounds_all",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article
# load(paste("larger_noise_closed_values_upper_bounds_all",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article

description_parameters = paste0(paste("sigma2",sigma_2,"p",p,"n",n,sep="_"),".RData")
pente_parameters = paste("final_upper_and_lower_bound_of_FDR",description_parameters,sep="_")
# pente_parameters = paste("final_upper_and_lower_bound_of_FDR_smaller_noise",description_parameters,sep="_")    # the configuration 2 of the scenario (ii)  of the article
# pente_parameters = paste("final_upper_and_lower_bound_of_FDR_larger_noise",description_parameters,sep="_")    # the configuration 3 of the scenario (ii)  of the article
setwd(paste(prefix,"/results",sep =""))
pdf(paste0(pente_parameters,".pdf"))

###### Plot of the FDR, theoretical upper bound and theoretical lower bound functions with respect to K
for (k in index_plot){
  quantile_asympt_95 <- 1.96      # for the confidence interval of the FDR empirical estimation. Application of the central limit theorem
  vect_empirical_FDR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))   # empirical FDR
  variance_estimation_FDR <- sapply(1:length(cste_mult_vect), function(it) (1/(nbr_it-1))*sum((PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]-mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))^2))    # empirical variance
  FDP_confidence_interval_inf <- vect_empirical_FDR - (variance_estimation_FDR*quantile_asympt_95)/sqrt(nbr_it)     # lower bound of the confidence interval for the FDR empirical estimation
  FDP_confidence_interval_sup <- vect_empirical_FDR + (variance_estimation_FDR*quantile_asympt_95)/sqrt(nbr_it)      # upper bound of the confidence interval for the FDR empirical estimation
  par(mfrow=c(1,1))
  plot(cste_mult_vect,FDP_confidence_interval_inf, col = 1,type = 'l',xlim = c(0.1,10),
       ylim=c(0,max(list_upper_bounds$upper_bound_list[[k]])),xlab = "K", ylab = NA, lty = 5, lwd=2, cex.lab = 1.3)
  lines(cste_mult_vect,FDP_confidence_interval_sup, col = 1, type = 'l', lty = 5, lwd=2)
  lines(cste_mult_vect,vect_empirical_FDR,col = 2, type="p", lty = 5, lwd=2)
  lines(cste_mult_vect,vect_empirical_FDR,col = 2,type="l", lwd=2)
  lines(cste_mult_vect,list_upper_bounds$upper_bound_list[[k]], col = 4, type="p", lwd=2)
  lines(cste_mult_vect,list_upper_bounds$upper_bound_list[[k]], col = 4, type="l", lwd=2)
  lines(cste_mult_vect,list_lower_bounds$lower_bound_list[[k]], col = 3, type = "p", lwd=2)
  lines(cste_mult_vect,list_lower_bounds$lower_bound_list[[k]], col = 3, type = "l", lwd=2)
  abline(h=0,col="black",lty=2)
  legend("topright", legend = c("empirical FDR","B","b"),
                                #"95% empirical confidence interval"
         col = c(2,4,3,1), lty = 1, cex=1.1, lwd = c(4,4,4,4,4))
}

# For only constants K >= 2
which_cste_sup2 <- which(cste_mult_vect>=2)
for (k in index_plot){
  quantile_asympt_95 <- 1.96      # for the confidence interval of the FDR empirical estimation. Application of the central limit theorem
  vect_empirical_FDR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))   # empirical FDR
  variance_estimation_FDR <- sapply(1:length(cste_mult_vect), function(it) (1/(nbr_it-1))*sum((PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]-mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))^2))    # empirical variance
  FDP_confidence_interval_inf <- vect_empirical_FDR - (variance_estimation_FDR*quantile_asympt_95)/sqrt(nbr_it)     # lower bound of the confidence interval for the FDR empirical estimation
  FDP_confidence_interval_sup <- vect_empirical_FDR + (variance_estimation_FDR*quantile_asympt_95)/sqrt(nbr_it)      # upper bound of the confidence interval for the FDR empirical estimation
  par(mfrow=c(1,1))
  plot(cste_mult_vect[which_cste_sup2],FDP_confidence_interval_inf[which_cste_sup2], col = 1,type = 'l',
       ylim = c(0,max(list_upper_bounds$upper_bound_list[[k]][which_cste_sup2])),
       xlab = "K",ylab = NA,lty = 5, lwd=2, cex.lab = 1.3)
  lines(cste_mult_vect[which_cste_sup2],FDP_confidence_interval_sup[which_cste_sup2], col = 1, type = 'l', lty = 5, lwd=2)
  lines(cste_mult_vect[which_cste_sup2],vect_empirical_FDR[which_cste_sup2],col = 2, type="p", lty = 5, lwd=2)
  lines(cste_mult_vect[which_cste_sup2],vect_empirical_FDR[which_cste_sup2],col = 2,type="l", lwd=2)
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="p", lwd=2)
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="l", lwd=2)
  lines(cste_mult_vect,list_lower_bounds$lower_bound_list[[k]], col = 3, type = "p", lwd=2)
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 3, type = "l", lwd=2)
  abline(h=0,col="black",lty=2)
  legend("topright", legend = c("empirical FDR","B","b"),
         #"95% empirical confidence interval"
         col = c(2,4,3,1), lty = 1, cex=1.1, lwd = c(4,4,4,4,4))
}

##### Plot of the Empirical FDR and PR values with respect to K
for (k in index_plot){
  vect_empirical_FDR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))   # empirical FDR
  vect_empirical_PR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$PR_hat[,it]))    # empirical PR
  par(mfrow=c(1,2))
  plot(cste_mult_vect,vect_empirical_FDR, col = 2,xlab = "K",ylim=c(0,1),ylab = NA,type= "p", lty = 5, cex.lab = 1.3, lwd = 3)
  lines(cste_mult_vect,vect_empirical_FDR, col = 2,type="l", lwd = 3)
  legend("topright", legend = c("empirical FDR"), col = c(2), lty = 1, cex=1.1, lwd = c(4))
  abline(v=3,lwd=2)
  plot(cste_mult_vect,vect_empirical_PR, col = 4, type= "p", lty = 5,xlab = "K",ylab = NA,cex.lab = 1.3, lwd = 3)
  lines(cste_mult_vect,vect_empirical_PR, col = 4,type="l", lwd = 3)
  legend("topright", legend = c("empirical PR"), col = c(4), lty = 1, cex=1, lwd = c(4))
  abline(v=3,lwd=2)
}

##### Values of the Empirical FDR and PR values with respect to K
for (k in index_plot){
  if (ordered_variable == TRUE){
    matrix_selection_1 <- matrix(0, nrow = 2, ncol = 11)
    colnames(matrix_selection_1) <- c("0.1","1","2","3","4","5","6","7","8","9","10")
    rownames(matrix_selection_1) <- c("empirical FDR(m_Hat)","empirical PR(m_Hat)")
    matrix_selection_1[1,] <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))[c(1,3,5,8,10,12,14,16,18,19,21)]   # empirical FDR
    matrix_selection_1[2,] <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$PR_hat[,it]))[c(1,3,5,8,10,12,14,16,18,19,21)]
    xtable(matrix_selection_1, type = "latex", file = "matrice.tex",digits=2)
  }
}


###### Plot of the estimated FDR and PR (calculated on only one dataset) with respect to K
sequence_temp <- seq(1,length(list_lower_bounds$lower_bound_list_hat_one_dataset_K2_list_slope[[1]][[1]]), by=2)
for (k in index_plot){
  vect_estimation_diff_PR <- PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope    # estimated difference in predictions
  par(mfrow=c(1,2))
  plot(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]),xlab = "K",ylab = "estimated upper bound of FDR on one dataset",type="p", col = 6)
  lines(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]),type="l", col = 6)
  legend("topright", legend = c("Hat_upper_bound(beta(Hat(m)(4)))"), col = c(6), lty = 1, cex=0.5, lwd = c(4))
  plot(cste_mult_vect,vect_estimation_diff_PR, col = 4,xlab = "K",ylab = "estimated PR on one dataset", type="p")
  lines(cste_mult_vect,vect_estimation_diff_PR, col = 4,type="l")
  legend("topright", legend = c("K -> estimated PR(Hat(m)(K))"), col = c(4), lty = 1, cex=0.5, lwd = c(4))
}


########### Histogram of the \Hat{\kappa} values corresponding to estimation values of the \sigma^2 parameter by the slope heuristic with the capushe R package
for (k in index_plot){
  par(mfrow=c(1,1))
  vect_kappa <- sapply(1:nbr_it_bound, function(it) list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[it]]$kappa_slope)
  hist(vect_kappa, freq = NULL, main = NULL, ylab = "Frequency", xlab = "Hat(sigma^2)")
}


########### plot of the estimated FDR and PR (calculated on only one dataset) and the empirical FDR and PR with respect to K and with \Tilde{K} = 4
cste_for_estimation <- 4     # \Tilde{K} = 4
for (k in index_plot){
  vect_empirical_FDR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))     # empirical FDR
  vect_empirical_PR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$PR_hat[,it]))    # empirical PR
  par(mfrow=c(1,2))
  plot(cste_mult_vect,vect_empirical_FDR, col = 2,xlab = "K",ylim=c(0,max(as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]))),ylab = "empirical FDR and estimated upper bound of FDR ",type= "p", lty = 5)
  lines(cste_mult_vect,vect_empirical_FDR, col = 2,type="l")
  lines(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]),ylim=c(0,1),xlab = "K",ylab = "estimated upper bound of FDR on one dataset",type="p", col = 6)
  lines(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]),type="l", col = 6)
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","Hat_upper_bound(beta(Hat(m)(4)))"), col = c(2,6), lty = 1, cex=0.7, lwd = c(4,4))
  min_diff_PR <- min(vect_empirical_PR,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope)
  max_diff_PR <- max(vect_empirical_PR,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope)
  plot(cste_mult_vect,vect_empirical_PR, col = 4, type="p", lty = 5,ylim = c(min_diff_PR,max_diff_PR),xlab = "K",ylab = "empirical PR and estimated PR")
  lines(cste_mult_vect,vect_empirical_PR, col = 4,type="l")
  lines(cste_mult_vect,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope, col = "#660099", type="p")
  lines(cste_mult_vect,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope, col = "#660099", type="l")
  legend("topright", legend = c("K -> empirical PR(Hat(m)(K))","K -> estimated PR(Hat(m)(K))"), col = c(4,"#660099"), lty = 1, cex=0.7, lwd = c(4,4))
}

# For only constants K >= 2
which_cste_sup2 <- which(cste_mult_vect>=2)
sequence_temp <- seq(1,length(list_lower_bounds$lower_bound_list_hat_one_dataset_K2_list_slope[[1]][[1]]), by=2)
for (k in index_plot){
  vect_empirical_FDR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))     # empirical FDR
  vect_empirical_PR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$PR_hat[,it]))     # empirical PR
  vect_empirical_diff_PR <- vect_empirical_PR - vect_empirical_PR[which(cste_mult_vect==2)]       # empirical difference in predictions
  vect_estimation_diff_PR <- PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope         # estimated difference in predictions
  par(mfrow=c(1,2))
  plot(cste_mult_vect[which_cste_sup2],vect_empirical_FDR[which_cste_sup2], col = 2,xlab = "K",ylim=c(0,max(as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2])),
       ylab = NA,type= "p", lty = 5, cex.lab = 1.2, lwd = 2)
  lines(cste_mult_vect[which_cste_sup2],vect_empirical_FDR[which_cste_sup2], col = 2, type = 'l', lty = 1, lwd=2, pch = 2)
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],ylim=c(0,1),xlab = "K",ylab = "estimated upper bound of FDR on one dataset",type="p", col = 6, lty = 1, lwd=2, pch = 3)
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],type="l", col = 6, lty = 1, lwd=2, pch = 3)
  legend("topright", legend = c("empirical FDR","estimated", "upper bound", "with 4"), col = c(2,6), lty = 1, cex=1.3, lwd = c(4,4,0,0,0),pch = c(2,3))
  min_diff_PR <- min(vect_empirical_diff_PR,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope)
  max_diff_PR <- max(vect_estimation_diff_PR,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope)
  plot(cste_mult_vect[which_cste_sup2],vect_empirical_diff_PR[which_cste_sup2], col = 4, xlab = "K", ylab = NA,type="p",
       lty = 5,ylim = c(min_diff_PR,max_diff_PR), cex.lab = 1.3, lwd = 2, pch = 2)
  lines(cste_mult_vect[which_cste_sup2],vect_empirical_diff_PR[which_cste_sup2], col = 4,type="l", lty = 1, lwd=2)
  lines(cste_mult_vect[which_cste_sup2],PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope[which_cste_sup2], col = "#660099", type="p", lty = 1, lwd=2, pch = 3)
  lines(cste_mult_vect[which_cste_sup2],PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope[which_cste_sup2], col = "#660099", type="l", lty = 1, lwd=2)
  legend("topleft", legend = c("empirical","difference in","predictions","estimated","difference in", "predictions"),
         col = c(4,1,1,"#660099",1,1), lty = 1, cex=1.2, lwd = c(4,0,0,4,0,0), pch = c(2,3))
}


########### plot of the empirical selected model dimension with respect to K by using pen = K*\sigma^2*D_m
for (k in index_plot){
  par(mfrow=c(1,1))
  vect_dim <- apply(PR_FDP_hat_K[[k]]$empirical$model_hat,2,mean)
  plot(cste_mult_vect,vect_dim, col = "#660099",xlab = "K",ylab = NA,type= "p", lty = 5, cex.lab = 1.2, lwd = 2)
  lines(cste_mult_vect,vect_dim, col = "#660099", type = 'l', lty = 1, lwd=2, pch = 2)
  legend("topright", legend = c("selected","model","dimension"), col = c("#660099",1,1), lty = 1, cex=1.2, lwd = c(4,0,0), pch = c(1,2))
}


########### plot of the estimated difference in predictions (calculated on only one dataset) and the empirical difference in predictions with respect to K and with \Tilde{K} = 4
for (k in index_plot){
  vect_empirical_PR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$PR_hat[,it]))  # empirical PR
  vect_empirical_diff_PR <- vect_empirical_PR - vect_empirical_PR[which(cste_mult_vect==2)]        # empirical difference in predictions
  par(mfrow=c(1,1))
  plot(PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]),
       pch = 4,xlab = "difference in predictions",ylab = "FDR",
       ylim=c(0,max(as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]))),
       xlim=c(min(PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope,vect_empirical_diff_PR),max(PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope,vect_empirical_diff_PR)),
       type= "p", lty = 5, lwd=2, cex.lab = 1.3, col = c(rep(1,2),rep(2,2),rep(3,3),rep(4,2),rep(5,2),rep(6,2),rep(7,2),rep(8,2),rep(9,2),rep(10,2)))
  lines(vect_empirical_diff_PR,vect_empirical_FDR, type="p", pch = 2, lty = 5, lwd=2, cex.lab = 1.3, col = c(rep(1,2),rep(2,2),rep(3,3),rep(4,2),rep(5,2),rep(6,2),rep(7,2),rep(8,2),rep(9,2),rep(10,2)))
  legend("top", legend = c("empirical","theoretical","K: 0-1","K: 1-2","K: 2-3","K: 3-4","K: 4-5","K: 5-6","K: 6-7","K: 7-8","K: 8-9","K: 9-10"),
         col = c(4,2,c(1,2,3,4,5,6,7,8,9,10)),
         lty = c(0,0,rep(1,10)), cex=1.1, lwd = c(4,4,rep(1,length(cste_mult_vect))), pch = c(4,2,rep(20,10)))
}

dev.off()



#####################################################################################################################
###### Plottings for evaluations of the upper and lower bounds
#####################################################################################################################


###### FDR, theoretical lower bound and estimated lower bounds functions with respect to K
which_cste_sup2 <- which(cste_mult_vect>=2)
sequence_temp <- seq(1,length(list_lower_bounds$lower_bound_list_hat_one_dataset_K2_list_slope[[1]][[1]]), by=2)
for (k in index_plot){
  vect_empirical_FDR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))   # empirical FDR
  par(mfrow=c(1,1))
  colfunc<-colorRampPalette(c("red","yellow","springgreen","#196F3D"))
  colors <- c(heat.colors(300)[1:100],topo.colors(300)[1:100],terrain.colors(300)[1:100])
  plot(cste_mult_vect[which_cste_sup2],vect_empirical_FDR[which_cste_sup2],col = colors[1],type="p",xlim = c(2,10),xlab = "K",ylab = NA, lty = 5,lwd=4,cex.lab = 1.8,pch=1)
  lines(cste_mult_vect[which_cste_sup2],vect_empirical_FDR[which_cste_sup2],lwd=4,col = colors[1],type="l",pch=1)
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2],pch=3,lwd=4, col = colors[100], type="p")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2],pch=3,lwd=4, col = colors[100], type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K1_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=4, col = colors[200], type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K1_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=4, col = colors[200], type = "l")
  #lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=4, col = colors[141], type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K2_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=4, cex = 1, col = colors[211], type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K2_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=4, cex = 1, col = colors[211], type = "l")
  #lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=4, col = colors[141], type = "l")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K1.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2, col = colors[211], type = "p")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K1.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2, col = colors[211], type = "l")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K2.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2, cex = 1, col = colors[231], type = "p")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K2.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2, cex = 1, col = colors[231], type = "l")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K3_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2, cex = 1, col = colors[241], type = "p")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K3_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2, cex = 1, col = colors[241], type = "l")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K3.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2,cex = 1, col = colors[251], type = "p")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K3.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2,cex = 1, col = colors[251], type = "l")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_Klog_n_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2,cex = 1, col = colors[261], type = "p")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_Klog_n_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2,cex = 1, col = colors[261], type = "l")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2,cex = 1, col = colors[271], type = "p")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2,cex = 1, col = colors[271], type = "l")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K4.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2,cex = 1, col = colors[281], type = "p")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K4.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2,cex = 1, col = colors[281], type = "l")
  abline(h=0,col="black",lty=2)
  legend("topright", legend = c("empirical FDR","b","estimated b with K=1", "estimated b with K in [1.5,5]"),
         col = c(colors[1],colors[100],colors[200],colors[211]), cex=1.8, lwd = c(4,4,4,4,4),
         pch = c(1,3,4,4),lty = c(1,1,1,1,1))
  
  ###### Only for $\Tilde{K} = 4$
  plot(cste_mult_vect,vect_empirical_FDR,col = 2,type="p",xlim = c(0.1,10),ylim=c(0,1),xlab = "K",ylab = NA,lty = 5, lwd=4, cex.lab = 1.3, pch = 1,cex.lab = 1.8)
  lines(cste_mult_vect,vect_empirical_FDR,col = 2,type="l",lty = 1, lwd=4, cex.lab = 1.8, pch = 1)
  lines(cste_mult_vect,list_lower_bounds$lower_bound_list[[k]],pch=3,lwd=4, col = colors[100], type="p")
  lines(cste_mult_vect,list_lower_bounds$lower_bound_list[[k]],pch=3,lwd=4, col = colors[100], type="l")
  lines(cste_mult_vect,as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]), col = "#993300", type = "p",lty = 1, lwd=4, cex.lab = 1.3, pch = 4)
  lines(cste_mult_vect,as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]), col = "#993300", type = "l",lty = 1, lwd=4, cex.lab = 1.3, pch = 4)
  abline(h=0,col="black",lty=2)
  legend("topright", legend = c("empirical FDR","b","estimated b with K=4"),
         col = c(2,colors[100],"#993300"), lty = c(1,1,1), cex=1.8, lwd = c(4,4,4),pch = c(1,3,4))
}



###### FDR, theoretical upper bound and estimated upper bounds functions with respect to K
which_cste_sup2 <- which(cste_mult_vect>=2)
sequence_temp <- seq(1,length(list_upper_bounds$upper_bound_list_hat_one_dataset_K2_list_slope[[1]][[1]]), by=2)
for (k in index_plot){
  vect_empirical_FDR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))   # empirical FDR
  colfunc<-colorRampPalette(c("red","yellow","springgreen","#196F3D"))
  colors <- c(heat.colors(300)[1:100],topo.colors(300)[1:100],terrain.colors(300)[1:100])
  par(mfrow=c(1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_empirical_FDR[which_cste_sup2],pch=1,col = colors[1],type="p",xlim = c(2,10),xlab = "K",ylab = NA, lty = 5,lwd=4,cex.lab = 1.8,ylim = c(0,0.25))
  lines(cste_mult_vect[which_cste_sup2],vect_empirical_FDR[which_cste_sup2],pch=1,lwd=4,col = colors[1],type="l")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2],pch=3,lwd=4, col = colors[100], type="p")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2],pch=3,lwd=4, col = colors[100], type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K1_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=4, col = colors[200], type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K1_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=4, col = colors[200], type = "l")
  #lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=4,cex = 1, col = colors[141], type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K2_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=4, cex = 1, col = colors[211], type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K2_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=4, cex = 1, col = colors[211], type = "l")
  #lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=4,cex = 1, col = colors[141], type = "l")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K1.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2, col = colors[211], type = "p")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K1.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2, col = colors[211], type = "l")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K2.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2, cex = 1, col = colors[231], type = "p")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K2.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2, cex = 1, col = colors[231], type = "l")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K3_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2, cex = 1, col = colors[241], type = "p")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K3_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2, cex = 1, col = colors[241], type = "l")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K3.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2,cex = 1, col = colors[251], type = "p")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K3.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2,cex = 1, col = colors[251], type = "l")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_Klog_n_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2,cex = 1, col = colors[261], type = "p")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_Klog_n_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2,cex = 1, col = colors[261], type = "l")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2,cex = 1, col = colors[271], type = "p")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2,cex = 1, col = colors[271], type = "l")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2,cex = 1, col = colors[281], type = "p")
  # lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2],pch=4,lwd=2,cex = 1, col = colors[281], type = "l")
  abline(h=0,col="black",lty=2)
  legend("topright", legend = c("empirical FDR","B","estimated B with K=1", "estimated b with K in [1.5,5]"),
         col = c(colors[1],colors[100],colors[200],colors[211]), cex=1.8, lwd = c(4,4,4,4,4),
         pch = c(1,3,4,4),lty = c(1,1,1,1,1))
  
  ###### Only for $\Tilde{K} = 4$
  plot(cste_mult_vect,vect_empirical_FDR,col = 2,type="p",xlim = c(0.1,10),xlab = "K", ylab = NA,lty = 5, lwd=4, cex.lab = 1.8, pch =1, ylim = c(0,max(as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]))))
  lines(cste_mult_vect,vect_empirical_FDR,col = 2,type="l",lty = 1, lwd=4, cex.lab = 1.3, pch = 1)
  lines(cste_mult_vect,list_upper_bounds$upper_bound_list[[k]], col = colors[101], type="p",lty = 1, lwd=4, cex.lab = 1.3, pch = 3)
  lines(cste_mult_vect,list_upper_bounds$upper_bound_list[[k]], col = colors[101], type="l",lty = 1, lwd=4, cex.lab = 1.3, pch = 3)
  lines(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]), col = "#993300", type = "p",lty = 1, lwd=4, cex.lab = 1.3, pch = 4)
  lines(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]), col = "#993300", type = "l",lty = 1, lwd=4, cex.lab = 1.3, pch = 4)
  abline(h=0,col="black",lty=2)
  legend("topright", legend = c("empirical FDR","B","estimated B with K=4"),
         col = c(2,colors[101],"#993300"), lty = c(1,1,1), cex=1.8, lwd = c(4,4,4),pch = c(1,3,4))
}

###### Relative change of the estimated lower bounds functions with respect to K
which_cste_sup2 <- which(cste_mult_vect>=2)
sequence_temp <- seq(1,length(list_lower_bounds$lower_bound_list_hat_one_dataset_K2_list_slope[[1]][[1]]), by=2)
colfunc<-colorRampPalette(c("red","yellow","springgreen","#196F3D"))
colors <- c(heat.colors(300)[1:100],topo.colors(300)[1:100],terrain.colors(300)[1:100])
for (k in index_plot){
  par(mfrow=c(1,1))
  relative_change_K1 <- (as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K1_list_slope[[k]][[1]][sequence_temp])-list_lower_bounds$lower_bound_list[[k]])/list_lower_bounds$lower_bound_list[[k]]
  relative_change_K1.5 <- (as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K1.5_list_slope[[k]][[1]][sequence_temp])-list_lower_bounds$lower_bound_list[[k]])/list_lower_bounds$lower_bound_list[[k]]
  relative_change_K2 <- (as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K2_list_slope[[k]][[1]][sequence_temp])-list_lower_bounds$lower_bound_list[[k]])/list_lower_bounds$lower_bound_list[[k]]
  relative_change_K2.5 <- (as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K2.5_list_slope[[k]][[1]][sequence_temp])-list_lower_bounds$lower_bound_list[[k]])/list_lower_bounds$lower_bound_list[[k]]
  relative_change_K3 <- (as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K3_list_slope[[k]][[1]][sequence_temp])-list_lower_bounds$lower_bound_list[[k]])/list_lower_bounds$lower_bound_list[[k]]
  relative_change_K3.5 <- (as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K3.5_list_slope[[k]][[1]][sequence_temp])-list_lower_bounds$lower_bound_list[[k]])/list_lower_bounds$lower_bound_list[[k]]
  relative_change_Klog_n <- (as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_Klog_n_list_slope[[k]][[1]][sequence_temp])-list_lower_bounds$lower_bound_list[[k]])/list_lower_bounds$lower_bound_list[[k]]
  relative_change_K4 <- (as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp])-list_lower_bounds$lower_bound_list[[k]])/list_lower_bounds$lower_bound_list[[k]]
  relative_change_K4.5 <- (as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K4.5_list_slope[[k]][[1]][sequence_temp])-list_lower_bounds$lower_bound_list[[k]])/list_lower_bounds$lower_bound_list[[k]]
  relative_change_K5 <- (as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K5_list_slope[[k]][[1]][sequence_temp])-list_lower_bounds$lower_bound_list[[k]])/list_lower_bounds$lower_bound_list[[k]]
  relative_change_K1[which(is.na(relative_change_K1))] = rep(0,length(which(is.na(relative_change_K1))))
  relative_change_K1.5[which(is.na(relative_change_K1.5))] = rep(0,length(which(is.na(relative_change_K1.5))))
  relative_change_K2[which(is.na(relative_change_K2))] = rep(0,length(which(is.na(relative_change_K2))))
  relative_change_K2.5[which(is.na(relative_change_K2.5))] = rep(0,length(which(is.na(relative_change_K2.5))))
  relative_change_K3[which(is.na(relative_change_K3))] = rep(0,length(which(is.na(relative_change_K3))))
  relative_change_K3.5[which(is.na(relative_change_K3.5))] = rep(0,length(which(is.na(relative_change_K3.5))))
  relative_change_Klog_n[which(is.na(relative_change_Klog_n))] = rep(0,length(which(is.na(relative_change_Klog_n))))
  relative_change_K4[which(is.na(relative_change_K4))] = rep(0,length(which(is.na(relative_change_K4))))
  relative_change_K4.5[which(is.na(relative_change_K4.5))] = rep(0,length(which(is.na(relative_change_K4.5))))
  relative_change_K5[which(is.na(relative_change_K5))] = rep(0,length(which(is.na(relative_change_K5))))
  ymax = max(relative_change_K2,relative_change_K3,relative_change_K2.5,relative_change_K3.5,relative_change_K1.5,
             relative_change_K1,relative_change_Klog_n,relative_change_K4,relative_change_K4.5,relative_change_K5)
  ymin = min(relative_change_K2,relative_change_K3,relative_change_K2.5,relative_change_K3.5,relative_change_K1.5,
             relative_change_K1,relative_change_Klog_n,relative_change_K4,relative_change_K4.5,relative_change_K5)
  par(mfrow=c(1,1))
  plot(cste_mult_vect,relative_change_K1,col = colors[200], xlim = c(0,10),ylim = c(ymin,ymax),xlab = "K", ylab = NA, lty = 5,lwd=4,type="p",cex=1,pch=1,cex.lab = 1.8)
  lines(cste_mult_vect,relative_change_K1,col = colors[200], type="l",cex=1,pch=1,cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_K1.5,col = colors[123],cex=1,pch=2, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_K1.5,col = colors[123],cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_K2,col = 3,cex=1,pch=3, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_K2,col = 3,cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_K2.5,col = colors[10],cex=1,pch=4, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_K2.5,col = colors[10],cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_K3,col = colors[160],cex=1,pch=5, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_K3,col = colors[160],cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_K3.5,col = colors[360],cex=1,pch=6, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_K3.5,col = colors[360],cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_K4,col = 6,cex=1,pch=7, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_K4,col = 6,cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_K4.5,col = 7,cex=1,pch=8, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_K4.5,col = 7,cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_Klog_n,col = 8,cex=1,pch=9, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_Klog_n,col = 8,cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_K5,cex=1,pch=10,col = "#660099",xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_K5,cex=1,pch=10,col = "#660099",type="l",cex.lab = 1.3,lty = 1,lwd=4)
  legend("topleft", legend = c(c("K=1","K=1.5","K=2","K=2.5","K=3","K=3.5","K=log(n)","K=4","K=4.5","K=5")),
         col = c(colors[200],colors[123],3,colors[10],colors[160],colors[360],6,7,8,"#660099"), cex=1.8, lwd = rep(4,10), pch = c(1,2,3,4,5,6,7,8,9,10), lty = rep(1,10))
}


###### Relative standard error of the estimated lower bounds functions with respect to K
for (k in index_plot){
  mean_lower_K1 <- c()
  mean_lower_K1.5 <- c()
  mean_lower_K2 <- c()
  mean_lower_K2.5 <- c()
  mean_lower_K3 <- c()
  mean_lower_K3.5 <- c()
  mean_lower_Klog_n <- c()
  mean_lower_K4 <- c()
  mean_lower_K4.5 <- c()
  mean_lower_K5 <- c()
  S_D_lower_K1 <- c()
  S_D_lower_K1.5 <- c()
  S_D_lower_K2 <- c()
  S_D_lower_K2.5 <- c()
  S_D_lower_K3 <- c()
  S_D_lower_K3.5 <- c()
  S_D_lower_Klog_n <- c()
  S_D_lower_K4 <- c()
  S_D_lower_K4.5 <- c()
  S_D_lower_K5 <- c()
  for (l in 1:length(sequence_temp)){
    it = sequence_temp[l]
    temp_K1 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K1_list_slope[[k]][[j]][it]))
    temp_K1.5 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K1.5_list_slope[[k]][[j]][it]))
    temp_K2 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K2_list_slope[[k]][[j]][it]))
    temp_K2.5 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K2.5_list_slope[[k]][[j]][it]))
    temp_K3 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K3_list_slope[[k]][[j]][it]))
    temp_K3.5 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K3.5_list_slope[[k]][[j]][it]))
    temp_Klog_n <- sapply(1:nbr_it_bound, function(j) as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_Klog_n_list_slope[[k]][[j]][it]))
    temp_K4 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K4_list_slope[[k]][[j]][it]))
    temp_K4.5 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K4.5_list_slope[[k]][[j]][it]))
    temp_K5 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K5_list_slope[[k]][[j]][it]))
    mean_lower_K1[l] <- mean(temp_K1)
    mean_lower_K1.5[l] <- mean(temp_K1.5)
    mean_lower_K2[l] <- mean(temp_K2)
    mean_lower_K2.5[l] <- mean(temp_K2.5)
    mean_lower_K3[l] <- mean(temp_K3)
    mean_lower_K3.5[l] <- mean(temp_K3.5)
    mean_lower_Klog_n[l] <- mean(temp_Klog_n)
    mean_lower_K4[l] <- mean(temp_K4)
    mean_lower_K4.5[l] <- mean(temp_K4.5)
    mean_lower_K5[l] <- mean(temp_K5)
    S_D_lower_K1[l] <- sqrt(var(temp_K1))
    S_D_lower_K1.5[l] <- sqrt(var(temp_K1.5))
    S_D_lower_K2[l] <- sqrt(var(temp_K2))
    S_D_lower_K2.5[l] <- sqrt(var(temp_K2.5))
    S_D_lower_K3[l] <- sqrt(var(temp_K3))
    S_D_lower_K3.5[l] <- sqrt(var(temp_K3.5))
    S_D_lower_Klog_n[l] <- sqrt(var(temp_Klog_n))
    S_D_lower_K4[l] <- sqrt(var(temp_K4))
    S_D_lower_K4.5[l] <- sqrt(var(temp_K4.5))
    S_D_lower_K5[l] <- sqrt(var(temp_K5))
  }
  relative_standard_error_K1 <- S_D_lower_K1/mean_lower_K1
  relative_standard_error_K1[which(is.na(relative_standard_error_K1))] = rep(0,length(which(is.na(relative_standard_error_K1))))
  relative_standard_error_K1.5 <- S_D_lower_K1.5/mean_lower_K1.5
  relative_standard_error_K1.5[which(is.na(relative_standard_error_K1.5))] = rep(0,length(which(is.na(relative_standard_error_K1.5))))
  relative_standard_error_K2 <- S_D_lower_K2/mean_lower_K2
  relative_standard_error_K2[which(is.na(relative_standard_error_K2))] = rep(0,length(which(is.na(relative_standard_error_K2))))
  relative_standard_error_K2.5 <- S_D_lower_K2.5/mean_lower_K2.5
  relative_standard_error_K2.5[which(is.na(relative_standard_error_K2.5))] = rep(0,length(which(is.na(relative_standard_error_K2.5))))
  relative_standard_error_K3 <- S_D_lower_K3/mean_lower_K3
  relative_standard_error_K3[which(is.na(relative_standard_error_K3))] = rep(0,length(which(is.na(relative_standard_error_K3))))
  relative_standard_error_K3.5 <- S_D_lower_K3.5/mean_lower_K3.5
  relative_standard_error_K3.5[which(is.na(relative_standard_error_K3.5))] = rep(0,length(which(is.na(relative_standard_error_K3.5))))
  relative_standard_error_Klog_n <- S_D_lower_Klog_n/mean_lower_Klog_n
  relative_standard_error_Klog_n[which(is.na(relative_standard_error_Klog_n))] = rep(0,length(which(is.na(relative_standard_error_Klog_n))))
  relative_standard_error_K4 <- S_D_lower_K4/mean_lower_K4
  relative_standard_error_K4[which(is.na(relative_standard_error_K4))] = rep(0,length(which(is.na(relative_standard_error_K4))))
  relative_standard_error_K4.5 <- S_D_lower_K4.5/mean_lower_K4.5
  relative_standard_error_K4.5[which(is.na(relative_standard_error_K4.5))] = rep(0,length(which(is.na(relative_standard_error_K4.5))))
  relative_standard_error_K5 <- S_D_lower_K5/mean_lower_K5
  relative_standard_error_K5[which(is.na(relative_standard_error_K5))] = rep(0,length(which(is.na(relative_standard_error_K5))))
  ymax = max(relative_standard_error_K1,relative_standard_error_K1.5,relative_standard_error_K2,
             relative_standard_error_K2.5,relative_standard_error_K3,relative_standard_error_K3.5,relative_standard_error_Klog_n,
             relative_standard_error_K4,relative_standard_error_K4.5,relative_standard_error_K5)
  par(mfrow=c(1,1))
  plot(cste_mult_vect,relative_standard_error_K1,col = colors[200], xlim = c(0,10),ylim = c(0,ymax),xlab = "K", ylab = NA, lty = 5,lwd=4,type="p",cex=1,pch=1,cex.lab = 1.8)
  lines(cste_mult_vect,relative_standard_error_K1,col = colors[200], type="l",cex=1,pch=1,cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_K1.5,col = colors[123],cex=1,pch=2, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_K1.5,col = colors[123],cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_K2,col = 3,cex=1,pch=3, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_K2,col = 3,cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_K2.5,col = colors[10],cex=1,pch=4, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_K2.5,col = colors[10],cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_K3,col = colors[160],cex=1,pch=5, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_K3,col = colors[160],cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_K3.5,col = colors[360],cex=1,pch=6, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_K3.5,col = colors[360],cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_K4,col = 6,cex=1,pch=7, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_K4,col = 6,cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_K4.5,col = 7,cex=1,pch=8, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_K4.5,col = 7,cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_Klog_n,col = 8,cex=1,pch=9, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_Klog_n,col = 8,cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_K5,cex=1,pch=10,col = "#660099",xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_K5,cex=1,pch=10,col = "#660099",type="l",cex.lab = 1.3,lty = 1,lwd=4)
  legend("topleft", legend = c(c("K=1","K=1.5","K=2","K=2.5","K=3","K=3.5","K=log(n)","K=4","K=4.5","K=5")),
         col = c(colors[200],colors[123],3,colors[10],colors[160],colors[360],6,7,8,"#660099"), cex=1.8, lwd = rep(4,10), pch = c(1,2,3,4,5,6,7,8,9,10), lty = rep(1,10))
}


###### Relative change of the estimated upper bounds functions with respect to K
which_cste_sup2 <- which(cste_mult_vect>=2)
sequence_temp <- seq(1,length(list_lower_bounds$lower_bound_list_hat_one_dataset_K2_list_slope[[1]][[1]]), by=2)
colfunc<-colorRampPalette(c("red","yellow","springgreen","#196F3D"))
colors <- c(heat.colors(300)[1:100],topo.colors(300)[1:100],terrain.colors(300)[1:100])
for (k in index_plot){
  par(mfrow=c(1,1))
  relative_change_K1 <- (as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K1_list_slope[[k]][[1]][sequence_temp])-list_upper_bounds$upper_bound_list[[k]])/list_upper_bounds$upper_bound_list[[k]]
  relative_change_K1.5 <- (as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K1.5_list_slope[[k]][[1]][sequence_temp])-list_upper_bounds$upper_bound_list[[k]])/list_upper_bounds$upper_bound_list[[k]]
  relative_change_K2 <- (as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K2_list_slope[[k]][[1]][sequence_temp])-list_upper_bounds$upper_bound_list[[k]])/list_upper_bounds$upper_bound_list[[k]]
  relative_change_K2.5 <- (as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K2.5_list_slope[[k]][[1]][sequence_temp])-list_upper_bounds$upper_bound_list[[k]])/list_upper_bounds$upper_bound_list[[k]]
  relative_change_K3 <- (as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K3_list_slope[[k]][[1]][sequence_temp])-list_upper_bounds$upper_bound_list[[k]])/list_upper_bounds$upper_bound_list[[k]]
  relative_change_K3.5 <- (as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K3.5_list_slope[[k]][[1]][sequence_temp])-list_upper_bounds$upper_bound_list[[k]])/list_upper_bounds$upper_bound_list[[k]]
  relative_change_Klog_n <- (as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_Klog_n_list_slope[[k]][[1]][sequence_temp])-list_upper_bounds$upper_bound_list[[k]])/list_upper_bounds$upper_bound_list[[k]]
  relative_change_K4 <- (as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp])-list_upper_bounds$upper_bound_list[[k]])/list_upper_bounds$upper_bound_list[[k]]
  relative_change_K4.5 <- (as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4.5_list_slope[[k]][[1]][sequence_temp])-list_upper_bounds$upper_bound_list[[k]])/list_upper_bounds$upper_bound_list[[k]]
  relative_change_K5 <- (as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K5_list_slope[[k]][[1]][sequence_temp])-list_upper_bounds$upper_bound_list[[k]])/list_upper_bounds$upper_bound_list[[k]]
  relative_change_K1[which(is.na(relative_change_K1))] = rep(0,length(which(is.na(relative_change_K1))))
  relative_change_K1.5[which(is.na(relative_change_K1.5))] = rep(0,length(which(is.na(relative_change_K1.5))))
  relative_change_K2[which(is.na(relative_change_K2))] = rep(0,length(which(is.na(relative_change_K2))))
  relative_change_K2.5[which(is.na(relative_change_K2.5))] = rep(0,length(which(is.na(relative_change_K2.5))))
  relative_change_K3[which(is.na(relative_change_K3))] = rep(0,length(which(is.na(relative_change_K3))))
  relative_change_K3.5[which(is.na(relative_change_K3.5))] = rep(0,length(which(is.na(relative_change_K3.5))))
  relative_change_Klog_n[which(is.na(relative_change_Klog_n))] = rep(0,length(which(is.na(relative_change_Klog_n))))
  relative_change_K4[which(is.na(relative_change_K4))] = rep(0,length(which(is.na(relative_change_K4))))
  relative_change_K4.5[which(is.na(relative_change_K4.5))] = rep(0,length(which(is.na(relative_change_K4.5))))
  relative_change_K5[which(is.na(relative_change_K5))] = rep(0,length(which(is.na(relative_change_K5))))
  ymax = max(relative_change_K2,relative_change_K3,relative_change_K2.5,relative_change_K3.5,relative_change_K1.5,
             relative_change_K1,relative_change_Klog_n,relative_change_K4,relative_change_K4.5,relative_change_K5)
  ymin = min(relative_change_K2,relative_change_K3,relative_change_K2.5,relative_change_K3.5,relative_change_K1.5,
             relative_change_K1,relative_change_Klog_n,relative_change_K4,relative_change_K4.5,relative_change_K5)
  par(mfrow=c(1,1))
  plot(cste_mult_vect,relative_change_K1,col = colors[200], xlim = c(0,10),ylim = c(ymin,ymax),xlab = "K", ylab = NA, lty = 5,lwd=4,type="p",cex=1,pch=1,cex.lab = 1.8)
  lines(cste_mult_vect,relative_change_K1,col = colors[200], type="l",cex=1,pch=1,cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_K1.5,col = colors[123],cex=1,pch=2, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_K1.5,col = colors[123],cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_K2,col = 3,cex=1,pch=3, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_K2,col = 3,cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_K2.5,col = colors[10],cex=1,pch=4, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_K2.5,col = colors[10],cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_K3,col = colors[160],cex=1,pch=5, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_K3,col = colors[160],cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_K3.5,col = colors[360],cex=1,pch=6, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_K3.5,col = colors[360],cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_K4,col = 6,cex=1,pch=7, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_K4,col = 6,cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_K4.5,col = 7,cex=1,pch=8, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_K4.5,col = 7,cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_Klog_n,col = 8,cex=1,pch=9, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_Klog_n,col = 8,cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_change_K5,cex=1,pch=10,col = "#660099",xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_change_K5,cex=1,pch=10,col = "#660099",type="l",cex.lab = 1.3,lty = 1,lwd=4)
  legend("topleft", legend = c(c("K=1","K=1.5","K=2","K=2.5","K=3","K=3.5","K=log(n)","K=4","K=4.5","K=5")),
         col = c(colors[200],colors[123],3,colors[10],colors[160],colors[360],6,7,8,"#660099"), cex=1.8, lwd = rep(4,10), pch = c(1,2,3,4,5,6,7,8,9,10), lty = rep(1,10))
}



###### Relative standard error of the estimated upper bounds functions with respect to K
for (k in index_plot){
  mean_upper_K1 <- c()
  mean_upper_K1.5 <- c()
  mean_upper_K2 <- c()
  mean_upper_K2.5 <- c()
  mean_upper_K3 <- c()
  mean_upper_K3.5 <- c()
  mean_upper_Klog_n <- c()
  mean_upper_K4 <- c()
  mean_upper_K4.5 <- c()
  mean_upper_K5 <- c()
  S_D_upper_K1 <- c()
  S_D_upper_K1.5 <- c()
  S_D_upper_K2 <- c()
  S_D_upper_K2.5 <- c()
  S_D_upper_K3 <- c()
  S_D_upper_K3.5 <- c()
  S_D_upper_Klog_n <- c()
  S_D_upper_K4 <- c()
  S_D_upper_K4.5 <- c()
  S_D_upper_K5 <- c()
  for (l in 1:length(sequence_temp)){
    it = sequence_temp[l]
    temp_K1 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K1_list_slope[[k]][[j]][it]))
    temp_K1.5 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K1.5_list_slope[[k]][[j]][it]))
    temp_K2 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K2_list_slope[[k]][[j]][it]))
    temp_K2.5 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K2.5_list_slope[[k]][[j]][it]))
    temp_K3 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K3_list_slope[[k]][[j]][it]))
    temp_K3.5 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K3.5_list_slope[[k]][[j]][it]))
    temp_Klog_n <- sapply(1:nbr_it_bound, function(j) as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_Klog_n_list_slope[[k]][[j]][it]))
    temp_K4 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[j]][it]))
    temp_K4.5 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4.5_list_slope[[k]][[j]][it]))
    temp_K5 <- sapply(1:nbr_it_bound, function(j) as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K5_list_slope[[k]][[j]][it]))
    mean_upper_K1[l] <- mean(temp_K1)
    mean_upper_K1.5[l] <- mean(temp_K1.5)
    mean_upper_K2[l] <- mean(temp_K2)
    mean_upper_K2.5[l] <- mean(temp_K2.5)
    mean_upper_K3[l] <- mean(temp_K3)
    mean_upper_K3.5[l] <- mean(temp_K3.5)
    mean_upper_Klog_n[l] <- mean(temp_Klog_n)
    mean_upper_K4[l] <- mean(temp_K4)
    mean_upper_K4.5[l] <- mean(temp_K4.5)
    mean_upper_K5[l] <- mean(temp_K5)
    S_D_upper_K1[l] <- sqrt(var(temp_K1))
    S_D_upper_K1.5[l] <- sqrt(var(temp_K1.5))
    S_D_upper_K2[l] <- sqrt(var(temp_K2))
    S_D_upper_K2.5[l] <- sqrt(var(temp_K2.5))
    S_D_upper_K3[l] <- sqrt(var(temp_K3))
    S_D_upper_K3.5[l] <- sqrt(var(temp_K3.5))
    S_D_upper_Klog_n[l] <- sqrt(var(temp_Klog_n))
    S_D_upper_K4[l] <- sqrt(var(temp_K4))
    S_D_upper_K4.5[l] <- sqrt(var(temp_K4.5))
    S_D_upper_K5[l] <- sqrt(var(temp_K5))
  }
  relative_standard_error_K1 <- S_D_upper_K1/mean_upper_K1
  relative_standard_error_K1[which(is.na(relative_standard_error_K1))] = rep(0,length(which(is.na(relative_standard_error_K1))))
  relative_standard_error_K1.5 <- S_D_upper_K1.5/mean_upper_K1.5
  relative_standard_error_K1.5[which(is.na(relative_standard_error_K1.5))] = rep(0,length(which(is.na(relative_standard_error_K1.5))))
  relative_standard_error_K2 <- S_D_upper_K2/mean_upper_K2
  relative_standard_error_K2[which(is.na(relative_standard_error_K2))] = rep(0,length(which(is.na(relative_standard_error_K2))))
  relative_standard_error_K2.5 <- S_D_upper_K2.5/mean_upper_K2.5
  relative_standard_error_K2.5[which(is.na(relative_standard_error_K2.5))] = rep(0,length(which(is.na(relative_standard_error_K2.5))))
  relative_standard_error_K3 <- S_D_upper_K3/mean_upper_K3
  relative_standard_error_K3[which(is.na(relative_standard_error_K3))] = rep(0,length(which(is.na(relative_standard_error_K3))))
  relative_standard_error_K3.5 <- S_D_upper_K3.5/mean_upper_K3.5
  relative_standard_error_K3.5[which(is.na(relative_standard_error_K3.5))] = rep(0,length(which(is.na(relative_standard_error_K3.5))))
  relative_standard_error_Klog_n <- S_D_upper_Klog_n/mean_upper_Klog_n
  relative_standard_error_Klog_n[which(is.na(relative_standard_error_Klog_n))] = rep(0,length(which(is.na(relative_standard_error_Klog_n))))
  relative_standard_error_K4 <- S_D_upper_K4/mean_upper_K4
  relative_standard_error_K4[which(is.na(relative_standard_error_K4))] = rep(0,length(which(is.na(relative_standard_error_K4))))
  relative_standard_error_K4.5 <- S_D_upper_K4.5/mean_upper_K4.5
  relative_standard_error_K4.5[which(is.na(relative_standard_error_K4.5))] = rep(0,length(which(is.na(relative_standard_error_K4.5))))
  relative_standard_error_K5 <- S_D_upper_K5/mean_upper_K5
  relative_standard_error_K5[which(is.na(relative_standard_error_K5))] = rep(0,length(which(is.na(relative_standard_error_K5))))
  ymax = max(relative_standard_error_K1,relative_standard_error_K1.5,relative_standard_error_K2,
             relative_standard_error_K2.5,relative_standard_error_K3,relative_standard_error_K3.5,relative_standard_error_Klog_n,
             relative_standard_error_K4,relative_standard_error_K4.5,relative_standard_error_K5)
  par(mfrow=c(1,1))
  plot(cste_mult_vect,relative_standard_error_K1,col = colors[200], xlim = c(0,10),ylim = c(0,ymax),xlab = "K", ylab = NA, lty = 5,lwd=4,type="p",cex=1,pch=1,cex.lab = 1.8)
  lines(cste_mult_vect,relative_standard_error_K1,col = colors[200], type="l",cex=1,pch=1,cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_K1.5,col = colors[123],cex=1,pch=2, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_K1.5,col = colors[123],cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_K2,col = 3,cex=1,pch=3, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_K2,col = 3,cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_K2.5,col = colors[10],cex=1,pch=4, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_K2.5,col = colors[10],cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_K3,col = colors[160],cex=1,pch=5, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_K3,col = colors[160],cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_K3.5,col = colors[360],cex=1,pch=6, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_K3.5,col = colors[360],cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_K4,col = 6,cex=1,pch=7, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_K4,col = 6,cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_K4.5,col = 7,cex=1,pch=8, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_K4.5,col = 7,cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_Klog_n,col = 8,cex=1,pch=9, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_Klog_n,col = 8,cex=1,pch=1, type="l",cex.lab = 1.3,lty = 1,lwd=4)
  lines(cste_mult_vect,relative_standard_error_K5,cex=1,pch=10,col = "#660099",xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p",cex.lab = 1.3)
  lines(cste_mult_vect,relative_standard_error_K5,cex=1,pch=10,col = "#660099",type="l",cex.lab = 1.3,lty = 1,lwd=4)
  legend("topleft", legend = c(c("K=1","K=1.5","K=2","K=2.5","K=3","K=3.5","K=log(n)","K=4","K=4.5","K=5")),
         col = c(colors[200],colors[123],3,colors[10],colors[160],colors[360],6,7,8,"#660099"), cex=1.8, lwd = rep(4,10), pch = c(1,2,3,4,5,6,7,8,9,10), lty = rep(1,10))
}



#####################################################################################################################
###### Plottings for deterministic non ordered model collections and random model collections
#####################################################################################################################

###### For deterministic non ordered model collections
cste_for_estimation <- 4
sequence_temp <- seq(1,length(list_lower_bounds$lower_bound_list_hat_one_dataset_K2_list_slope[[1]][[1]]), by=2)
# Plot of the estimated upper bound of FDR, theoretical upper bound of FDR and empirical FDR values with respect to K
par(mfrow=c(1,1))
k = 1
quantile_asympt_95 <- 1.96      # for the confidence interval of the FDR empirical estimation. Application of the central limit theorem
vect_empirical_FDR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))     # empirical FDR
variance_estimation_FDR <- sapply(1:length(cste_mult_vect), function(it) (1/(nbr_it-1))*sum((PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]-mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))^2))    # empirical variance
FDP_confidence_interval_inf <- vect_empirical_FDR - (variance_estimation_FDR*quantile_asympt_95)/sqrt(nbr_it)     # lower bound of the confidence interval for the FDR empirical estimation
FDP_confidence_interval_sup <- vect_empirical_FDR + (variance_estimation_FDR*quantile_asympt_95)/sqrt(nbr_it)      # upper bound of the confidence interval for the FDR empirical estimation
plot(cste_mult_vect,FDP_confidence_interval_inf, col = 1,type = 'l',xlim = c(0.1,10),
     ylim=c(0,max(list_upper_bounds$upper_bound_list[[1]],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]))),ylab = NA, xlab = "K",
     lty = 5, lwd=1, cex.lab = 1.3)
lines(cste_mult_vect,FDP_confidence_interval_sup, col = 1, type = 'l', lty = 1, lwd=1)
lines(cste_mult_vect,vect_empirical_FDR,col = 2, type="p", lty = 1, lwd=1,pch=1)
lines(cste_mult_vect,vect_empirical_FDR,col = 2,type="l", lty = 1, lwd=1,pch=1)
lines(cste_mult_vect,list_upper_bounds$upper_bound_list[[k]], col = 4, type="p", lty = 1, lwd=1,pch=1)
lines(cste_mult_vect,list_upper_bounds$upper_bound_list[[k]], col = 4, type="l", lty = 1, lwd=1,pch=1)
lines(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]),ylim=c(0,1),xlab = "K",ylab = "estimated upper bound of FDR on one dataset",type="p", col = 6, lty = 1, lwd=1,pch=1)
lines(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]),type="l", col = 6, lty = 1, lwd=1,pch=1)
for (k in index_plot[2:length(index_plot)]){
  lines(cste_mult_vect,FDP_confidence_interval_inf, col = 1,type = 'l',xlim = c(0.1,10),
        ylim=c(0,max(list_upper_bounds$upper_bound_list[[k]])),ylab = NA, xlab = "K",lty = 5, lwd=1, cex.lab = 1.3)
  lines(cste_mult_vect,FDP_confidence_interval_sup, col = 1, type = 'l', lty = 1, lwd=1)
  lines(cste_mult_vect,vect_empirical_FDR,col = 2, type="p", lty = 1, lwd=1,pch=k)
  lines(cste_mult_vect,vect_empirical_FDR,col = 2,type="l", lty = 1, lwd=1,pch=k)
  lines(cste_mult_vect,list_upper_bounds$upper_bound_list[[k]], col = 4, type="p", lty = 1, lwd=1,pch=k)
  lines(cste_mult_vect,list_upper_bounds$upper_bound_list[[k]], col = 4, type="l", lty = 1, lwd=1,pch=k)
  lines(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]),ylim=c(0,1),xlab = "K",ylab = "estimated upper bound of FDR on one dataset",type="p", col = 6, lty = 1, lwd=1,pch=k)
  lines(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]),type="l", col = 6, lty = 1, lwd=1,pch=k)
}
legend("topright", legend = c("empirical FDR","B","estimated B","with 4","permutation_1_10","exchanging9_10_11_12","exchanging1_2_16_20","perm_10","perm_12","perm_15"),
       col = c(2,4,6,1,1,1,1), cex=1.1, lwd = c(2,2,2,0,2,2,2,2,2,2),
       pch = c(26,26,26,26,1:6),lty = c(1,1,1,0,0,0,0,0,0,0))


# Plot of the estimated difference in predictions and empirical difference in predictions with respect to K
par(mfrow=c(1,1))
k = 1
vect_empirical_PR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$PR_hat[,it]))    # empirical PR 
vect_empirical_diff_PR <- vect_empirical_PR - vect_empirical_PR[which(cste_mult_vect==2)]       # empirical difference in predictions
vect_estimation_diff_PR <- PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope       # estimated difference in predictions
min_diff_PR <- min(vect_empirical_diff_PR,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope)
max_diff_PR <- max(vect_empirical_diff_PR,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope)
plot(cste_mult_vect,vect_empirical_diff_PR, col = 4, type="p",ylim = c(min_diff_PR,max_diff_PR),xlab = "K",ylab=NA,lty = 5, lwd=1, cex.lab = 1.3,pch=1)
lines(cste_mult_vect,vect_empirical_diff_PR, col = 4,type="l", lty = 1, lwd=1,pch=1)
lines(cste_mult_vect,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope, col = "#660099", type="p", lty = 1, lwd=1,pch=1)
lines(cste_mult_vect,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope, col = "#660099", type="l", lty = 1, lwd=1,pch=1)
for (k in index_plot[2:length(index_plot)]){
  vect_empirical_PR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$PR_hat[,it]))
  vect_empirical_PR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$PR_hat[,it]))
  vect_empirical_diff_PR <- vect_empirical_PR - vect_empirical_PR[which(cste_mult_vect==2)]
  vect_estimation_diff_PR <- PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope
  min_diff_PR <- min(vect_empirical_diff_PR,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope)
  max_diff_PR <- max(vect_empirical_diff_PR,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope)
  lines(cste_mult_vect,vect_empirical_diff_PR, col = 4, type="p",ylim = c(min_diff_PR,max_diff_PR),xlab = "K",ylab=NA,lty = 5, lwd=1, cex.lab = 1.3,pch=3)
  lines(cste_mult_vect,vect_empirical_diff_PR, col = 4,type="l", lty = 1, lwd=1,pch=3)
  lines(cste_mult_vect,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope, col = "#660099", type="p", lty = 1, lwd=1,pch=3)
  lines(cste_mult_vect,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope, col = "#660099", type="l", lty = 1, lwd=1,pch=3)
}
legend("topright", legend = c("empirical","difference in","predictions","estimated","difference in", "predictions","permutation_1_10","exchanging9_10_11_12","exchanging1_2_16_20","perm_10","perm_12","perm_15"),
       col = c(4,4,4,"#660099","#660099","#660099",1,1,1,1,1,1), cex=0.7, lwd = c(2,0,0,2,0,0,2,2,2,2,2,2),
       pch = c(26,26,26,26,26,26,1:6),lty = c(1,0,0,1,0,0,0,0,0,0,0,0))


# Plot of the empirical FDR and the empirical difference in predictions with respect to K
for (k in index_plot){
  quantile_asympt_95 <- 1.96      # for the confidence interval of the FDR empirical estimation. Application of the central limit theorem
  vect_empirical_FDR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))     # empirical FDR
  variance_estimation_FDR <- sapply(1:length(cste_mult_vect), function(it) (1/(nbr_it-1))*sum((PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]-mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))^2))    # empirical variance
  FDP_confidence_interval_inf <- vect_empirical_FDR - (variance_estimation_FDR*quantile_asympt_95)/sqrt(nbr_it)     # lower bound of the confidence interval for the FDR empirical estimation
  FDP_confidence_interval_sup <- vect_empirical_FDR + (variance_estimation_FDR*quantile_asympt_95)/sqrt(nbr_it)      # upper bound of the confidence interval for the FDR empirical estimation
  par(mfrow=c(1,2))
  plot(cste_mult_vect,FDP_confidence_interval_inf, col = 1,type = 'l',xlim = c(0.1,10),ylim=c(0,max(list_upper_bounds$upper_bound_list[[k]])),xlab = "K",ylab = "empirical FDR and bounds", lty = 5)
  lines(cste_mult_vect,FDP_confidence_interval_sup, col = 1, type = 'l', lty = 5)
  lines(cste_mult_vect,vect_empirical_FDR,col = 2, type="p", lty = 5)
  lines(cste_mult_vect,vect_empirical_FDR,col = 2,type="l")
  lines(cste_mult_vect,list_upper_bounds$upper_bound_list[[k]], col = 4, type="p")
  lines(cste_mult_vect,list_upper_bounds$upper_bound_list[[k]], col = 4, type="l")
  lines(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]),ylim=c(0,1),xlab = "K",ylab = "estimated upper bound of FDR on one dataset",type="p", col = 6)
  lines(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]),type="l", col = 6)
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","Hat_upper_bound(beta(Hat(m)(4)))"), col = c(2,6), lty = 1, cex=0.7, lwd = c(1,1))
  vect_empirical_PR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$PR_hat[,it]))     # empirical PR
  vect_empirical_diff_PR <- vect_empirical_PR - vect_empirical_PR[which(cste_mult_vect==2)]       # empirical difference in predictions
  min_diff_PR <- min(vect_empirical_diff_PR,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope)
  max_diff_PR <- max(vect_empirical_diff_PR,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope)
  plot(cste_mult_vect,vect_empirical_diff_PR, col = 4, type="p", lty = 5,ylim = c(min_diff_PR,max_diff_PR),xlab = "K",ylab = "empirical PR and estimated PR")
  lines(cste_mult_vect,vect_empirical_diff_PR, col = 4,type="l")
  lines(cste_mult_vect,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope, col = "#660099", type="p")
  lines(cste_mult_vect,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope, col = "#660099", type="l")
  legend("topright", legend = c("K -> diff_empirical PR(Hat(m)(K))","K -> diff_estimated PR(Hat(m)(K))"), col = c(4,"#660099"), lty = 1, cex=0.7, lwd = c(1,1))
}




# #####################################################################################################################
# ###### Study of Algorithm 1 of the article to calibrate the hyperparameter K
# #####################################################################################################################
## To calibrate the hyperparameter K providing a trade-off betweeen the PR and the FDR cost functions, we use \Tilde{K} = 4 for the estimated upper bound of the FDR function.
## This code below provides the intervals of K with the lowest values of the estimated upper bound of the FDR and of the estimated difference function of the PR,
## as well as the intersection of the two intervals and the smallest of the values.
## To evaluate the obtained \Hat{K}, the algorithm below provides the intervals of K of the lowest values of the empirical PR and FDR, as well as the intersection of both intervals.

# estimated FDR and PR with \Tilde{K} = 4
cste_for_estimation <- 4
alpha_FDR_our_algo <- 0.001
gamma_PR <- 0.1
I_1 <- list()
I_2 <- list()
intersection <- list()
min_intersection <- c()
for (k in index_plot){
  vect_estimation_FDR <- as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][seq(1,length(list_lower_bounds$lower_bound_list_hat_one_dataset_K2_list_slope[[1]][[1]]), by=2)])    # estimated FDR
  vect_estimation_diff_PR <- PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope       # estimated difference in predictions
  I_1[[k]] <- cste_mult_vect[which(vect_estimation_FDR<alpha_FDR_our_algo)]
  I_2[[k]] <- cste_mult_vect[which(abs(vect_estimation_diff_PR)<gamma_PR*list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]]$kappa_slope)]
  intersection[[k]] <- intersect(I_1[[k]],I_2[[k]])
  if (length(intersection[[k]])>0){
    min_intersection <- c(min_intersection,min(intersect(I_1[[k]],I_2[[k]])))
  }else{
    min_intersection <- c(min_intersection,NA)
  }
}
# empirical FDR and difference in predictions
I_1_theo <- list()
I_2_theo <- list()
intersection_theo <- list()
for (k in index_plot){
  vect_empirical_FDR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))      # empirical FDR  
  vect_empirical_PR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$PR_hat[,it]))        # empirical PR
  vect_empirical_diff_PR <- vect_empirical_PR - vect_empirical_PR[which(cste_mult_vect==2)]        # empirical difference in predictions
  I_1_theo[[k]] <- cste_mult_vect[which(vect_empirical_FDR<alpha_FDR_our_algo)]
  I_2_theo[[k]] <- cste_mult_vect[which(abs(vect_empirical_diff_PR)<gamma_PR*sigma_2)]
  intersection_theo[[k]] <- intersect(I_1_theo[[k]],I_2_theo[[k]])
}



####################################################################################################################
##### Comparaisons of performances of our algorithm and some existing variable selection methods
####################################################################################################################
##### Application of the variable selection procedures
# To compare the model selection provided by our algorithm, we process a comparaison with some existing model selection procedures:
# 1) LinSelect with the R function tuneLasso of the R package LINselect
# 2) K-fold cross-validation with R function train of the R package caret. We use the option method = "lm" and 50 Kfold
# 3) The knockoffs method with the R function knockoff.filter of the R package knockoff
# 4) our algorithm

l=11         # the dimension of the true model is fixed to 10
cste_for_estimation <- 4    # for our algorithm: \Tilde{K} = 4     
gamma_PR <- 0.1    # threshold for our algorithm
alpha_FDR <- 0.05    # threshold for our algorithm
Y_bounds <- Y_bounds_matrix[,l]
result_LinSelect <- list()
result_KCV <- list()
result_Knockoffs <- list()
result_our_algorithm <- list()
if (random_collection == TRUE){
  if (min(n,p)<p){     # knockoff random model collections do not exist
    path <- list(path_collection[[1]],path_collection[[2]],path_collection[[3]])
  }else{
    path <- path_collection
  }
}
if (ordered_variable == TRUE){
  path <- list()
  path[[1]] <- path_collection_bounds[[l]]
}
for (x in 1:length(path)){
  result_LinSelect[[x]] <- list()
  result_KCV[[x]] <- list()
  result_Knockoffs[[x]] <- list()
  pb = txtProgressBar(min = 1, max = nbr_it_bound, initial = 1)
  for (k in seq(1,n*nbr_it_bound,by=n)){
    j = which(seq(1,n*nbr_it_bound,by=n) ==k)
    ## LinSelect
    Y_s <- Y_select_matrix[,l][k:(k+(n-1))]
    data_s <- data_select
    res <- path[[x]][[j]]
    result_LinSelect[[x]][[j]] <- LinSelect_function(res,n,p,dmax_LIN,max.steps_LIN,coeff_mult_proj_LIN,coeff_mult_pen_LIN,Y_s,data_s)
    ## K-fold CV
    if (ordered_variable == TRUE){
      order_index <- 1:p
      support_num <- list()    # list for the supports of variable index in the collection
      for(y in 2:(p+1)){
        support_num[[y]] <- order_index[1:(y-1)]     # construction of the model support
      }
      data_path <- lapply(support_num, function(y) data_select[,y])
    }
    if (random_collection == TRUE){
      data_path <- lapply(path[[x]][[j]]$val5, function(y) data_select[,y])
    }
    index_temp <- 1
    if (length(data_path[[1]])==0){
      index_temp <- index_temp+1
    }
    if (length(data_path[[2]])==n){
      data_path[[2]] <- matrix(data_path[[2]],nrow = n, ncol = 1)
      if (random_collection == TRUE){
        colnames(data_path[[2]]) <- paste("covariable",path[[x]][[j]]$val5[[2]])
      }
      if (ordered_variable == TRUE){
        colnames(data_path[[2]]) <- paste("covariable",support_num[[2]])
      }
    }
    res_CV <- lapply(index_temp:length(data_path), function(i) train(y = Y_bounds[k:(k+(n-1))], x = data_path[[i]], method = "lm", trControl = trainControl(method = "cv", number = K_fold)))
    result_KCV[[x]][[j]] <- which.min(sapply(res_CV, function(y) y$results$RMSE))
    ## Knockoffs
    if (min(n,p) >= p){
      result_Knockoffs[[x]][[j]] <- knockoff.filter(X = data_select, y = Y_bounds[k:(k+(n-1))], knockoffs = create.second_order,statistic = stat.lasso_lambdasmax, fdr = alpha_FDR)
    }
    setTxtProgressBar(pb,j)
  }
}
## Our algorithm
result_our_algorithm <- trade_off_PR_FDR(cste_for_estimation,alpha_FDR,gamma_PR,list_upper_bounds,nbr_it_bound,choice_random_collection = 1:length(path),PR_FDP_hat_K,cste_mult_vect,sigma_2)
results_selection <- list(LinSelect = result_LinSelect, KCV = result_KCV, Knockoffs = result_Knockoffs,Our_algorithm = result_our_algorithm)
save(results_selection, file = paste("results_selection",config_data,sep="_"))
# save(results_selection, file = paste("results_selection_smaller_noise",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# save(results_selection, file = paste("results_selection_larger_noise_closed_values",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article



##### Comparison of the variable selection procedures
load(paste("results_selection",config_data,sep="_"))
# load(paste("results_selection_smaller_noise",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("results_selection_larger_noise_closed_values",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article

if (ordered_variable == TRUE){
  index_temp <- 1    # only one model collection: the deterministic nested model collection
}
if (random_collection == TRUE){
  if (min(n,p) >= p){
    index_temp <- 1:4    # 4 random model collections: from Bolasso, SLOPE and RF
  }else{
    index_temp <- 1:3      # knockoff random model collections do not exist
  }
}

# LinSelect for the different model collections: Bolasso, Slope, RF, Knockoffs
PR_FDR <- list()
if (ordered_variable == TRUE){
  PR_FDR[[1]] <-  PR_FDP_hat_K[[11]]
}
if (random_collection == TRUE){
  PR_FDR <- PR_FDP_hat_K
}
mhat_LinSelect <- sapply(index_temp, function(x) mean(unlist(results_selection$LinSelect[[x]])))
PR_LinSelect <- sapply(index_temp, function(x) mean(sapply(1:nbr_it_bound, function(y) mean(sapply(1:nbr_it, function(j) PR_FDR[[x]]$path$PR_path[[j]][results_selection$LinSelect[[x]][[y]]])))))
FDR_LinSelect <- sapply(index_temp, function(x) mean(sapply(1:nbr_it_bound, function(y) mean(sapply(1:nbr_it, function(j) PR_FDR[[x]]$path$FDP_path[[j]][results_selection$LinSelect[[x]][[y]]])))))

# KCV for the different model collections: Bolasso, Slope, RF, Knockoffs
mhat_KCV <- sapply(index_temp, function(x) mean(unlist(results_selection$KCV[[x]])))
PR_KCV <- sapply(index_temp, function(x) mean(sapply(1:nbr_it_bound, function(y) mean(sapply(1:nbr_it, function(j) PR_FDR[[x]]$path$PR_path[[j]][results_selection$KCV[[x]][[y]]])))))
FDR_KCV <- sapply(index_temp, function(x) mean(sapply(1:nbr_it_bound, function(y) mean(sapply(1:nbr_it, function(j) PR_FDR[[x]]$path$FDP_path[[j]][results_selection$KCV[[x]][[y]]])))))

# Knockoffs for the different model collections: Bolasso, Slope, RF, Knockoffs
if (min(n,p) >= p){
  mhat_Knockoffs <- sapply(index_temp, function(x) mean(sapply(1:100, function(y) length(results_selection$Knockoffs[[x]][[y]]$selected)+1)))
  PR_Knockoffs <- sapply(index_temp, function(x) mean(sapply(1:nbr_it_bound, function(y) mean(sapply(1:nbr_it, function(j) PR_FDR[[x]]$path$PR_path[[j]][length(results_selection$Knockoffs[[x]][[y]]$selected)+1])))))
  FDR_Knockoffs <- sapply(index_temp, function(x) mean(sapply(1:nbr_it_bound, function(y) mean(sapply(1:nbr_it, function(j) PR_FDR[[x]]$path$FDP_path[[j]][length(results_selection$Knockoffs[[x]][[y]]$selected)+1])))))
}else{      # the knockoff method is not adapted to the n < p case
  mhat_Knockoffs <- rep(NA,3)
  PR_Knockoffs <- rep(NA,3)
  FDR_Knockoffs <- rep(NA,3)
}

# Our algorithm for the different model collections: Bolasso, Slope, RF, Knockoffs
mhat_our_algorithm <- sapply(index_temp, function(x) mean(sapply(1:nbr_it_bound, function(y) mean(PR_FDR[[x]]$empirical$model_hat[,which(results_selection$Our_algorithm$K_hat[[x]][[y]]==cste_mult_vect)]))))
PR_our_algorithm <- sapply(index_temp, function(x) mean(sapply(1:nbr_it_bound, function(y) mean(PR_FDR[[x]]$empirical$PR_hat[,which(results_selection$Our_algorithm$K_hat[[x]][[y]]==cste_mult_vect)]))))
FDR_our_algorithm <- sapply(index_temp, function(x) mean(sapply(1:nbr_it_bound, function(y) mean(PR_FDR[[x]]$empirical$FDP_hat[,which(results_selection$Our_algorithm$K_hat[[x]][[y]]==cste_mult_vect)]))))


##### Generation of tables of results
if (ordered_variable == TRUE){
  matrix_selection_1 <- matrix(0, nrow = 4, ncol = 3)
  colnames(matrix_selection_1) <- c("m_Hat","PR(m_Hat)","FDR(m_Hat)")
  rownames(matrix_selection_1) <- c("LinSelect","KCV","Knockoff","Our_algorithm")
  matrix_selection_1[,1] <- c(mhat_LinSelect[1],mhat_KCV[1],mhat_Knockoffs[1],mhat_our_algorithm[1])
  matrix_selection_1[,2] <- c(PR_LinSelect[1],PR_KCV[1],PR_Knockoffs[1],PR_our_algorithm[1])
  matrix_selection_1[,3] <- c(FDR_LinSelect[1],FDR_KCV[1],FDR_Knockoffs[1],FDR_our_algorithm[1])
  xtable(matrix_selection_1, type = "latex", file = "matrice.tex",digits=2)
}

if (random_collection == TRUE){
  # collection Bolasso
  matrix_selection_1 <- matrix(0, nrow = 4, ncol = 3)
  colnames(matrix_selection_1) <- c("m_Hat","PR(m_Hat)","FDR(m_Hat)")
  rownames(matrix_selection_1) <- c("LinSelect","KCV","Knockoff","Our_algorithm")
  matrix_selection_1[,1] <- c(mhat_LinSelect[1],mhat_KCV[1],mhat_Knockoffs[1],mhat_our_algorithm[1])
  matrix_selection_1[,2] <- c(PR_LinSelect[1],PR_KCV[1],PR_Knockoffs[1],PR_our_algorithm[1])
  matrix_selection_1[,3] <- c(FDR_LinSelect[1],FDR_KCV[1],FDR_Knockoffs[1],FDR_our_algorithm[1])

  # collection Slope
  matrix_selection_2 <- matrix(0, nrow = 4, ncol = 3)
  colnames(matrix_selection_2) <- c("m_Hat","PR(m_Hat)","FDR(m_Hat)")
  rownames(matrix_selection_2) <- c("LinSelect","KCV","Knockoff","Our_algorithm")
  matrix_selection_2[,1] <- c(mhat_LinSelect[2],mhat_KCV[2],mhat_Knockoffs[2],mhat_our_algorithm[2])
  matrix_selection_2[,2] <- c(PR_LinSelect[2],PR_KCV[2],PR_Knockoffs[2],PR_our_algorithm[2])
  matrix_selection_2[,3] <- c(FDR_LinSelect[2],FDR_KCV[2],FDR_Knockoffs[2],FDR_our_algorithm[2])

  # collection RF
  matrix_selection_3 <- matrix(0, nrow = 4, ncol = 3)
  colnames(matrix_selection_3) <- c("m_Hat","PR(m_Hat)","FDR(m_Hat)")
  rownames(matrix_selection_3) <- c("LinSelect","KCV","Knockoff","Our_algorithm")
  matrix_selection_3[,1] <- c(mhat_LinSelect[3],mhat_KCV[3],mhat_Knockoffs[3],mhat_our_algorithm[3])
  matrix_selection_3[,2] <- c(PR_LinSelect[3],PR_KCV[3],PR_Knockoffs[3],PR_our_algorithm[3])
  matrix_selection_3[,3] <- c(FDR_LinSelect[3],FDR_KCV[3],FDR_Knockoffs[3],FDR_our_algorithm[3])

  # collection Knockoff
  matrix_selection_4 <- matrix(0, nrow = 4, ncol = 3)
  colnames(matrix_selection_4) <- c("m_Hat","PR(m_Hat)","FDR(m_Hat)")
  rownames(matrix_selection_4) <- c("LinSelect","KCV","Knockoff","Our_algorithm")
  matrix_selection_4[,1] <- c(mhat_LinSelect[4],mhat_KCV[4],mhat_Knockoffs[4],mhat_our_algorithm[4])
  matrix_selection_4[,2] <- c(PR_LinSelect[4],PR_KCV[4],PR_Knockoffs[4],PR_our_algorithm[4])
  matrix_selection_4[,3] <- c(FDR_LinSelect[4],FDR_KCV[4],FDR_Knockoffs[4],FDR_our_algorithm[4])

  xtable(matrix_selection_1, type = "latex", file = "matrice.tex",digits=2)
  xtable(matrix_selection_2, type = "latex", file = "matrice.tex",digits=2)
  xtable(matrix_selection_3, type = "latex", file = "matrice.tex",digits=2)
  xtable(matrix_selection_4, type = "latex", file = "matrice.tex",digits=2)
}