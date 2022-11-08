## script to launch simulations
## author : perrine lacroix
## date : March, 10 2021

## This code is the main R script. 
## From the created data set, the PR and the FDR are firstly computed for each K and for each iteration. 
## Then, the P_{2r} term in the FDR expression is estimated by an empirical process.
## Lower bounds and upper bounds are computed, for each K and for all the \Tilde{K} considered in the article
## The rest of the R script stands for plottings. 


rm(list=objects())

prefix = "~/Documents/GitHub/Trade_off_FDR_PR"

# Packages loading
library(xtable)
library(ggplot2)
library(capushe)

## Loading of coded function
setwd(paste(prefix,"/source",sep =""))
source("Git_cste_mult_variation_PR_FDP.R")
source("Git_empirical_estimation.R")
source("Git_P_2r_estimation.R")
source("Git_lower_bounds.R")
source("Git_upper_bounds.R")
source("Git_lower_bound_estimated_one_dataset.R")
source("Git_upper_bound_estimated_one_dataset.R")

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

# Loading the generated data sets
config_data = paste0(paste("sigma2",sigma_2,"p",p,"n",n,"nbr_it",nbr_it,sep="_"),".RData")
load(paste("matrix_X",config_data,sep="_"))
config_data = paste0(paste("beta_true_0_20","sigma2",sigma_2,"p",p,"n",n,"nbr_it",nbr_it,sep="_"),".RData")
load(paste("vector_Y",config_data,sep="_"))
# load(paste("smaller_noise_vector_Y",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("larger_noise_closed_values_vector_Y",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article
load(paste("path_collection",config_data,sep="_"))
# load(paste("smaller_noise_path_collection",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("larger_noise_closed_values_path_collection",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article



#####################################################################################################################
###### Computation of the PR and FDP of the selected model for each K and for each iteration
#####################################################################################################################

#### The following line gives the PR and the FDP of the selected model for each considered K and on all the nbr_it data sets
#### PR_FDP_Hat_K is a matrix of size nbr_it x length(cste_mult_vect)
PR_FDP_hat_K <- empirical_estimation(cste_mult_vect,nbr_it,n,p,sigma_2,path_collection,Y_metric_matrix,data_select,data_metric)
setwd(paste(prefix,"/data",sep =""))
save(PR_FDP_hat_K, file = paste("PR_FDP_hat_K",config_data,sep="_"))
# save(PR_FDP_hat_K, file = paste("smaller_noise_PR_FDP_hat_K",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# save(PR_FDP_hat_K, file = paste("larger_noise_closed_values_PR_FDP_hat_K",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article
##########################################################


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
list_lower_bounds <- lower_bounds(detail,true_beta_matrix,n,p,cste_mult_vect,sigma_2,path_collection,list_P_2r)
setwd(paste(prefix,"/data",sep =""))
save(list_lower_bounds,file=paste("lower_bounds_all",config_data,sep="_"))
# save(list_lower_bounds,file=paste("smaller_noise_lower_bounds_all",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# save(list_lower_bounds,file=paste("larger_noise_closed_values_lower_bounds_all",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article

############# UPPER BOUNDS
setwd(paste(prefix,"/data",sep =""))
load(paste0(paste("estimation_of_P_2r_for_all_D_m_in_0_q","p",p,"n",n,"asympt_number",asympt_number,sep="_"),".RData"))
load(paste("lower_bounds_all",config_data,sep="_"))
# load(paste("smaller_noise_lower_bounds_all",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("larger_noise_closed_values_lower_bounds_all",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article
detail <- FALSE
comparison <- TRUE
with_hat <- TRUE
list_upper_bounds <- upper_bounds(detail,true_beta_matrix,n,p,cste_mult_vect,sigma_2,path_collection,list_P_2r)
setwd(paste(prefix,"/data",sep =""))
save(list_upper_bounds,file=paste("upper_bounds_all",config_data,sep="_"))
# save(list_upper_bounds,file=paste("smaller_noise_upper_bounds_all",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# save(list_upper_bounds,file=paste("larger_noise_closed_values_upper_bounds_all",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article



#####################################################################################################################
###### Plotting
#####################################################################################################################
setwd(paste(prefix,"/data",sep =""))
load(paste0(paste("estimation_of_P_2r_for_all_D_m_in_0_q","p",p,"n",n,"asympt_number",asympt_number,sep="_"),".RData"))
load(paste("lower_bounds_all",config_data,sep="_"))
load(paste("upper_bounds_all",config_data,sep="_"))
# load(paste("smaller_noise_lower_bounds_all",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("smaller_noise_upper_bounds_all",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("larger_noise_closed_values_lower_bounds_all",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article
# load(paste("larger_noise_closed_values_upper_bounds_all",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article


###### FDR, theoretical upper bound and theoretical lower bound functions with respect to K
description_parameters = paste0(paste("sigma2",sigma_2,"p",p,"n",n,sep="_"),".RData")
pente_parameters = paste("final_upper_and_lower_bound_of_FDR",description_parameters,sep="_")
# pente_parameters = paste("final_upper_and_lower_bound_of_FDR_smaller_noise",description_parameters,sep="_")    # the configuration 2 of the scenario (ii)  of the article
# pente_parameters = paste("final_upper_and_lower_bound_of_FDR_larger_noise",description_parameters,sep="_")    # the configuration 3 of the scenario (ii)  of the article
setwd(paste(prefix,"/results",sep =""))
pdf(paste0(pente_parameters,".pdf"))
for (k in 1:ncol(true_beta_matrix)){
  quantile_asympt_95 <- 1.96      # for the confidence interval of the FDR empirical estimation. Application of the central limit theorem
  vect_estimation_FDR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))   # empirical mean
  variance_estimation_FDR <- sapply(1:length(cste_mult_vect), function(it) (1/(nbr_it-1))*sum((PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]-mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))^2))    # empirical variance
  FDP_confidence_interval_inf <- vect_estimation_FDR - (variance_estimation_FDR*quantile_asympt_95)/sqrt(nbr_it)     # lower bound of the confidence interval for the FDR empirical estimation
  FDP_confidence_interval_sup <- vect_estimation_FDR + (variance_estimation_FDR*quantile_asympt_95)/sqrt(nbr_it)      # upper bound of the confidence interval for the FDR empirical estimation
  par(mfrow=c(1,1))
  plot(cste_mult_vect,FDP_confidence_interval_inf, col = 1,type = 'l',xlim = c(0.1,10),ylim=c(0,max(list_upper_bounds$upper_bound_list[[k]])),xlab = "K",ylab = "empirical FDR and bounds", lty = 5)
  lines(cste_mult_vect,FDP_confidence_interval_sup, col = 1, type = 'l', lty = 5)
  lines(cste_mult_vect,vect_estimation_FDR,col = 2, type="p", lty = 5)
  lines(cste_mult_vect,vect_estimation_FDR,col = 2,type="l")
  lines(cste_mult_vect,list_upper_bounds$upper_bound_list[[k]], col = 4, type="p")
  lines(cste_mult_vect,list_upper_bounds$upper_bound_list[[k]], col = 4, type="l")
  lines(cste_mult_vect,list_lower_bounds$lower_bound_list[[k]], col = 3, type = "p")
  lines(cste_mult_vect,list_lower_bounds$lower_bound_list[[k]], col = 3, type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","theoretical_upper_bound","theoretical_lower_bound","95% empirical confidence interval"),
         col = c(2,4,3,1), lty = 1, cex=1.3, lwd = c(1,1,1,1,1))
}

###### FDR, theoretical lower bound and estimated lower bounds functions with respect to K
which_cste_sup2 <- which(cste_mult_vect>=2)
sequence_temp <- seq(1,length(list_lower_bounds$lower_bound_list_hat_one_dataset_K2_list_slope[[1]][[1]]), by=2)
for (k in 1:ncol(true_beta_matrix)){
  vect_estimation_FDR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))   # empirical mean
  par(mfrow=c(2,3))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and lower bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K1_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K1_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","lower_bound_with_beta^*",
                                "Hat_lower_bound(beta(Hat(m)(1)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and lower bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K1.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K1.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","lower_bound_with_beta^*",
                                "Hat_lower_bound(beta(Hat(m)(1.5)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and lower bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K2_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K2_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","lower_bound_with_beta^*",
                                "Hat_lower_bound(beta(Hat(m)(2)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and lower bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K2.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K2.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","lower_bound_with_beta^*",
                                "Hat_lower_bound(beta(Hat(m)(2.5)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and lower bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K3_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K3_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","lower_bound_with_beta^*",
                                "Hat_lower_bound(beta(Hat(m)(3)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and lower bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K3.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K3.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","lower_bound_with_beta^*",
                                "Hat_lower_bound(beta(Hat(m)(3.5)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  par(mfrow=c(2,3))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and lower bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_Klog_n_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_Klog_n_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","lower_bound_with_beta^*",
                                "Hat_lower_bound(beta(Hat(m)(log(n))))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and lower bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","lower_bound_with_beta^*",
                                "Hat_lower_bound(beta(Hat(m)(4)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and lower bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K4.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K4.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","lower_bound_with_beta^*",
                                "Hat_lower_bound(beta(Hat(m)(4.5)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and lower bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_lower_bounds$lower_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","lower_bound_with_beta^*",
                                "Hat_lower_bound(beta(Hat(m)(5)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot.new()
  plot(cste_mult_vect,vect_estimation_FDR,col = 2,type="p",xlim = c(0.1,10),ylim=c(0,1),xlab = "K",ylab = "empirical FDR and lower bounds", lty = 5)
  lines(cste_mult_vect,vect_estimation_FDR,col = 2,type="l")
  lines(cste_mult_vect,list_lower_bounds$lower_bound_list[[k]], col = 4, type="p")
  lines(cste_mult_vect,list_lower_bounds$lower_bound_list[[k]], col = 4, type="l")
  lines(cste_mult_vect,as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]), col = "#993300", type = "p")
  lines(cste_mult_vect,as.numeric(list_lower_bounds$lower_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]), col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","lower_bound_with_beta^*",
                                "Hat_lower_bound(beta(Hat(m)(4)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
}

###### FDR, theoretical upper bound and estimated upper bounds functions with respect to K
which_cste_sup2 <- which(cste_mult_vect>=2)
for (k in 1:ncol(true_beta_matrix)){
  vect_estimation_FDR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))   # empirical mean
  par(mfrow=c(2,3))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and upper bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K1_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K1_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","upper_bound_with_beta^*",
                                "Hat_upper_bound(beta(Hat(m)(1)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and upper bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K1.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K1.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","upper_bound_with_beta^*",
                                "Hat_upper_bound(beta(Hat(m)(1.5)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and upper bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K2_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K2_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","upper_bound_with_beta^*",
                                "Hat_upper_bound(beta(Hat(m)(2)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and upper bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K2.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K2.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","upper_bound_with_beta^*",
                                "Hat_upper_bound(beta(Hat(m)(2.5)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and upper bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K3_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K3_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","upper_bound_with_beta^*",
                                "Hat_upper_bound(beta(Hat(m)(3)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and upper bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K3.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K3.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","upper_bound_with_beta^*",
                                "Hat_upper_bound(beta(Hat(m)(3.5)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  par(mfrow=c(2,3))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and upper bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_Klog_n_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_Klog_n_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","upper_bound_with_beta^*",
                                "Hat_upper_bound(beta(Hat(m)(log(n))))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and upper bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","upper_bound_with_beta^*",
                                "Hat_upper_bound(beta(Hat(m)(4)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and upper bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4.5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","upper_bound_with_beta^*",
                                "Hat_upper_bound(beta(Hat(m)(4.5)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="p",xlim = c(2,10),xlab = "K",ylab = "empirical FDR and upper bounds", lty = 5)
  lines(cste_mult_vect[which_cste_sup2],vect_estimation_FDR[which_cste_sup2],col = 2,type="l")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="p")
  lines(cste_mult_vect[which_cste_sup2],list_upper_bounds$upper_bound_list[[k]][which_cste_sup2], col = 4, type="l")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "p")
  lines(cste_mult_vect[which_cste_sup2],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K5_list_slope[[k]][[1]][sequence_temp])[which_cste_sup2], col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","upper_bound_with_beta^*",
                                "Hat_upper_bound(beta(Hat(m)(5)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
  plot.new()
  plot(cste_mult_vect,vect_estimation_FDR,col = 2,type="p",xlim = c(0.1,10),ylim=c(0,max(list_upper_bounds$upper_bound_list[[k]],as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]),vect_estimation_FDR)),xlab = "K",ylab = "empirical FDR and upper bounds", lty = 5)
  lines(cste_mult_vect,vect_estimation_FDR,col = 2,type="l")
  lines(cste_mult_vect,list_upper_bounds$upper_bound_list[[k]], col = 4, type="p")
  lines(cste_mult_vect,list_upper_bounds$upper_bound_list[[k]], col = 4, type="l")
  lines(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]), col = "#993300", type = "p")
  lines(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]), col = "#993300", type = "l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","upper_bound_with_beta^*",
                                "Hat_upper_bound(beta(Hat(m)(4)))"),
         col = c(2,4,"#993300"), lty = 1,cex=1.4, lwd = c(1,1,1))
}

###### Relative change of the estimated lower bounds functions with respect to K
for (k in 1:ncol(true_beta_matrix)){
  par(mfrow=c(2,3))
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
  plot(cste_mult_vect,relative_change_K1,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for lower bounds",type="p")
  lines(cste_mult_vect,relative_change_K1,col = 1, type="l")
  legend("topright", legend = c("Hat_lower_bound(beta(Hat(m)(1)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_change_K1.5,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for lower bounds",type="p")
  lines(cste_mult_vect,relative_change_K1.5,col = 1, type="l")
  legend("topleft", legend = c("Hat_lower_bound(beta(Hat(m)(1.5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_change_K2,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for lower bounds",type="p")
  lines(cste_mult_vect,relative_change_K2,col = 1, type="l")
  legend("topleft", legend = c("Hat_lower_bound(beta(Hat(m)(2)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_change_K2.5,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for lower bounds",type="p")
  lines(cste_mult_vect,relative_change_K2.5,col = 1, type="l")
  legend("topleft", legend = c("Hat_lower_bound(beta(Hat(m)(2.5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_change_K3,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for lower bounds",type="p")
  lines(cste_mult_vect,relative_change_K3,col = 1, type="l")
  legend("topleft", legend = c("Hat_lower_bound(beta(Hat(m)(3)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_change_K3.5,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for lower bounds",type="p")
  lines(cste_mult_vect,relative_change_K3.5,col = 1, type="l")
  legend("topleft", legend = c("Hat_lower_bound(beta(Hat(m)(3.5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  par(mfrow=c(2,3))
  plot(cste_mult_vect,relative_change_Klog_n,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for lower bounds",type="p")
  lines(cste_mult_vect,relative_change_Klog_n,col = 1, type="l")
  legend("topleft", legend = c("Hat_lower_bound(beta(Hat(m)(log(n))))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_change_K4,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for lower bounds",type="p")
  lines(cste_mult_vect,relative_change_K4,col = 1, type="l")
  legend("topleft", legend = c("Hat_lower_bound(beta(Hat(m)(4)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_change_K4.5,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for lower bounds",type="p")
  lines(cste_mult_vect,relative_change_K4.5,col = 1, type="l")
  legend("topleft", legend = c("Hat_lower_bound(beta(Hat(m)(4.5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_change_K5,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for lower bounds",type="p")
  lines(cste_mult_vect,relative_change_K5,col = 1, type="l")
  legend("topleft", legend = c("Hat_lower_bound(beta(Hat(m)(5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
}

###### Relative standard error of the estimated lower bounds functions with respect to K
for (k in 1:ncol(true_beta_matrix)){
  ###### Variance of estimators
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
  par(mfrow=c(2,3))
  plot(cste_mult_vect,relative_standard_error_K1,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K1,col = 1, type="l")
  legend("bottomright", legend = c("Hat_lower_bound(beta(Hat(m)(1)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_standard_error_K1.5,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K1.5,col = 1, type="l")
  legend("topright", legend = c("Hat_lower_bound(beta(Hat(m)(1.5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_standard_error_K2,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K2,col = 1, type="l")
  legend("topright", legend = c("Hat_lower_bound(beta(Hat(m)(2)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_standard_error_K2.5,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K2.5,col = 1, type="l")
  legend("topright", legend = c("Hat_lower_bound(beta(Hat(m)(2.5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_standard_error_K3,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K3,col = 1, type="l")
  legend("topright", legend = c("Hat_lower_bound(beta(Hat(m)(3)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_standard_error_K3.5,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K3.5,col = 1, type="l")
  legend("topright", legend = c("Hat_lower_bound(beta(Hat(m)(3.5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  par(mfrow=c(2,3))
  plot(cste_mult_vect,relative_standard_error_Klog_n,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_Klog_n,col = 1, type="l")
  legend("topright", legend = c("Hat_lower_bound(beta(Hat(m)(log(n))))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_standard_error_K4,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K4,col = 1, type="l")
  legend("topright", legend = c("Hat_lower_bound(beta(Hat(m)(4)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_standard_error_K4.5,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K4.5,col = 1, type="l")
  legend("topright", legend = c("Hat_lower_bound(beta(Hat(m)(4.5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_standard_error_K5,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated lower bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K5,col = 1, type="l")
  legend("topright", legend = c("Hat_lower_bound(beta(Hat(m)(5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
}

###### Relative change of the estimated upper bounds functions with respect to K
for (k in 1:ncol(true_beta_matrix)){
  par(mfrow=c(2,3))
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
  plot(cste_mult_vect,relative_change_K1,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for upper bounds",type="p")
  lines(cste_mult_vect,relative_change_K1,col = 1, type="l")
  legend("topright", legend = c("Hat_upper_bound(beta(Hat(m)(1)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_change_K1.5,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for upper bounds",type="p")
  lines(cste_mult_vect,relative_change_K1.5,col = 1, type="l")
  legend("bottomright", legend = c("Hat_upper_bound(beta(Hat(m)(1.5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_change_K2,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for upper bounds",type="p")
  lines(cste_mult_vect,relative_change_K2,col = 1, type="l")
  legend("bottomright", legend = c("Hat_upper_bound(beta(Hat(m)(2)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_change_K2.5,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for upper bounds",type="p")
  lines(cste_mult_vect,relative_change_K2.5,col = 1, type="l")
  legend("bottomright", legend = c("Hat_upper_bound(beta(Hat(m)(2.5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_change_K3,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for upper bounds",type="p")
  lines(cste_mult_vect,relative_change_K3,col = 1, type="l")
  legend("bottomright", legend = c("Hat_upper_bound(beta(Hat(m)(3)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_change_K3.5,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for upper bounds",type="p")
  lines(cste_mult_vect,relative_change_K3.5,col = 1, type="l")
  legend("bottomright", legend = c("Hat_upper_bound(beta(Hat(m)(3.5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  par(mfrow=c(2,3))
  plot(cste_mult_vect,relative_change_Klog_n,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for upper bounds",type="p")
  lines(cste_mult_vect,relative_change_Klog_n,col = 1, type="l")
  legend("bottomright", legend = c("Hat_upper_bound(beta(Hat(m)(log(n))))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_change_K4,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for upper bounds",type="p")
  lines(cste_mult_vect,relative_change_K4,col = 1, type="l")
  legend("bottomright", legend = c("Hat_upper_bound(beta(Hat(m)(4)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_change_K4.5,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for upper bounds",type="p")
  lines(cste_mult_vect,relative_change_K4.5,col = 1, type="l")
  legend("bottomright", legend = c("Hat_upper_bound(beta(Hat(m)(4.5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_change_K5,col = 1, xlim = c(0,10),xlab = "K",ylab = "relative difference for upper bounds",type="p")
  lines(cste_mult_vect,relative_change_K5,col = 1, type="l")
  legend("bottomright", legend = c("Hat_upper_bound(beta(Hat(m)(5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
}

###### Relative standard error of the estimated upper bounds functions with respect to K
for (k in 1:ncol(true_beta_matrix)){
###### Variance of estimators
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
  par(mfrow=c(2,3))
  plot(cste_mult_vect,relative_standard_error_K1,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K1,col = 1, type="l")
  legend("bottomright", legend = c("Hat_upper_bound(beta(Hat(m)(1)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_standard_error_K1.5,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K1.5,col = 1, type="l")
  legend("topright", legend = c("Hat_upper_bound(beta(Hat(m)(1.5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_standard_error_K2,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K2,col = 1, type="l")
  legend("topright", legend = c("Hat_upper_bound(beta(Hat(m)(2)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_standard_error_K2.5,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K2.5,col = 1, type="l")
  legend("topright", legend = c("Hat_upper_bound(beta(Hat(m)(2.5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_standard_error_K3,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K3,col = 1, type="l")
  legend("topright", legend = c("Hat_upper_bound(beta(Hat(m)(3)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_standard_error_K3.5,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K3.5,col = 1, type="l")
  legend("topright", legend = c("Hat_upper_bound(beta(Hat(m)(3.5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  par(mfrow=c(2,3))
  plot(cste_mult_vect,relative_standard_error_Klog_n,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_Klog_n,col = 1, type="l")
  legend("topright", legend = c("Hat_upper_bound(beta(Hat(m)(log(n))))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_standard_error_K4,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K4,col = 1, type="l")
  legend("topright", legend = c("Hat_upper_bound(beta(Hat(m)(4)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_standard_error_K4.5,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K4.5,col = 1, type="l")
  legend("topright", legend = c("Hat_upper_bound(beta(Hat(m)(4.5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
  plot(cste_mult_vect,relative_standard_error_K5,col = 1, xlim = c(0,10),ylim = c(0,ymax),xlab = "K",ylab = "Relative standard deviation of the estimated upper bounds",type="p")
  lines(cste_mult_vect,relative_standard_error_K5,col = 1, type="l")
  legend("topright", legend = c("Hat_upper_bound(beta(Hat(m)(5)))"),
         col = c(1), lty = 1,cex=1.4, lwd = c(1))
}

###### Empirical FDR and PR
for (k in 1:dim(true_beta_matrix)[2]){
  vect_estimation_FDR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))
  vect_estimation_PR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$PR_hat[,it]))
  par(mfrow=c(1,2))
  plot(cste_mult_vect,vect_estimation_FDR, col = 2,xlab = "K",ylim=c(0,1),ylab = "empirical FDR",type= "p", lty = 5)
  lines(cste_mult_vect,vect_estimation_FDR, col = 2,type="l")
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))"), col = c(2), lty = 1, cex=0.8, lwd = c(1))
  plot(cste_mult_vect,vect_estimation_PR, col = 4, type= "p", lty = 5,xlab = "K",ylab = "empirical PR")
  lines(cste_mult_vect,vect_estimation_PR, col = 4,type="l")
  legend("topright", legend = c("K -> empirical PR(Hat(m)(K))"), col = c(4), lty = 1, cex=0.8, lwd = c(1))
}

###### FDR and PR on one dataset
for (k in 1:dim(true_beta_matrix)[2]){
  vect_estimation_diff_PR <- PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope
  par(mfrow=c(1,2))
  plot(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]),xlab = "K",ylab = "estimated upper bound of FDR on one dataset",type="p", col = 6)
  lines(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]),type="l", col = 6)
  legend("topright", legend = c("Hat_upper_bound(beta(Hat(m)(4)))"), col = c(6), lty = 1, cex=0.5, lwd = c(1))
  plot(cste_mult_vect,vect_estimation_diff_PR, col = 4,xlab = "K",ylab = "estimated PR on one dataset", type="p")
  lines(cste_mult_vect,vect_estimation_diff_PR, col = 4,type="l")
  legend("topright", legend = c("K -> estimated PR(Hat(m)(K))"), col = c(4), lty = 1, cex=0.5, lwd = c(1))
}

########### \Hat{\kappa}: estimation of \sigma^2 by the slope heuristic with the capushe R package
for (k in 1:dim(true_beta_matrix)[2]){
  par(mfrow=c(1,1))
  vect_kappa <- sapply(1:nbr_it_bound, function(it) list_upper_bounds$upper_bound_list_hat_one_dataset_K3_list_slope[[k]][[it]]$kappa_slope)
  hist(vect_kappa, freq = NULL, main = NULL, ylab = "Frequency", xlab = "Hat(sigma^2)")
}

########### FDR and PR on one dataset with empirical FDR and PR for \Tilde{K} = 4
cste_for_estimation <- 4
for (k in 1:ncol(true_beta_matrix)){
  vect_estimation_FDR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,it]))
  vect_estimation_PR <- sapply(1:length(cste_mult_vect), function(it) mean(PR_FDP_hat_K[[k]]$empirical$PR_hat[,it]))
  par(mfrow=c(1,2))
  plot(cste_mult_vect,vect_estimation_FDR, col = 2,xlab = "K",ylim=c(0,max(as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]))),ylab = "empirical FDR and estimated upper bound of FDR ",type= "p", lty = 5)
  lines(cste_mult_vect,vect_estimation_FDR, col = 2,type="l")
  lines(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]),ylim=c(0,1),xlab = "K",ylab = "estimated upper bound of FDR on one dataset",type="p", col = 6)
  lines(cste_mult_vect,as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[1]][sequence_temp]),type="l", col = 6)
  legend("topright", legend = c("K -> empirical FDR(Hat(m)(K))","Hat_upper_bound(beta(Hat(m)(4)))"), col = c(2,6), lty = 1, cex=0.7, lwd = c(1,1))
  min_diff_PR <- min(vect_estimation_PR,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope)
  max_diff_PR <- max(vect_estimation_PR,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope)
  plot(cste_mult_vect,vect_estimation_PR, col = 4, type="p", lty = 5,ylim = c(min_diff_PR,max_diff_PR),xlab = "K",ylab = "empirical PR and estimated PR")
  lines(cste_mult_vect,vect_estimation_PR, col = 4,type="l")
  lines(cste_mult_vect,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope, col = "#660099", type="p")
  lines(cste_mult_vect,PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope, col = "#660099", type="l")
  legend("topright", legend = c("K -> empirical PR(Hat(m)(K))","K -> estimated PR(Hat(m)(K))"), col = c(4,"#660099"), lty = 1, cex=0.7, lwd = c(1,1))
}
dev.off()
