# This function source applies the R function source "cste_mult_variation_LS_PR_FDR" on each available model collection 
# when variables are non ordered and when D_m^* = 10. 
empirical_estimation_non_ordered <- function(cste_mult_vect,nbr_it,n,p,sigma_2,path_collection,Y_metric_matrix,data_select,data_metric){
  PR_FDP_hat_K <- list()
  time.in = Sys.time() 
  pb = txtProgressBar(min = 1, max = length(path_collection), initial = 1) 
  k = 11     # D_m^* = 10
  Y_metric <- Y_metric_matrix[,k]
  Y_select <- Y_select_matrix[,k]
  beta_true <- paste('beta_true',k-1,sep = "_")
  size_beta_true <- k-1
  for(l in 1:length(path_collection)){
    PR_FDP_hat_K[[l]] <- list()
    path <- path_collection[[l]]
    PR_FDP_hat_K[[l]] <- cste_mult_variation_LS_PR_FDR(cste_mult_vect,nbr_it,path,data_select,data_metric,Y_metric,size_beta_true,n,p,sigma_2)
    setTxtProgressBar(pb,l)
  }
  return(PR_FDP_hat_K)
}