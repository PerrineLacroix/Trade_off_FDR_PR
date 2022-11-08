# This function source is devoted to evaluate the theoretical upper bound with the available dataset and the used procedure only.
# If the variance is unknown, it is estimated by the slope heuristics principle with R package capushe. 
# The substitution of beta^* is given as an entry of the function. 
upper_bound_estimated_one_dataset <- function(cste_mult,path,n,q,chosen_cste,sigma_2,m,cste_mult_vect,list_P_2r,known_variance){
  set.seed(1234)
  cste_K <- which(cste_mult_vect==cste_mult)
  if (known_variance == FALSE){     # case of unknown variance: it is estimated by the slope heuristics principle with R package capushe
    forcap <- cbind(path[[m]]$val1$dim,path[[m]]$val1$dim/n,path[[m]]$val1$dim,path[[m]]$val1$LS)
    ResCapushe <- capushe(data = forcap)
    variance <- as.numeric(ResCapushe@DDSE@graph$reg$coefficients[2])
  }else{
    variance <- sigma_2
  }
  model_hat <- which.min(chosen_cste*variance*(path[[m]]$val1$dim/n) + path[[m]]$val1$LS)
  D_m <- path[[m]]$val1$dim[model_hat]
  beta_m <- as.numeric(path[[m]]$val2[model_hat,])
  temp <- 0
  temp_K_hat <- 0
  if (D_m == q){   # no FP
    upper_bound_hat <- 0
  }else{
    if (D_m < (q-1)){
      for (r in ((D_m+1):(q-1))){    # if r = q, the event is empty,      # Computation of the reccurence for each term of the sum
        temp <- (r-D_m)/r
        if (D_m == 0){
          max_without_beta_part <- max(sapply(1:(r-D_m), function(j) pchisq((j*cste_mult), j, ncp = 0)))
          temp_1_r <- 1-max_without_beta_part
        }else{
          max_without_beta_part <- max(sapply(1:(r-D_m), function(j) pchisq((j*cste_mult), j, ncp = 0)))
          max_with_beta_part <- max(sapply(0:(D_m-1), function(j) pchisq((((r-j)*cste_mult)/2)-sum(sapply((j+1):D_m, function(l) (beta_m[l]^2)/variance)),r-j, ncp = 0)))
          temp_1_r <- 1-max(max_without_beta_part,max_with_beta_part)
        }
        temp_K_hat <- c(temp_K_hat,temp*temp_1_r*list_P_2r[[D_m+1]][[cste_K]][r-D_m])
      }
    }
    if (D_m == 0){
      max_without_beta_part <- max(sapply(1:(q-D_m), function(j) pchisq((j*cste_mult), j, ncp = 0)))
      temp_1_r <- 1-max_without_beta_part
    }else{
      max_without_beta_part <- max(sapply(1:(q-D_m), function(j) pchisq((j*cste_mult), j, ncp = 0)))
      max_with_beta_part <- max(sapply(0:(D_m-1), function(j) pchisq((((q-j)*cste_mult)/2)-sum(sapply((j+1):D_m, function(l) (beta_m[l]^2)/variance)),q-j, ncp = 0)))
      temp_1_r <- 1-max(max_without_beta_part,max_with_beta_part)
    }
    temp_K_hat <- c(temp_K_hat,((q-D_m)/q)*temp_1_r)   # the term for r=q is added,   if D_m_star == q-1, there is only the term for r=q
    upper_bound_hat <- sum(temp_K_hat)
  }
  if (known_variance == FALSE){
    return(list(upper_bound_hat = I(upper_bound_hat), kappa_slope = I(variance)))
  }else{
    return(upper_bound_hat)
  }
}