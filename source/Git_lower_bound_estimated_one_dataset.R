# This function source is devoted to evaluate the theoretical lower bound with the available dataset and the used procedure only.
# If the variance is unknown, it is estimated by the slope heuristics principle with R package capushe. 
# The substitution of beta^* is given as an entry of the function. 
lower_bound_estimated_one_dataset <- function(cste_mult,path,n,q,chosen_cste,sigma_2,m,cste_mult_vect,list_P_2r,known_variance){
  set.seed(1234)
  cste_K <- which(cste_mult_vect==cste_mult)
  if (known_variance == FALSE){   # case of unknown variance: it is estimated by the slope heuristics principle with R package capushe
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
    lower_bound_hat <- 0
  }else{
    for (r in ((D_m+1):q)){     # Computation of the reccurence for each term of the sum
      temp <- (r-D_m)/r
      A <- c()
      B <- c()
      if (D_m == 0){
        for (j in (D_m+1):r){
          A <- c(A,2*(1-pnorm(q = sqrt(j*cste_mult), mean = 0, sd = 1)))
        }
        if (r>=2){
          for (j in (D_m+2):r){
            B <- c(B,2*(pnorm(q = sqrt(j*cste_mult), mean = 0, sd = 1)-pnorm(q = sqrt(cste_mult), mean = 0, sd = 1)))
          }
        }
      }else{
        #A
        for (j in (1:D_m)){
          A <- c(A,2-(pnorm(q = sqrt(j*cste_mult)-beta_m[j]/sqrt(variance), mean = 0, sd = 1)+pnorm(q = sqrt(j*cste_mult)+beta_m[j]/sqrt(variance), mean = 0, sd = 1)))
        }
        if (D_m>=2){
          #B
          for (j in 2:D_m){
            B <- c(B,pnorm(q = sqrt(j*cste_mult)-beta_m[j]/sqrt(variance), mean = 0, sd = 1)+pnorm(q = sqrt(j*cste_mult)+beta_m[j]/sqrt(variance), mean = 0, sd = 1)
                   -(pnorm(q = sqrt(cste_mult)-beta_m[j]/sqrt(variance), mean = 0, sd = 1)+pnorm(q = sqrt(cste_mult)+beta_m[j]/sqrt(variance), mean = 0, sd = 1)))
          }
        }
        for (j in (D_m+1):r){
          A <- c(A,2*(1-pnorm(q = sqrt(j*cste_mult), mean = 0, sd = 1)))
          B <- c(B,2*(pnorm(q = sqrt(j*cste_mult), mean = 0, sd = 1)-pnorm(q = sqrt(cste_mult), mean = 0, sd = 1)))
        }
      }
      temp_A_B <- A[1]
      if (length(B) > 0){    # If r !=1 (only possible if D_m_star =0)
        for (l in 1:length(B)){
          temp_A_B <- A[l+1]+B[l]*temp_A_B
        }
      }
      temp_1_r <- temp_A_B
      temp_K_hat <- c(temp_K_hat,temp*temp_1_r*list_P_2r[[D_m+1]][[cste_K]][r-D_m])
    }
    lower_bound_hat <- sum(temp_K_hat)
  }
  if (known_variance == FALSE){
    return(list(lower_bound_hat = I(lower_bound_hat), kappa_slope = I(variance)))
  }else{
    return(lower_bound_hat)
  }
}