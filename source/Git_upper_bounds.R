# This function source is devoted to compute the theoretical upper bound on the FDR function with respect to the considered grid of K.
# The upper bounds of the P_1r terms are computed, both with unknown parameters (theoretical bound) and with estimated parameters (data-driven bounds). 
# For the last ones, the variance is either known or unknown (estimated by the slope heuristic principle).
# Then, the upper bound on the FDR is computed, both with unknown parameters (theoretical bound) and with estimated parameters (data-driven bounds).
upper_bounds <- function(detail,true_beta_matrix,n,p,cste_mult_vect,sigma_2,path_collection,path_collection_bounds,list_P_2r){
  q = min(n,p) 
  set.seed(1234)
  pb = txtProgressBar(min = 1, max = ncol(true_beta_matrix), initial = 0) 
  upper_bound_list <- list()
  upper_bound_list_hat_one_dataset_K2_list <- list()  # with \Tilde{K} = 2 with a known variance
  upper_bound_list_hat_one_dataset_K3_list <- list()  # with \Tilde{K} = 3 with a known variance
  upper_bound_list_hat_one_dataset_K2.5_list <- list()  # with \Tilde{K} = 2.5 with a known variance
  upper_bound_list_hat_one_dataset_K3.5_list <- list()  # with \Tilde{K} = 3.5 with a known variance
  upper_bound_list_hat_one_dataset_K1.5_list <- list()  # with \Tilde{K} = 1.5 with a known variance
  upper_bound_list_hat_one_dataset_K1_list <- list()  # with \Tilde{K} = 1 with a known variance
  upper_bound_list_hat_one_dataset_Klog_n_list <- list()  # with \Tilde{K} = log(n) with a known variance
  upper_bound_list_hat_one_dataset_K4_list <- list()  # with \Tilde{K} = 4 with a known variance
  upper_bound_list_hat_one_dataset_K4.5_list <- list()  # with \Tilde{K} = 4.5 with a known variance
  upper_bound_list_hat_one_dataset_K5_list <- list()  # with \Tilde{K} = 5 with a known variance
  upper_bound_list_hat_one_dataset_K2_list_slope <- list()  # with \Tilde{K} = 2 with an unknown variance
  upper_bound_list_hat_one_dataset_K3_list_slope <- list()  # with \Tilde{K} = 3 with an unknown variance
  upper_bound_list_hat_one_dataset_K2.5_list_slope <- list()  # with \Tilde{K} = 2.5 with an unknown variance
  upper_bound_list_hat_one_dataset_K3.5_list_slope <- list()  # with \Tilde{K} = 3.5 with an unknown variance
  upper_bound_list_hat_one_dataset_K1.5_list_slope <- list()  # with \Tilde{K} = 1.5 with an unknown variance
  upper_bound_list_hat_one_dataset_K1_list_slope <- list()  # with \Tilde{K} = 1 with an unknown variance
  upper_bound_list_hat_one_dataset_Klog_n_list_slope <- list()  # with \Tilde{K} = log(n) with an unknown variance
  upper_bound_list_hat_one_dataset_K4_list_slope <- list()  # with \Tilde{K} = 4 with an unknown variance
  upper_bound_list_hat_one_dataset_K4.5_list_slope <- list()  # with \Tilde{K} = 4.5 with an unknown variance
  upper_bound_list_hat_one_dataset_K5_list_slope <- list()  # with \Tilde{K} = 5 with an unknown variance
  for (k in 1:ncol(true_beta_matrix)){
    D_m_star <- length(which(true_beta_matrix[,k] !=0))
    upper_bound <- c()
    term_1_r_K <- list()
    upper_bound_hat_one_dataset_K2 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 2 with a known variance
    upper_bound_hat_one_dataset_K3 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 3 with a known variance
    upper_bound_hat_one_dataset_K2.5 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 2.5 with a known variance
    upper_bound_hat_one_dataset_K3.5 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 3.5 with a known variance
    upper_bound_hat_one_dataset_K1.5 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 1.5 with a known variance
    upper_bound_hat_one_dataset_K1 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 1 with a known variance
    upper_bound_hat_one_dataset_Klog_n <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = log(n) with a known variance
    upper_bound_hat_one_dataset_K4 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 4 with a known variance
    upper_bound_hat_one_dataset_K4.5 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 4.5 with a known variance
    upper_bound_hat_one_dataset_K5 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 5 with a known variance
    upper_bound_hat_one_dataset_K2_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 2 with an unknown variance
    upper_bound_hat_one_dataset_K3_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 3 with an unknown variance
    upper_bound_hat_one_dataset_K2.5_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 2.5 with an unknown variance
    upper_bound_hat_one_dataset_K3.5_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 3.5 with an unknown variance
    upper_bound_hat_one_dataset_K1.5_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 1.5 with an unknown variance
    upper_bound_hat_one_dataset_K1_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 1 with an unknown variance
    upper_bound_hat_one_dataset_Klog_n_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = log(n) with an unknown variance
    upper_bound_hat_one_dataset_K4_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 4 with an unknown variance
    upper_bound_hat_one_dataset_K4.5_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 4.5 with an unknown variance
    upper_bound_hat_one_dataset_K5_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 5 with an unknown variance
    ##### The theoretical upper bound
    for (K in cste_mult_vect){
      m = which(cste_mult_vect==K)
      temp_K <- c()
      temp_1_r_K <- c()
      if (D_m_star == q){  # no FP
        upper_bound <- 0
        term_1_r_K[[m]] <- 0
      }else{
        if (D_m_star < (q-1)){
          for (r in ((D_m_star+1):(q-1))){    # if r = q, the event is empty,      # Computation of the reccurence for each term of the sum
            temp <- (r-D_m_star)/r
            if (D_m_star == 0){
              max_without_beta_part <- max(sapply(1:(r-D_m_star), function(j) pchisq((j*K), j, ncp = 0)))
              temp_1_r <- 1-max_without_beta_part
            }else{
              max_without_beta_part <- max(sapply(1:(r-D_m_star), function(j) pchisq((j*K), j, ncp = 0)))
              max_with_beta_part <- max(sapply(0:(D_m_star-1), function(j) pchisq((((r-j)*K)/2)-sum(sapply((j+1):D_m_star, function(l) (as.numeric(true_beta_matrix[l,k])^2)/sigma_2)),r-j, ncp = 0)))
              temp_1_r <- 1-max(max_without_beta_part,max_with_beta_part)
            }
            temp_1_r_K <- c(temp_1_r_K,temp_1_r)
            temp_K <- c(temp_K,temp*temp_1_r*list_P_2r[[D_m_star+1]][[m]][r-D_m_star])
          }
        }
        if (D_m_star == 0){
          max_without_beta_part <- max(sapply(1:(q-D_m_star), function(j) pchisq((j*K), j, ncp = 0)))
          temp_1_r <- 1-max_without_beta_part
        }else{
          max_without_beta_part <- max(sapply(1:(q-D_m_star), function(j) pchisq((j*K), j, ncp = 0)))
          max_with_beta_part <- max(sapply(0:(D_m_star-1), function(j) pchisq((((q-j)*K)/2)-sum(sapply((j+1):D_m_star, function(l) (as.numeric(true_beta_matrix[l,k])^2)/sigma_2)),q-j, ncp = 0)))
          temp_1_r <- 1-max(max_without_beta_part,max_with_beta_part)
        }
        temp_1_r_K <- c(temp_1_r_K,temp_1_r)
        temp_K <- c(temp_K,((q-D_m_star)/q)*temp_1_r)  # the term for r=q is added,   if D_m_star == q-1, there is only the term for r=q
        upper_bound <- c(upper_bound,sum(temp_K))
        term_1_r_K[[m]] <- temp_1_r_K
      }
      ##### Evaluation of the theoretical upper bound with the available dataset and the used procedure only.
      for (it in 1:nbr_it_bound){
        # known variance
        upper_bound_hat_one_dataset_K2[[it]] <- c(upper_bound_hat_one_dataset_K2[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =2,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        upper_bound_hat_one_dataset_K3[[it]] <- c(upper_bound_hat_one_dataset_K3[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =3,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        upper_bound_hat_one_dataset_K2.5[[it]] <- c(upper_bound_hat_one_dataset_K2.5[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =2.5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        upper_bound_hat_one_dataset_K3.5[[it]] <- c(upper_bound_hat_one_dataset_K3.5[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =3.5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        upper_bound_hat_one_dataset_K1.5[[it]] <- c(upper_bound_hat_one_dataset_K1.5[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =1.5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        upper_bound_hat_one_dataset_K1[[it]] <- c(upper_bound_hat_one_dataset_K1[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =1,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        upper_bound_hat_one_dataset_Klog_n[[it]] <- c(upper_bound_hat_one_dataset_Klog_n[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =log(n),sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        upper_bound_hat_one_dataset_K4[[it]] <- c(upper_bound_hat_one_dataset_K4[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =4,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        upper_bound_hat_one_dataset_K4.5[[it]] <- c(upper_bound_hat_one_dataset_K4.5[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =4.5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        upper_bound_hat_one_dataset_K5[[it]] <- c(upper_bound_hat_one_dataset_K5[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        # unknown variance
        upper_bound_hat_one_dataset_K2_slope[[it]] <- c(upper_bound_hat_one_dataset_K2_slope[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =2,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        upper_bound_hat_one_dataset_K3_slope[[it]] <- c(upper_bound_hat_one_dataset_K3_slope[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =3,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        upper_bound_hat_one_dataset_K2.5_slope[[it]] <- c(upper_bound_hat_one_dataset_K2.5_slope[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =2.5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        upper_bound_hat_one_dataset_K3.5_slope[[it]] <- c(upper_bound_hat_one_dataset_K3.5_slope[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =3.5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        upper_bound_hat_one_dataset_K1.5_slope[[it]] <- c(upper_bound_hat_one_dataset_K1.5_slope[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =1.5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        upper_bound_hat_one_dataset_K1_slope[[it]] <- c(upper_bound_hat_one_dataset_K1_slope[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =1,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        upper_bound_hat_one_dataset_Klog_n_slope[[it]] <- c(upper_bound_hat_one_dataset_Klog_n_slope[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =log(n),sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        upper_bound_hat_one_dataset_K4_slope[[it]] <- c(upper_bound_hat_one_dataset_K4_slope[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =4,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        upper_bound_hat_one_dataset_K4.5_slope[[it]] <- c(upper_bound_hat_one_dataset_K4.5_slope[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =4.5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        upper_bound_hat_one_dataset_K5_slope[[it]] <- c(upper_bound_hat_one_dataset_K5_slope[[it]],upper_bound_estimated_one_dataset(cste_mult = K,path = path_collection_bounds[[k]],n,q,chosen_cste =5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
      }
    }
    #### For plotting the details of the theoretical upper bound. 
    if (detail == TRUE){
      ## We compare upper bound for i between 0 and (r-1)
      crit_test_inf_r <- list()  # We approach P(pour tout i<r, crit(m_r)<crit(m_i)) wich is a term in the sum of the FDR
      crit_test_sup_r <- list()  # We approach P(pour tout i>r, crit(m_r)<crit(m_i)) wich is a term in the sum of the FDR
      path <- path_collection[[k]]
      for (K in cste_mult_vect){
        m = which(cste_mult_vect==K)
        number_crit_r_inf_crit_i_inf_r <- c()
        number_crit_r_inf_crit_i_sup_r <- c()
        for (r in (D_m_star+1):q){
          crit_r <- sapply(1:nbr_it, function(it) path[[it]]$val1$LS[r+1] + K*sigma_2*(r/n))
          crit_i_inf_r <- list()
          crit_i_sup_r <- list()
          for (i in 0:(r-1)){
            crit_i_inf_r[[i+1]] <- sapply(1:nbr_it, function(it) path[[it]]$val1$LS[i+1] + K*sigma_2*(i/n))
          }
          for (i in (r+1):q){
            crit_i_sup_r[[i+1]] <- sapply(1:nbr_it, function(it) path[[it]]$val1$LS[i+1] + K*sigma_2*(i/n))
          }
          temp_crit_inf_r <- mean(sapply(1:nbr_it, function(it) all(sapply(0:(r-1), function(l) crit_r[it] < crit_i_inf_r[[l+1]][it]))))
          temp_crit_sup_r <- mean(sapply(1:nbr_it, function(it) all(sapply((r+1):q, function(l) crit_r[it] < crit_i_sup_r[[l+1]][it]))))
          number_crit_r_inf_crit_i_inf_r <- c(number_crit_r_inf_crit_i_inf_r,temp_crit_inf_r)
          number_crit_r_inf_crit_i_sup_r <- c(number_crit_r_inf_crit_i_sup_r,temp_crit_sup_r)
        }
        number_crit_r_inf_crit_i_sup_r[length(number_crit_r_inf_crit_i_sup_r)] <- 1 # We change the case for r=q
        crit_test_inf_r[[m]] <- number_crit_r_inf_crit_i_inf_r
        crit_test_sup_r[[m]] <- number_crit_r_inf_crit_i_sup_r
      }
      par(mfrow=c(1,1))
      for (r in (D_m_star+1):q){
        plot(cste_mult_vect, sapply(1:length(cste_mult_vect), function(K) crit_test_inf_r[[K]][r-D_m_star]), ylim = c(0,1), col = 4, type = 'p',
             xlab = "K",ylab = "estimation and bound of the [0,r-1]-FDR part")
        lines(cste_mult_vect, sapply(1:length(cste_mult_vect), function(K) term_1_r_K[[K]][r-D_m_star]), col = 2)
        legend("topright", legend = c(paste("|beta^*| = ",D_m_star),paste("r = ",r),"Proba of crit(m_r)<inf_{i in [0,r-1]}(crit(m_i))","theoretical_upper_bound_0_r-1"),
               lwd = c(0,0,1,1), col = c(3,3,4,2), lty = 1, cex=0.7)
      }
      for (r in (D_m_star+1):q){
        plot(cste_mult_vect, sapply(1:length(cste_mult_vect), function(K) crit_test_sup_r[[K]][r-D_m_star]), ylim = c(0,1), col = 4, type = 'p',
             xlab = "K",ylab = "estimation and bound of the [r+1,q]-FDR part")
        lines(cste_mult_vect, sapply(1:length(cste_mult_vect), function(K) list_P_2r[[D_m_star+1]][[K]][r-D_m_star]), col = 2)
        legend("bottomright", legend = c(paste("|beta^*| = ",D_m_star),paste("r = ",r),"Proba of crit(m_r)<inf_{i in [r+1,q]}(crit(m_i))","theoretical_upper_bound_r+1_q"),
               lwd = c(0,0,1,1), col = c(3,3,4,2), lty = 1, cex=0.7)
      }
    }
    ##############
    setTxtProgressBar(pb,k)
    upper_bound_list[[k]] <- upper_bound
    upper_bound_list_hat_one_dataset_K2_list[[k]] <- upper_bound_hat_one_dataset_K2
    upper_bound_list_hat_one_dataset_K3_list[[k]] <- upper_bound_hat_one_dataset_K3
    upper_bound_list_hat_one_dataset_K2.5_list[[k]] <- upper_bound_hat_one_dataset_K2.5
    upper_bound_list_hat_one_dataset_K3.5_list[[k]] <- upper_bound_hat_one_dataset_K3.5
    upper_bound_list_hat_one_dataset_K1.5_list[[k]] <- upper_bound_hat_one_dataset_K1.5
    upper_bound_list_hat_one_dataset_K1_list[[k]] <- upper_bound_hat_one_dataset_K1
    upper_bound_list_hat_one_dataset_Klog_n_list[[k]] <- upper_bound_hat_one_dataset_Klog_n
    upper_bound_list_hat_one_dataset_K4_list[[k]] <- upper_bound_hat_one_dataset_K4
    upper_bound_list_hat_one_dataset_K4.5_list[[k]] <- upper_bound_hat_one_dataset_K4.5
    upper_bound_list_hat_one_dataset_K5_list[[k]] <- upper_bound_hat_one_dataset_K5
    upper_bound_list_hat_one_dataset_K2_list_slope[[k]] <- upper_bound_hat_one_dataset_K2_slope
    upper_bound_list_hat_one_dataset_K3_list_slope[[k]] <- upper_bound_hat_one_dataset_K3_slope
    upper_bound_list_hat_one_dataset_K2.5_list_slope[[k]] <- upper_bound_hat_one_dataset_K2.5_slope
    upper_bound_list_hat_one_dataset_K3.5_list_slope[[k]] <- upper_bound_hat_one_dataset_K3.5_slope
    upper_bound_list_hat_one_dataset_K1.5_list_slope[[k]] <- upper_bound_hat_one_dataset_K1.5_slope
    upper_bound_list_hat_one_dataset_K1_list_slope[[k]] <- upper_bound_hat_one_dataset_K1_slope
    upper_bound_list_hat_one_dataset_Klog_n_list_slope[[k]] <- upper_bound_hat_one_dataset_Klog_n_slope
    upper_bound_list_hat_one_dataset_K4_list_slope[[k]] <- upper_bound_hat_one_dataset_K4_slope
    upper_bound_list_hat_one_dataset_K4.5_list_slope[[k]] <- upper_bound_hat_one_dataset_K4.5_slope
    upper_bound_list_hat_one_dataset_K5_list_slope[[k]] <- upper_bound_hat_one_dataset_K5_slope
    print(k)
  }
  return(data.frame(upper_bound_list=I(upper_bound_list),
                    upper_bound_list_hat_one_dataset_K2_list=I(upper_bound_list_hat_one_dataset_K2_list),upper_bound_list_hat_one_dataset_K3_list=I(upper_bound_list_hat_one_dataset_K3_list),
                    upper_bound_list_hat_one_dataset_K2.5_list=I(upper_bound_list_hat_one_dataset_K2.5_list),upper_bound_list_hat_one_dataset_K3.5_list=I(upper_bound_list_hat_one_dataset_K3.5_list),
                    upper_bound_list_hat_one_dataset_K1.5_list=I(upper_bound_list_hat_one_dataset_K1.5_list),upper_bound_list_hat_one_dataset_K1_list=I(upper_bound_list_hat_one_dataset_K1_list),
                    upper_bound_list_hat_one_dataset_Klog_n_list=I(upper_bound_list_hat_one_dataset_Klog_n_list),upper_bound_list_hat_one_dataset_K4_list = I(upper_bound_list_hat_one_dataset_K4_list),
                    upper_bound_list_hat_one_dataset_K4.5_list = I(upper_bound_list_hat_one_dataset_K4.5_list),upper_bound_list_hat_one_dataset_K5_list = I(upper_bound_list_hat_one_dataset_K5_list),
                    upper_bound_list_hat_one_dataset_K2_list_slope=I(upper_bound_list_hat_one_dataset_K2_list_slope),upper_bound_list_hat_one_dataset_K3_list_slope=I(upper_bound_list_hat_one_dataset_K3_list_slope),
                    upper_bound_list_hat_one_dataset_K2.5_list_slope=I(upper_bound_list_hat_one_dataset_K2.5_list_slope),upper_bound_list_hat_one_dataset_K3.5_list_slope=I(upper_bound_list_hat_one_dataset_K3.5_list_slope),
                    upper_bound_list_hat_one_dataset_K1.5_list_slope=I(upper_bound_list_hat_one_dataset_K1.5_list_slope),upper_bound_list_hat_one_dataset_K1_list_slope=I(upper_bound_list_hat_one_dataset_K1_list_slope),
                    upper_bound_list_hat_one_dataset_Klog_n_list_slope=I(upper_bound_list_hat_one_dataset_Klog_n_list_slope),upper_bound_list_hat_one_dataset_K4_list_slope = I(upper_bound_list_hat_one_dataset_K4_list_slope),
                    upper_bound_list_hat_one_dataset_K4.5_list_slope = I(upper_bound_list_hat_one_dataset_K4.5_list_slope),upper_bound_list_hat_one_dataset_K5_list_slope = I(upper_bound_list_hat_one_dataset_K5_list_slope)))  
}