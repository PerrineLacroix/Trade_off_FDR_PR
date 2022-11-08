# This function source is devoted to compute the theoretical lower bound on the FDR function with respect to the considered grid of K.
# The lower bounds of the P_1r terms are computed, both with unknown parameters (theoretical bound) and with estimated parameters (data-driven bounds).
# For the last ones, the variance is either known or unknown (estimated by the slope heuristic principle).
# Then, the lower bound on the FDR is computed, both with unknown parameters (theoretical bound) and with estimated parameters (data-driven bounds).
lower_bounds <- function(detail,true_beta_matrix,n,p,cste_mult_vect,sigma_2,path_collection,list_P_2r){
  q = min(n,p)
  set.seed(1234)
  pb = txtProgressBar(min = 1, max = ncol(true_beta_matrix), initial = 0) 
  lower_bound_list <- list()
  lower_bound_list_hat_one_dataset_K2_list <- list()  # with \Tilde{K} = 2 with a known variance
  lower_bound_list_hat_one_dataset_K3_list <- list()  # with \Tilde{K} = 3 with a known variance
  lower_bound_list_hat_one_dataset_K2.5_list <- list()  # with \Tilde{K} = 2.5 with a known variance
  lower_bound_list_hat_one_dataset_K3.5_list <- list()  # with \Tilde{K} = 3.5 with a known variance
  lower_bound_list_hat_one_dataset_K1.5_list <- list()  # with \Tilde{K} = 1.5 with a known variance
  lower_bound_list_hat_one_dataset_K1_list <- list()  # with \Tilde{K} = 1 with a known variance
  lower_bound_list_hat_one_dataset_Klog_n_list <- list()  # with \Tilde{K} = log(n) with a known variance
  lower_bound_list_hat_one_dataset_K4_list <- list()  # with \Tilde{K} = 4 with a known variance
  lower_bound_list_hat_one_dataset_K4.5_list <- list()  # with \Tilde{K} = 4.5 with a known variance
  lower_bound_list_hat_one_dataset_K5_list <- list()  # with \Tilde{K} = 5 with a known variance
  lower_bound_list_hat_one_dataset_K2_list_slope <- list()  # with \Tilde{K} = 2 with an unknown variance
  lower_bound_list_hat_one_dataset_K3_list_slope <- list()  # with \Tilde{K} = 3 with an unknown variance
  lower_bound_list_hat_one_dataset_K2.5_list_slope <- list()  # with \Tilde{K} = 2.5 with an unknown variance
  lower_bound_list_hat_one_dataset_K3.5_list_slope <- list()  # with \Tilde{K} = 3.5 with an unknown variance
  lower_bound_list_hat_one_dataset_K1.5_list_slope <- list()  # with \Tilde{K} = 1.5 with an unknown variance
  lower_bound_list_hat_one_dataset_K1_list_slope <- list()  # with \Tilde{K} = 1 with an unknown variance
  lower_bound_list_hat_one_dataset_Klog_n_list_slope <- list()  # with \Tilde{K} = log(n) with an unknown variance
  lower_bound_list_hat_one_dataset_K4_list_slope <- list()  # with \Tilde{K} = 4 with an unknown variance
  lower_bound_list_hat_one_dataset_K4.5_list_slope <- list()  # with \Tilde{K} = 4.5 with an unknown variance
  lower_bound_list_hat_one_dataset_K5_list_slope <- list()  # with \Tilde{K} = 5 with an unknown variance
  for (k in 1:ncol(true_beta_matrix)){
    D_m_star <- length(which(true_beta_matrix[,k] !=0))
    lower_bound <- c()
    term_1_r_K <- list()
    lower_bound_hat_one_dataset_K2 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 2 with a known variance
    lower_bound_hat_one_dataset_K3 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 3 with a known variance
    lower_bound_hat_one_dataset_K2.5 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 2.5 with a known variance
    lower_bound_hat_one_dataset_K3.5 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 3.5 with a known variance
    lower_bound_hat_one_dataset_K1.5 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 1.5 with a known variance
    lower_bound_hat_one_dataset_K1 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 1 with a known variance
    lower_bound_hat_one_dataset_Klog_n <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = log(n) with a known variance
    lower_bound_hat_one_dataset_K4 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 4 with a known variance
    lower_bound_hat_one_dataset_K4.5 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 4.5 with a known variance
    lower_bound_hat_one_dataset_K5 <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 5 with a known variance
    lower_bound_hat_one_dataset_K2_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 2 with an unknown variance
    lower_bound_hat_one_dataset_K3_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 3 with an unknown variance
    lower_bound_hat_one_dataset_K2.5_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 2.5 with an unknown variance
    lower_bound_hat_one_dataset_K3.5_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 3.5 with an unknown variance
    lower_bound_hat_one_dataset_K1.5_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 1.5 with an unknown variance
    lower_bound_hat_one_dataset_K1_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 1 with an unknown variance
    lower_bound_hat_one_dataset_Klog_n_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = log(n) with an unknown variance
    lower_bound_hat_one_dataset_K4_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 4 with an unknown variance
    lower_bound_hat_one_dataset_K4.5_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 4.5 with an unknown variance
    lower_bound_hat_one_dataset_K5_slope <- lapply(1:nbr_it_bound, function(it) c())  # with \Tilde{K} = 5 with an unknown variance
    ##### The theoretical lower bound
    for (K in cste_mult_vect){
      m = which(cste_mult_vect==K)
      temp_K <- c()
      temp_1_r_K <- c()
      if (D_m_star == q){   # no FP
        lower_bound <- 0
        term_1_r_K[[m]] <- 0
      }else{
        for (r in ((D_m_star+1):q)){    # Computation of the reccurence for each term of the sum
          temp <- (r-D_m_star)/r
          A <- c()
          B <- c()
          if (D_m_star == 0){
            for (j in (D_m_star+1):r){
              A <- c(A,2*(1-pnorm(q = sqrt(j*K), mean = 0, sd = 1)))
            }
            if (r>=2){
              for (j in (D_m_star+2):r){
                B <- c(B,2*(pnorm(q = sqrt(j*K), mean = 0, sd = 1)-pnorm(q = sqrt(K), mean = 0, sd = 1)))
              }
            }
          }else{
            #A
            for (j in (1:D_m_star)){
              A <- c(A,2-(pnorm(q = sqrt(j*K)-as.numeric(true_beta_matrix[j,k])/sqrt(sigma_2), mean = 0, sd = 1)+pnorm(q = sqrt(j*K)+as.numeric(true_beta_matrix[j,k])/sqrt(sigma_2), mean = 0, sd = 1)))
            }
            if (D_m_star>=2){
              #B
              for (j in 2:D_m_star){
                B <- c(B,pnorm(q = sqrt(j*K)-as.numeric(true_beta_matrix[j,k])/sqrt(sigma_2), mean = 0, sd = 1)+pnorm(q = sqrt(j*K)+as.numeric(true_beta_matrix[j,k])/sqrt(sigma_2), mean = 0, sd = 1)
                       -(pnorm(q = sqrt(K)-as.numeric(true_beta_matrix[j,k])/sqrt(sigma_2), mean = 0, sd = 1)+pnorm(q = sqrt(K)+as.numeric(true_beta_matrix[j,k])/sqrt(sigma_2), mean = 0, sd = 1)))
              }
            }
            for (j in (D_m_star+1):r){
              A <- c(A,2*(1-pnorm(q = sqrt(j*K), mean = 0, sd = 1)))
              B <- c(B,2*(pnorm(q = sqrt(j*K), mean = 0, sd = 1)-pnorm(q = sqrt(K), mean = 0, sd = 1)))
            }
          }
          temp_A_B <- A[1]
          if (length(B) > 0){    # If r !=1 (only possible if D_m_star =0)
            for (l in 1:length(B)){
              temp_A_B <- A[l+1]+B[l]*temp_A_B
            }
          }
          temp_1_r <- temp_A_B
          temp_1_r_K <- c(temp_1_r_K,temp_1_r)
          temp_K <- c(temp_K,temp*temp_1_r*list_P_2r[[D_m_star+1]][[m]][r-D_m_star])
        }
        lower_bound <- c(lower_bound,sum(temp_K))
        term_1_r_K[[m]] <- temp_1_r_K
      }
      ##### Evaluation of the theoretical lower bound with the available dataset and the used procedure only.
      for (it in 1:nbr_it_bound){
        # known variance
        lower_bound_hat_one_dataset_K2[[it]] <- c(lower_bound_hat_one_dataset_K2[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =2,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        lower_bound_hat_one_dataset_K3[[it]] <- c(lower_bound_hat_one_dataset_K3[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =3,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        lower_bound_hat_one_dataset_K2.5[[it]] <-  c(lower_bound_hat_one_dataset_K2.5[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =2.5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        lower_bound_hat_one_dataset_K3.5[[it]] <- c(lower_bound_hat_one_dataset_K3.5[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =3.5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        lower_bound_hat_one_dataset_K1.5[[it]] <- c(lower_bound_hat_one_dataset_K1.5[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =1.5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        lower_bound_hat_one_dataset_K1[[it]] <- c(lower_bound_hat_one_dataset_K1[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =1,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        lower_bound_hat_one_dataset_Klog_n[[it]] <- c(lower_bound_hat_one_dataset_Klog_n[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =log(n),sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        lower_bound_hat_one_dataset_K4[[it]] <- c(lower_bound_hat_one_dataset_K4[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =4,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        lower_bound_hat_one_dataset_K4.5[[it]] <- c(lower_bound_hat_one_dataset_K4.5[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =4.5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        lower_bound_hat_one_dataset_K5[[it]] <- c(lower_bound_hat_one_dataset_K5[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=TRUE))
        # unknown variance
        lower_bound_hat_one_dataset_K2_slope[[it]] <- c(lower_bound_hat_one_dataset_K2_slope[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =2,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        lower_bound_hat_one_dataset_K3_slope[[it]] <- c(lower_bound_hat_one_dataset_K3_slope[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =3,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        lower_bound_hat_one_dataset_K2.5_slope[[it]] <-  c(lower_bound_hat_one_dataset_K2.5_slope[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =2.5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        lower_bound_hat_one_dataset_K3.5_slope[[it]] <- c(lower_bound_hat_one_dataset_K3.5_slope[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =3.5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        lower_bound_hat_one_dataset_K1.5_slope[[it]] <- c(lower_bound_hat_one_dataset_K1.5_slope[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =1.5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        lower_bound_hat_one_dataset_K1_slope[[it]] <- c(lower_bound_hat_one_dataset_K1_slope[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =1,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        lower_bound_hat_one_dataset_Klog_n_slope[[it]] <- c(lower_bound_hat_one_dataset_Klog_n_slope[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =log(n),sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        lower_bound_hat_one_dataset_K4_slope[[it]] <- c(lower_bound_hat_one_dataset_K4_slope[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =4,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        lower_bound_hat_one_dataset_K4.5_slope[[it]] <- c(lower_bound_hat_one_dataset_K4.5_slope[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =4.5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
        lower_bound_hat_one_dataset_K5_slope[[it]] <- c(lower_bound_hat_one_dataset_K5_slope[[it]],lower_bound_estimated_one_dataset(cste_mult = K,path = path_collection[[k]],n,q,chosen_cste =5,sigma_2,m=it,cste_mult_vect,list_P_2r,known_variance=FALSE))
      }
    }
    #### For plotting the details of the theoretical lower bound. 
    if (detail == TRUE){
      ## We compare lower bound for i between 0 and (r-1)
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
             xlab = "K",ylab = "estimation of the [0,r-1]-FDR part")
        lines(cste_mult_vect, sapply(1:length(cste_mult_vect), function(K) term_1_r_K[[K]][r-D_m_star]), col = 2)
        legend("topright", legend = c(paste("|beta^*| = ",D_m_star),paste("r = ",r),"Proba of crit(m_r)<inf_{i in [0,r-1]}(crit(m_i))","theoretical_lower_bound_0_r-1"),
               lwd = c(0,0,1,1), col = c(3,3,4,2), lty = 1, cex=0.7)
      }
      for (r in (D_m_star+1):q){
        plot(cste_mult_vect, sapply(1:length(cste_mult_vect), function(K) crit_test_sup_r[[K]][r-D_m_star]), ylim = c(0,1), col = 4, type = 'p',
             xlab = "K",ylab = "estimation of the [r+1,q]-FDR part")
        lines(cste_mult_vect, sapply(1:length(cste_mult_vect), function(K) list_P_2r[[D_m_star+1]][[K]][r-D_m_star]), col = 2)
        legend("bottomright", legend = c(paste("|beta^*| = ",D_m_star),paste("r = ",r),"Proba of crit(m_r)<inf_{i in [r+1,q]}(crit(m_i))","theoretical_lower_bound_r+1_q"),
               lwd = c(0,0,1,1), col = c(3,3,4,2), lty = 1, cex=0.7)
      }
    }
    ##############
    setTxtProgressBar(pb,k)
    lower_bound_list[[k]] <- lower_bound
    lower_bound_list_hat_one_dataset_K2_list[[k]] <- lower_bound_hat_one_dataset_K2
    lower_bound_list_hat_one_dataset_K3_list[[k]] <- lower_bound_hat_one_dataset_K3
    lower_bound_list_hat_one_dataset_K2.5_list[[k]] <- lower_bound_hat_one_dataset_K2.5
    lower_bound_list_hat_one_dataset_K3.5_list[[k]] <- lower_bound_hat_one_dataset_K3.5
    lower_bound_list_hat_one_dataset_K1.5_list[[k]] <- lower_bound_hat_one_dataset_K1.5
    lower_bound_list_hat_one_dataset_K1_list[[k]] <- lower_bound_hat_one_dataset_K1
    lower_bound_list_hat_one_dataset_Klog_n_list[[k]] <- lower_bound_hat_one_dataset_Klog_n
    lower_bound_list_hat_one_dataset_K4_list[[k]] <- lower_bound_hat_one_dataset_K4
    lower_bound_list_hat_one_dataset_K4.5_list[[k]] <- lower_bound_hat_one_dataset_K4.5
    lower_bound_list_hat_one_dataset_K5_list[[k]] <- lower_bound_hat_one_dataset_K5
    lower_bound_list_hat_one_dataset_K2_list_slope[[k]] <- lower_bound_hat_one_dataset_K2_slope
    lower_bound_list_hat_one_dataset_K3_list_slope[[k]] <- lower_bound_hat_one_dataset_K3_slope
    lower_bound_list_hat_one_dataset_K2.5_list_slope[[k]] <- lower_bound_hat_one_dataset_K2.5_slope
    lower_bound_list_hat_one_dataset_K3.5_list_slope[[k]] <- lower_bound_hat_one_dataset_K3.5_slope
    lower_bound_list_hat_one_dataset_K1.5_list_slope[[k]] <- lower_bound_hat_one_dataset_K1.5_slope
    lower_bound_list_hat_one_dataset_K1_list_slope[[k]] <- lower_bound_hat_one_dataset_K1_slope
    lower_bound_list_hat_one_dataset_Klog_n_list_slope[[k]] <- lower_bound_hat_one_dataset_Klog_n_slope
    lower_bound_list_hat_one_dataset_K4_list_slope[[k]] <- lower_bound_hat_one_dataset_K4_slope
    lower_bound_list_hat_one_dataset_K4.5_list_slope[[k]] <- lower_bound_hat_one_dataset_K4.5_slope
    lower_bound_list_hat_one_dataset_K5_list_slope[[k]] <- lower_bound_hat_one_dataset_K5_slope
    print(k)
  }
  return(data.frame(lower_bound_list=I(lower_bound_list),
                    lower_bound_list_hat_one_dataset_K2_list=I(lower_bound_list_hat_one_dataset_K2_list),lower_bound_list_hat_one_dataset_K3_list=I(lower_bound_list_hat_one_dataset_K3_list),
                    lower_bound_list_hat_one_dataset_K2.5_list=I(lower_bound_list_hat_one_dataset_K2.5_list),lower_bound_list_hat_one_dataset_K3.5_list=I(lower_bound_list_hat_one_dataset_K3.5_list),
                    lower_bound_list_hat_one_dataset_K1.5_list=I(lower_bound_list_hat_one_dataset_K1.5_list),lower_bound_list_hat_one_dataset_K1_list=I(lower_bound_list_hat_one_dataset_K1_list),
                    lower_bound_list_hat_one_dataset_Klog_n_list=I(lower_bound_list_hat_one_dataset_Klog_n_list),lower_bound_list_hat_one_dataset_K4_list = I(lower_bound_list_hat_one_dataset_K4_list),
                    lower_bound_list_hat_one_dataset_K4.5_list = I(lower_bound_list_hat_one_dataset_K4.5_list),lower_bound_list_hat_one_dataset_K5_list = I(lower_bound_list_hat_one_dataset_K5_list),
                    lower_bound_list_hat_one_dataset_K2_list_slope=I(lower_bound_list_hat_one_dataset_K2_list_slope),lower_bound_list_hat_one_dataset_K3_list_slope=I(lower_bound_list_hat_one_dataset_K3_list_slope),
                    lower_bound_list_hat_one_dataset_K2.5_list_slope=I(lower_bound_list_hat_one_dataset_K2.5_list_slope),lower_bound_list_hat_one_dataset_K3.5_list_slope=I(lower_bound_list_hat_one_dataset_K3.5_list_slope),
                    lower_bound_list_hat_one_dataset_K1.5_list_slope=I(lower_bound_list_hat_one_dataset_K1.5_list_slope),lower_bound_list_hat_one_dataset_K1_list_slope=I(lower_bound_list_hat_one_dataset_K1_list_slope),
                    lower_bound_list_hat_one_dataset_Klog_n_list_slope=I(lower_bound_list_hat_one_dataset_Klog_n_list_slope),lower_bound_list_hat_one_dataset_K4_list_slope = I(lower_bound_list_hat_one_dataset_K4_list_slope),
                    lower_bound_list_hat_one_dataset_K4.5_list_slope = I(lower_bound_list_hat_one_dataset_K4.5_list_slope),lower_bound_list_hat_one_dataset_K5_list_slope = I(lower_bound_list_hat_one_dataset_K5_list_slope)))
}