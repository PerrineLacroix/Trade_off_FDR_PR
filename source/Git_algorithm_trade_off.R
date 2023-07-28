# This function source encodes the algorithm to select the multiplicative constant \Hat{K}
# achieving a trade-off between PR and FDR.
# Thresholds on PR and FDR are given by the user.
trade_off_PR_FDR <- function(cste_for_estimation,alpha_FDR,gamma_PR,list_upper_bounds,nbr_it_bound,choice_random_collection,
                             PR_FDP_hat_K,cste_mult_vect,sigma_2){
  I_1 <- list()     # interval of K such that FDR_estimated < alpha_FDR
  I_2 <- list()       # interval of K such that PR_estimated < gamma_FDR
  intersection <- list()        # intersection of interval of K such that FDR_estimated < alpha_FDR and interval of K such that PR_estimated < gamma_FDR
  min_intersection <- list()     # min of the previous intersection
  K_hat <- list()       # K_hat is either the min of the previous intersection, or min of the interval of K such that FDR_estimated < alpha_FDR (if the first intersection is empty)
  I_1_theo <- list()    # interval of K such that FDR_empirical < alpha_FDR
  I_2_theo <- list()        # interval of K such that PR_empirical < alpha_FDR
  intersection_theo <- list()    # intersection of interval of K such that FDR_empirical < alpha_FDR and interval of K such that PR_empirical < alpha_FDR
  # estimated FDR and PR with \Tilde{K} = 4
  for (k in choice_random_collection){
    I_1[[k]] <- list()
    I_2[[k]] <- list()
    intersection[[k]] <- list()
    min_intersection[[k]] <- list()
    K_hat[[k]] <- list()
    for (it in 1:nbr_it_bound){
      vect_estimation_FDR <- as.numeric(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[it]][seq(1,length(list_upper_bounds$upper_bound_list_hat_one_dataset_K2_list_slope[[k]][[it]]), by=2)])
      vect_estimation_diff_PR <- PR_FDP_hat_K[[k]]$unknown_variance$diff_PR_one_dataset_slope
      I_1[[k]][[it]] <- cste_mult_vect[which(vect_estimation_FDR<alpha_FDR)]
      I_2[[k]][[it]] <- cste_mult_vect[which(abs(vect_estimation_diff_PR)<gamma_PR*sqrt(list_upper_bounds$upper_bound_list_hat_one_dataset_K4_list_slope[[k]][[it]]$kappa_slope))]
      intersection[[k]][[it]] <- intersect(I_1[[k]][[it]],I_2[[k]][[it]])
      if (length(intersection[[k]][[it]])>0){
        min_intersection[[k]][[it]] <- min(intersect(I_1[[k]][[it]],I_2[[k]][[it]]))
        K_hat[[k]][[it]] <- min(intersect(I_1[[k]][[it]],I_2[[k]][[it]]))
      }else{
        min_intersection[[k]][[it]] <- NA
        K_hat[[k]][[it]] <- min(I_1[[k]][[it]])
      }
    }
  }
  # empirical FDR and PR
  for (k in choice_random_collection){
    vect_empirical_FDR <- sapply(1:length(cste_mult_vect), function(l) mean(PR_FDP_hat_K[[k]]$empirical$FDP_hat[,l]))
    vect_empirical_PR <- sapply(1:length(cste_mult_vect), function(l) mean(PR_FDP_hat_K[[k]]$empirical$PR_hat[,l]))
    vect_empirical_diff_PR <- vect_empirical_PR - vect_empirical_PR[which(cste_mult_vect==2)]
    I_1_theo[[k]] <- cste_mult_vect[which(vect_empirical_FDR<alpha_FDR)]
    I_2_theo[[k]] <- cste_mult_vect[which(abs(vect_empirical_diff_PR)<gamma_PR*sigma_2)]
    intersection_theo[[k]] <- intersect(I_1_theo[[k]],I_2_theo[[k]])
  }
 return(list(FDR_estime = I_1,PR_estime = I_2, FDR_theo = I_1_theo, PR_theo = I_2_theo, min_intersection = min_intersection, K_hat = K_hat))
}
