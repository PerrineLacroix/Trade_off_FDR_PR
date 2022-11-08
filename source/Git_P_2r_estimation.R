# This function source is devoted to the P_2r terms estimation for each K of the considered grid.
P_2r_estimation <- function(n,p,asympt_number,cste_mult_vect){
  q = min(n,p) 
  P_2r_list <- list()
  set.seed(1234)
  list_standard_gaussian <- lapply(1:asympt_number, function(l) rnorm(q,0,1))   # a Gaussian vector
  for (D_m_star in 0:q){    # the P_2r term is computed for all r in {0,...,q}
    term_r_plus_1_q_K <- list()
    for (K in cste_mult_vect){
      m = which(cste_mult_vect==K)
      temp_r_plus_1_q_K <- c()
      if (D_m_star == q){  # no FP
        term_r_plus_1_q_K[[m]] <- 1  # Probability of an empty set is 1 per default
      }else{
        if (D_m_star < (q-1)){
          for (r in ((D_m_star+1):(q-1))){    # si r = q, l'evenement est vide
            sum_gaussian <- lapply(1:asympt_number, function(it) sapply((r+1):q, function(l) sum((list_standard_gaussian[[it]][(r+1):l])^2)))   # the chi_squared vector from the Gaussian one
            vect_chi_2 <- sapply((r+1):q, function(l) K*(l-r))
            temp_r_plus_1_q <- mean(sapply(1:asympt_number, function(it) all(sapply(1:length((r+1):q), function(l) sum_gaussian[[it]][l] < vect_chi_2[l]))))    # estimation of the P_2r term
            temp_r_plus_1_q_K <- c(temp_r_plus_1_q_K,temp_r_plus_1_q)
          }
        }
        temp_r_plus_1_q_K <- c(temp_r_plus_1_q_K,1)   # ensemble vide
        term_r_plus_1_q_K[[m]] <- temp_r_plus_1_q_K
      }
    }
    P_2r_list[[D_m_star+1]] <- term_r_plus_1_q_K
    print(D_m_star)
  }
  return(P_2r_list)
}