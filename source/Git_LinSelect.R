# This function source applies the model selection procedure with the LinSelect penalty function. 
# It applies point by point the R function tuneLasso of the R package LINselect.
LinSelect_function <- function(res,n,p,dmax_LIN,max.steps_LIN,coeff_mult_proj_LIN,coeff_mult_pen_LIN,Y_s,data_s){
  calc.proj <- FALSE
  calc.pen <- FALSE
  result <- NULL
  # maximum number of steps il the lasso procedure
  if (is.null(max.steps_LIN)){
    max.steps_LIN = 2 * min(p, n)   # by default
  }
  # maximum number of variables in the lasso estimator
  if (p >= n){
    dmax_LIN <- floor(min(c(3 * p/4, n - 5, dmax_LIN)))    # default for high dimension
  }
  if (p < n){
    dmax_LIN <- floor(min(c(p, n - 5, dmax_LIN)))      # default for non high dimension
  }
  # complexity
  if (max(res$val1$dim) <= dmax_LIN){
    Nmod.lasso <- length(res$val1$dim)
  }
  if (max(res$val1$dim) > dmax_LIN){
    Nmod.lasso <- (1:length(res$val1$dim))[(res$val1$dim >= dmax_LIN)][1]
  }
  # Obtaining the X beta_chap on the Lambda grid
  f.lasso <- matrix(0, nrow = n, ncol = Nmod.lasso)
  beta_estimateur <- as.matrix(res$val2)
  mu = mean(Y_s)
  X_mean <- apply(data_s,2,mean)
  Intercept <- mu - beta_estimateur[1:Nmod.lasso,] %*% X_mean
  for (l in 1:Nmod.lasso){
    f.lasso[, l] <- rep(Intercept[l],n) + data_s %*% beta_estimateur[l,]
  }
  ## 'Lasso' method
  result <- list(lasso = NULL)
  calc.proj <- TRUE
  D <- 0:dmax_LIN
  dm <- pmin(D, rep(p/2, dmax_LIN + 1))
  Delta <- lgamma(p+1) - lgamma(dm + 1) - lgamma(p - dm + 1) + 2 * log(D + 1)
  pen <- penalty(Delta, n, p, coeff_mult_pen_LIN)
  calc.pen <- TRUE
  I.lasso <- list(NULL)
  ProjMod <- array(0, c(n, n, Nmod.lasso))
  A <- array(Inf, c(Nmod.lasso, Nmod.lasso))
  B = A
  SCR <- rep(0, Nmod.lasso)
  penSCR <- rep(0, Nmod.lasso)
  sumY2 <- sum((Y_s - mu)^2)
  un <- rep(1, n)
  ## Linselect criterion calculation to minimize
  for (l in 1:Nmod.lasso) {
    I.lasso[[l]] <- (1:p)[beta_estimateur[l, ] != 0]
    if (length(I.lasso[[l]]) > 0) {
      ProjM <- LINselect:::Proj(cbind(un, data_s[, I.lasso[[l]]]), length(I.lasso[[l]]) + 1)
      ProjMod[, , l] <- ProjM$Proj
      SCR[l] <- sum((Y_s - ProjMod[, , l] %*% Y_s)^2)
      penSCR[l] <- SCR[l] * pen[ProjM$rg + 1]/(n - ProjM$rg)
    }
    if (length(I.lasso[[l]]) == 0) {
      ProjMod[, , l] <- un %*% t(un)/n
      ProjM <- list(Proj = ProjMod[, , l], rg = 1)
      SCR[l] <- sumY2
      penSCR[l] <- sumY2 * pen[ProjM$rg + 1]/(n - ProjM$rg)
    }
  }
  Ind <- rep(1:Nmod.lasso, rep(Nmod.lasso, Nmod.lasso))
  YY <- Y_s %*% t(rep(1, Nmod.lasso))
  for (m in 1:Nmod.lasso) {
    Proj.f <- ProjMod[, , m] %*% f.lasso
    B[m, ] <- apply((YY - Proj.f)^2, 2, sum) + coeff_mult_proj_LIN * apply((f.lasso - Proj.f)^2, 2, sum)
    A[m, ] <- B[m, ] + penSCR[m]
  }
  l.lasso <- Ind[which.min(A)]
  return(l.lasso)
}