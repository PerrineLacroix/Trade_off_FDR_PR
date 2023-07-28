## script to launch simulations
## author : Perrine Lacroix
## date : March, 10 2021

## This code creates the data sets and the deterministic ordered model collections. 
## Two groups of data sets are generated and each countains nbr_it independant data sets. Data sets for nbr_it = 1 is used for the estimated quantities, all of data sets are used for the empirical quantity estimations. 
## The first group of data sets is used to get the estimated FDR and the estimated PR functions, as well as the empirical estimations of FDR and PR
## The second group of data sets is used to test the FDR bounds construction with their variability. 
## For each data set, the fixed design matrix is firstly built and a Gaussian noise is generated. 
## Then, Y is obtained after the beta_0 construction providing an order of relevance on the X_j. 
## Lastly, the nested model collection respecting the order of relevance on the X_j is generated and characteristics of each model are computed. 


rm(list=objects())

prefix = "~/Documents/GitHub/Trade_off_FDR_PR"
setwd(paste(prefix,"/data",sep =""))
# Parameters of simulation
nbr_it = 1000   # Data sets number for the empirical estimations of FDR and PR. 
nbr_it_bound = 100   # Data sets number to compute the FDR bounds and to study the variability of the FDR bounds.
n = 50  # number of observations
p = 50   # number of variables
n_temps = 2*n*nbr_it   # Required total iterations number to get nbr_it train sets for the empirical estimation of FDR, with nbr_it validation sets for the empirical estimation of PR. 
n_temps_bounds = n*nbr_it_bound   # Required total iterations number to get nbr_it_bound data sets for the FDR bounds and the variability studies of the FDR bounds.
sigma_2 <- 1   # Variance of the noise in the Gaussian linear regression

## parameters of the created data for the saving
config_data = paste0(paste("sigma2",sigma_2,"p",p,"n",n,"nbr_it",nbr_it,sep="_"),".RData")



#####################################################################################################################
###### Creation of sets used for estimated and empirical estimations. 
#####################################################################################################################

###### Creation of the matrix X and noise: X is fixed for all iterations since the design is supposed to be fixed in the model.
## Here, the design matrix is fixed and we choose the usual orthonormal matrix for each iteration. 
## Only the noise is random.
set.seed(1234)
## X matrix
data_select <- matrix(0,nrow=n,ncol=p)   # the train matrix
diag(data_select) <- rep(1,length(diag(data_select)))
data_metric <- matrix(0,nrow=n,ncol=p)  # the validation matrix
diag(data_metric) <- rep(1,length(diag(data_metric)))
colnames(data_select) <-  paste(rep("covariable",p),1:p)
rownames(data_select) <-  paste(rep("observation",n),1:n)
colnames(data_metric) <-  colnames(data_select)
rownames(data_metric) <-  rownames(data_select)
## noise
noise <- rnorm(2*nbr_it*n,0,sqrt(sigma_2))
## Saving the data sets
save(data_select,data_metric,noise,file=paste("matrix_X",config_data,sep="_"))

###### Creation of the beta_0 vectors and the response variables Y for the empirical estimations 
## beta_0 construction
## We study all support sizes D between 0 and 20
## For each size D, the beta coefficients are genereated as follow:
## - The coordinate number D is choosen to be higher than the noise (or not)
## - The coordinate j<D is choosen (much or not) higher than the j+1 one
setwd(paste(prefix,"/data",sep =""))
load(paste("matrix_X",config_data,sep="_"))
config_data = paste0(paste("beta_true_0_20","sigma2",sigma_2,
                           "p",p,"n",n,"nbr_it",nbr_it,sep="_"),".RData")
vrai_beta_matrix <- matrix(0,nrow = p,ncol =21)
colnames(vrai_beta_matrix) <- sapply(0:20, function(j) paste("size_beta_0",j,sep = "_"))
vrai_beta_matrix[1,2] <- 2*sigma_2
# vrai_beta_matrix[1,2] <- sigma_2/10    # the configuration 2 of the scenario (ii)  of the article
# vrai_beta_matrix[1,2] <- 2*sigma_2     # the configuration 3 of the scenario (ii)  of the article
for (k in 2:20){  # size of beta_0 support
  index <- k+1
  vrai_beta_matrix[(index-1),index] <- 2
  # vrai_beta_matrix[(index-1),index] <- sigma_2/10      # the configuration 2 of the scenario (ii)  of the article
  for (j in (index-2):1){
    vrai_beta_matrix[j,index] <- round(runif(1,as.numeric(vrai_beta_matrix[j+1,index]+0.5),as.numeric(vrai_beta_matrix[j+1,index])+1.5),3)
    # vrai_beta_matrix[j,index] <- round(runif(1,as.numeric(vrai_beta_matrix[j+1,index]+0.05),as.numeric(vrai_beta_matrix[j+1,index])+0.15),3)   # the configurations 2 and 3 of the scenario (ii)  of the article
  }
}
## Y construction
Y_select_matrix <- matrix(NA,nrow = n*nbr_it, ncol = ncol(vrai_beta_matrix))  # the train sets of Y
Y_metric_matrix <- matrix(NA,nrow = n*nbr_it, ncol = ncol(vrai_beta_matrix))    # the validation sets of Y
colnames(Y_select_matrix) <- sapply(0:20, function(j) paste("size_beta_0",j,sep = "_"))
colnames(Y_metric_matrix) <- sapply(0:20, function(j) paste("size_beta_0",j,sep = "_"))
for (k in seq(1,n*nbr_it,by=n)){
  Y_select_matrix[k:(k+(n-1)),] <- sapply(1:ncol(vrai_beta_matrix), function(j) data_select %*% vrai_beta_matrix[,j] + noise[k:(k+(n-1))])
  Y_metric_matrix[k:(k+(n-1)),] <- sapply(1:ncol(vrai_beta_matrix), function(j) data_metric %*% vrai_beta_matrix[,j] + noise[(n*nbr_it+k):(n*nbr_it+k+(n-1))])
}
## Saving the data sets
save(Y_select_matrix,Y_metric_matrix,vrai_beta_matrix,file=paste("vector_Y",config_data,sep="_"))
# save(Y_select_matrix,Y_metric_matrix,vrai_beta_matrix,file=paste("smaller_noise_vector_Y",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# save(Y_select_matrix,Y_metric_matrix,vrai_beta_matrix,file=paste("larger_noise_closed_values_vector_Y",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article


#####################################################################################################################
###### Creation of sets used for FDR bounds and the study of their variability.
#####################################################################################################################
# Creation of the response variables Y for the FDR bounds and the variability studies of the FDR bounds.
# X and beta_0 are fixed. Only the noise varies.
setwd(paste(prefix,"/data",sep =""))
config_data = paste0(paste("sigma2",sigma_2,
                           "p",p,"n",n,"nbr_it",nbr_it,sep="_"),".RData")
load(paste("matrix_X",config_data,sep="_"))
config_data = paste0(paste("beta_true_0_20","sigma2",sigma_2,"p",p,"n",n,"nbr_it",nbr_it,sep="_"),".RData")
load(paste("vector_Y",config_data,sep="_"))
# load(paste("smaller_noise_vector_Y",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("larger_noise_closed_values_vector_Y",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article

## noise
noise_bounds <- rnorm(nbr_it_bound*n,0,sqrt(sigma_2))
## Y_bounds
Y_bounds_matrix <- matrix(NA,nrow = n*nbr_it_bound, ncol = ncol(vrai_beta_matrix))
colnames(Y_bounds_matrix) <- sapply(0:20, function(j) paste("size_beta_0",j,sep = "_"))
for (k in seq(1,n*nbr_it_bound,by=n)){
  Y_bounds_matrix[k:(k+(n-1)),] <- sapply(1:ncol(vrai_beta_matrix), function(j) data_select %*% vrai_beta_matrix[,j] + noise_bounds[k:(k+(n-1))])
}
save(Y_bounds_matrix,noise_bounds,file=paste("vector_Y_noise_bounds",config_data,sep="_"))
# save(Y_bounds_matrix,noise_bounds,file=paste("vector_Y_smaller_noise_noise_bounds",config_data,sep="_"))   # the configuration 2 of the scenario (ii)  of the article
# save(Y_bounds_matrix,noise_bounds,file=paste("vector_Y_larger_noise_closed_values_bounds",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article


#####################################################################################################################
###### Creation of the model collections
#####################################################################################################################
setwd(paste(prefix,"/data",sep =""))
config_data = paste0(paste("sigma2",sigma_2,
                           "p",p,"n",n,"nbr_it",nbr_it,sep="_"),".RData")
load(paste("matrix_X",config_data,sep="_"))
config_data = paste0(paste("beta_true_0_20","sigma2",sigma_2,"p",p,"n",n,"nbr_it",nbr_it,sep="_"),".RData")
load(paste("vector_Y",config_data,sep="_"))
# load(paste("smaller_noise_vector_Y",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("larger_noise_closed_values_vector_Y",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article
load(paste("vector_Y_noise_bounds",config_data,sep="_"))
# load(paste("vector_Y_smaller_noise_bounds",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# load(paste("vector_Y_larger_noise_closed_values_bounds",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article


###### Creation of the model collections used for empirical estimations. 
## According to the construction of beta_0, X_1 is the most relevant variable, then X_2, then X_3, ...
## Models are nested in the collection. The first model m_1 is the empty model, m_2 = Vect(X_1), m_3 = Vect(X_1,X_2), ...
## Collection for the empirical estimations sets
dimension <- 0:p   
path_collection <- list()
support <- list()   # list for the supports of variable in the collection
support_num <- list()    # list for the supports of variable index in the collection
for(k in 2:(p+1)){
  support[[k]] <- colnames(data_metric)[1:(k-1)]
  support_num[[k]] <- 1:(k-1)
}
pb = txtProgressBar(min = 1, max = ncol(vrai_beta_matrix), initial = 1)
for (l in 1:(ncol(Y_select_matrix))){
  Y_select <- Y_select_matrix[,l]
  path_collection[[l]] <- list()
  for (k in seq(1,n*nbr_it,by=n)){
    j = which(seq(1,n*nbr_it,by=n) ==k)
    set.seed(j+1234)
    index <- 2:(p+1)
    re_estimation <- list()   
    for (i in index){
      re_estimation[[i]] <- lm(Y_select[k:(k+(n-1))] ~ data_select[,support_num[[i]]] -1)  # No intercept, mean square estimators onto each model
            # The beta coefficients are re-estimated onto each model of the collection by the MSE minization
    }
    beta_new <- sapply(re_estimation, function(a) as.vector(a$coefficients))  
    for (i in index){
      names(beta_new[[i]]) <- support[[i]]
    }
    beta_estimator <- matrix(0, nrow = length(support), ncol = p)
    for (i in 2:length(support)){
      beta_estimator[i,support_num[[i]]] <- as.numeric(beta_new[[i]])
    }
    LS = sapply(1:length(dimension), function(i)
      (1/nrow(data_select))*sum((Y_select[k:(k+(n-1))] - data_select %*% beta_estimator[i,])^2))   # The least squared values are computed for each beta of the collection
    data_frame_temp_1 <- data.frame(LS=LS,dim=dimension,complexite=lchoose(p,dimension))      # characteristics of the model collection
    data_frame_temp_2 <- data.frame(beta=beta_estimator)
    path_collection[[l]][[j]]=list(val1 = data_frame_temp_1, val2 = data_frame_temp_2)
  }
  setTxtProgressBar(pb,l)
}
save(path_collection,file=paste("path_collection",config_data,sep="_"))
# save(path_collection,file=paste("path_collection_smaller_noise",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# save(path_collection,file=paste("path_collection_larger_noise_closed_values",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article


###### Creation of the model collections used for FDR bound evaluations. 
path_collection_bounds <- list()
support_bounds <- list()    # list for the supports of variable in the collection
support_num_bounds <- list()    # list for the supports of variable index in the collection
for(k in 2:(p+1)){
  support_bounds[[k]] <- colnames(data_metric)[1:(k-1)]
  support_num_bounds[[k]] <- 1:(k-1)
}
pb = txtProgressBar(min = 1, max = ncol(vrai_beta_matrix), initial = 1) 
for (l in 1:(ncol(Y_bounds_matrix))){
  Y_bounds <- Y_bounds_matrix[,l]
  path_collection_bounds[[l]] <- list()
  for (k in seq(1,n*nbr_it_bound,by=n)){
    j = which(seq(1,n*nbr_it_bound,by=n) ==k)
    set.seed(j+1234)
    index <- 2:(p+1)
    re_estimation <- list()
    for (i in index){
      re_estimation[[i]] <- lm(Y_bounds[k:(k+(n-1))] ~ data_select[,support_num_bounds[[i]]] -1)  # No intercept
                  # The beta coefficients are re-estimated onto each model of the collection by the MSE minization
    }
    beta_new <- sapply(re_estimation, function(a) as.vector(a$coefficients))
    for (i in index){
      names(beta_new[[i]]) <- support_bounds[[i]]
    }
    beta_estimator <- matrix(0, nrow = length(support_bounds), ncol = p)
    for (i in 2:length(support_bounds)){
      beta_estimator[i,support_num_bounds[[i]]] <- as.numeric(beta_new[[i]])
    }
    LS = sapply(1:length(dimension), function(i)
      (1/nrow(data_select))*sum((Y_bounds[k:(k+(n-1))] - data_select %*% beta_estimator[i,])^2))   # The least squared values are computed for each beta of the collection
    data_frame_temp_1 <- data.frame(LS=LS,dim=dimension,complexite=lchoose(p,dimension))    # characteristics of the model collection
    data_frame_temp_2 <- data.frame(beta=beta_estimator)
    path_collection_bounds[[l]][[j]]=list(val1 = data_frame_temp_1, val2 = data_frame_temp_2)
  }
  setTxtProgressBar(pb,l)
}
save(path_collection_bounds,file=paste("path_collection_bounds",config_data,sep="_"))
# save(path_collection_bounds,file=paste("path_collection_bounds_smaller_noise",config_data,sep="_"))    # the configuration 2 of the scenario (ii)  of the article
# save(path_collection_bounds,file=paste("path_collection_bounds_larger_noise_closed_values",config_data,sep="_"))    # the configuration 3 of the scenario (ii)  of the article