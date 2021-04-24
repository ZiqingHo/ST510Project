library(dbscan)
library(dplyr)
library(MASS)
library(mixAK)
library(mclust)
library(ggplot2)
library(ggpubr)
library(parallel)
library(RSKC)
library(sparcl)

# Generate synthetic dataset with outliers and noise variables
# we should know which variables are noise variables and 
# which observations are considered as outliers in signal and noise variables, respectively
# Follow the notation used in the WRSK paper
# K - number of clusters
# cluster_size - size of each cluster
# p_inf - number of informative variables
# p_noise - number of noise variables
# p_inf_out - number of informative variables contaminated by either scattered outliers or uniformly distributed outliers
# out_type - uniformly distributed outliers (TRUE) or scattered outliers (FALSE)
# inf_out - the percentage of contaminated observations in informative variables 
# p_noise_out - number of noise variables contaminated by uniformly distributed outliers
# noise_out - the percentage of contaminated observations in noise variables 
# unif_interval - interval of uniform distribution to generate uniformly distributed outliers
# mu_interval - interval of uniform distribution to generate mu of each cluster
# scatter_interval - interval of uniform distribution to generate sigma for the generation of scattered outliers
# rho_interval - interval of unifrom distribution to generate rho
# mu_cluster - mean vector of the clusters
# cov_cluster - covariance matrix of each cluster
# n_inf_out - number of observations that have contaminated informative variables in each cluster
# s_sigma - covariance structure used to generate scattered outliers

# Generate a dataset with contaminated informative variables
# which includes the true cluster labels before contamination
# and the outlier labels
# so we will know the outliers belong to which cluster

SimData <- function(K, cluster_size = NULL,
                    p_inf, p_noise = 0,
                    p_inf_out = 0, out_type = TRUE, inf_out = 0.10,
                    p_noise_out = 0, noise_out = 0.10,
                    unif_interval = c(-12,-6,6,12),
                    mu_interval = c(-6,-3,3,6),
                    scatter_interval = c(3,10),
                    rho_interval = c(0.1,0.9)){
  
  # dimension of dataset
  p <- ifelse(p_noise==0, p_inf, p_inf + p_noise)
  
  # size of each cluster is set to be 50 
  if(is.null(cluster_size)){
    cluster_size <- rep(50, K)
  }
  
  if(K != length(cluster_size)){
    stop("K must be equal to the length(cluster_size)!")
  }
  
  # total number of observations
  n <- sum(cluster_size)
  
  # Follow Eq (18): Q%*%rho_matrix%*%t(Q)
  # Generate covariance matrix for each cluster
  # Generate off-diagonal elements of rho_matrix, 
  # follow uniform distribution
  off_diag_Rho <- runif(K, min = rho_interval[1], max = rho_interval[2])
  
  # Generate K different covariance matrices
  # use list to store them
  cov_cluster <- list()
  
  for(i in 1:K){
    # Off-diagonal elements of rho_matrix = rho_i
    Rho <- matrix(off_diag_Rho[i], ncol = p_inf, nrow = p_inf)
    # All the diagonal elements of rho_matrix = 1
    diag(Rho) <- 1
    # Q is a random rotation matrix satisfying t(Q) = inv(Q)
    Q <- rRotationMatrix(n = 1, dim = p_inf)
    # Eq (18)
    cov_cluster[[i]] <- Q %*% Rho %*% t(Q)
  }
  
  # number of contaminated observations in informative variables
  n_inf_out <- c(round(cluster_size*inf_out))
  
  # Generate scattered outliers
  # sigma of scattered outliers follows Uniform distribution
  s_sigma <- runif(K, min = scatter_interval[1], 
                   max = scatter_interval[2])
  
  # Generate the dataset that contaminated in the informative variables
  # arbitrary vectors for true cluster label before contamination and outlier labels
  # arbitrary vector for synthetic dataset (or use arbitrary matrix)
  y <- label <- x <- c()
  
  # mean vector for each cluster, see Eq. 17
  # MU is a matrix with the mean vector for each cluster as the column vector
  # Eg: K = 3, p_inf = 5, mu_1 = (mu,0,0,mu,0), 
  #     mu_2 = (0,mu,0,0,mu), mu_3 = (0,0,mu,0,0)
  
  # Generate the mean vector for each cluster
  # follow uniform distribution
  mu_cluster <- rep(sample(c(runif(K, min = mu_interval[1], max = mu_interval[2]),
                             runif(K, min = mu_interval[3], max = mu_interval[4])),1),
                    K)
  
  if(p_inf<K){
    MU <- matrix(0, nrow = K, ncol = K)
    diag(MU) <- mu_cluster
    MU <- MU[1:p_inf,]
  }else{
    # number of repetitions
    n_repeat <- p_inf%/%K
    # get the first K by K matrix then duplicate MU matrix n_repeat times
    MU <- matrix(0, nrow = K, ncol = K)
    diag(MU) <- mu_cluster
    MU <- matrix(rep(MU, n_repeat), ncol = K, byrow = T)
    
    # if remainder of p_inf/K != 0, we need to add mu to get a full mean matrix
    if(nrow(MU) != p_inf){
      n_remainder <- p_inf%%K
      MU_remainder <- matrix(0,ncol = K, nrow = n_remainder)
      diag(MU_remainder) <- mu_cluster[1:n_remainder]
      MU <- rbind(MU,MU_remainder)
    }
  }
  
  # Generate the dataset by cluster
  for(i in 1:K){
    # true cluster labels before contamination
    y <- c(y,rep(i,cluster_size[i]))
    
    # covariance matrix of i-th cluster
    cov_structure <- cov_cluster[[i]]
    
    # labels for one specific cluster, 
    # later some cluster labels will be replaced by outlier labels (0)
    label_cluster <- rep(i,cluster_size[i])
    
    # Each cluster is explained by informative variables and 
    # is generated by Gaussian model with mean vector of i-th cluster 
    # and its covariance structure
    X <- mvrnorm(n = cluster_size[i], mu = MU[,i], Sigma = cov_structure)
    
    # Generate noise variables by using univariate N(0,1)
    if(p_noise != 0){
      X <- cbind(X, matrix(rnorm(p_noise*cluster_size[i]), ncol = p_noise))
      # Now X is the dataset before contamination
    }
    
    # Generate outliers in informative variables:
    # contaminate the first n_inf_out[i] observations in the informative variables
    if(n_inf_out[i] != 0){
      if(p_inf_out != 0){
        contaminated_inf <- 1:p_inf_out
      }else{
        stop("The number of informative variables to be contaminated is missing!")
      }
      
      # number of contaminated informative variables
      nc = length(contaminated_inf)
      
      # covariance structure of scatter outliers = sigma*identity matrix
      S_out <- matrix(0, ncol = nc, nrow = nc)
      diag(S_out) <- s_sigma[i]
      
      # Generate uniformly distributed outliers
      # follow uniform distribution
      if(out_type == T){
        X[1:n_inf_out[i],contaminated_inf] <- matrix(sample(c(runif(n_inf_out[i]*nc,
                                                                    min = unif_interval[1], 
                                                                    max = unif_interval[2]),
                                                              runif(n_inf_out[i]*nc,
                                                                    min = unif_interval[3], 
                                                                    max = unif_interval[4])),
                                                            n_inf_out[i]*nc),
                                                     nrow = n_inf_out[i], 
                                                     ncol = nc)
        # the labels of the first n_inf_out[i] observations are replaced by outlier labels (0)
        label_cluster[1:n_inf_out[i]] <- 0
        
      }else{ 
        # Generate scattered outliers
        # follow Gaussian model with the same mean vector 
        # but different covariance structure
        X[1:n_inf_out[i],contaminated_inf] <- mvrnorm(n_inf_out[i], mu = MU[contaminated_inf,i], Sigma = S_out)
        # the labels of the first n_inf_out[i] observations are replaced by outlier labels (0)
        label_cluster[1:n_inf_out[i]] <- 0
      }
      
    }
    # combine the dataset and labels of each cluster repeatedly
    x <- rbind(x,X)
    label <- c(label,label_cluster)
  }
  
  # Generate outliers in noise variables
  # this process doesn't make sense if there is no noise variable
  if(p_noise != 0){
    if(p_noise_out != 0){
      # randomly choose the noise variables to be contaminated
      contaminated_noise <- sample((p_inf + 1):p, p_noise_out)
      
      # number of contaminated observations in noise variables
      n_noise_out <- c(round(cluster_size*noise_out))
      
      # the observations are randomly selected
      # they are differ from those contaminated in the informative variables      
      # select equal amount of observations contaminated in 
      # noise variables from each cluster
      first_index <- c(1 + c(n_inf_out, 0) + c(0, cumsum(cluster_size)))[1:K]
      
      # construct a matrix consists of the first (1st column) and 
      # last indices (2nd column) of observation from each cluster 
      # that can be contaminated in noise variables
      # last column is the number of contaminated observations
      c_noise_index <- cbind(first_index,cumsum(cluster_size),n_noise_out)
      
      # randomly select the observations to be contaminated in noise variables
      x_noise_index <- apply(c_noise_index,1, function(x){
        sample(c(x[1]:x[2]),x[3])})
      
      # get a vector of indices of selected observations
      x_noise_index <- unlist(c(x_noise_index))
      nc_noise <- length(x_noise_index)
      
      # Generate outliers in noise variables:
      # contaminate the observations in the noise variables
      # follow uniform distribution, no scattered outliers
      x[x_noise_index, contaminated_noise] <- matrix(sample(c(runif(nc_noise*p_noise_out,
                                                                    min = unif_interval[1],
                                                                    max = unif_interval[2]),
                                                              runif(p_noise_out*nc_noise,
                                                                    min = unif_interval[3],
                                                                    max = unif_interval[4])),
                                                            p_noise_out*nc_noise),
                                                     nrow = nc_noise,
                                                     ncol = p_noise_out)
      
      # assign outlier labels
      label[x_noise_index] <- 0
      
      # a vector storing the cluster labels that uncontaminated and outlier labels
      # in informative variables (0) 
      # and noise variables (K + 1)
      label_out <- label
      label_out[x_noise_index] <- K+1
    }
    else{
      # no contamination in noise variables
      label_out <- label
    }
  }
  
  # assign column names
  colnames(x) <- paste0("x", rep(1:p))
  # assign row names
  rownames(x) <- rep(1:n)
  
  # X - synthetic dataset, y - true cluster labels
  # label - cluster labels and outlier labels
  #         without distinguishing the contamination in informative 
  #         and noise variables
  # outliers - cluster labels and outlier labels in informative and noise variables
  return(list(X = x, y = y,
              label = label, outliers = label_out))
  
}

# L1 norm - to select sparsity parameter, s
L1norm <- function(x){
  return(sum(abs(x)))
}

# L2 norm
L2norm <- function(x){
  return(sqrt(sum(x^2)))
}

# Soft-thresholding
soft_thres <- function(x, y){
  return(pmax(0, abs(x) - y) * sign(x))
}

# use bisection method to find delta
# diff_var = total sum of squares_variable - within cluster sum of squares_variable
# s - sparsity parameter
# n_iter - maximum number of iterations
# crit - stopping criterion (follow Tibshirani's sparse k-means)
bisection <- function(diff_var, s, n_iter = 10, crit = 1e-4){
  # L2 norm
  L2 <- L2norm(diff_var)
  
  # Lasso penalty
  # this condition follows Step 2b in Tibshirani's sparse k-means
  if(L2 == 0 || L1norm(diff_var)/L2 <= s){
    print(0)
  } 
  
  # bisection method
  lower <- 0
  upper <- max(abs(diff_var)) - 1e-5
  iter <- 1
  
  while(iter <= n_iter && (upper - lower) > crit){
    new_delta <- (lower + upper)/2
    new_diff <- soft_thres(diff_var, new_delta)
    
    if(L1norm(new_diff)/L2norm(new_diff) < s){
      upper <- new_delta
    } else {
      lower <- new_delta
    }
    
    iter <- iter + 1
  }
  
  return(new_delta)
}

# ROBIN Initialization
# Use ROBIN to find the initial centroids by using Local Outliear Factor (LOF)
# avoid using outliers as the centroids
# Key idea: The centroids are the observations located in the most dense region
# and far away from each other, LOF(x)>>1, x is a possible outlier
# if the LOF of observation is approx equal to 1, it is chosen as centroid

# 2 parameters: K - number of clusters; mp - number of neighbors
# X - data matrix, threshold - cuf-off value for LOF
# distance - distance measure
# Follow the ROBIN algorithm to write the function
ROBIN <- function(X, K, mp = 10, threshold = 1.05, distance = "euclidean"){
  # create a list storing the centroids 
  # (actually we store the indices instead of the observations)
  centroids <- rep(0,K)
  
  # calculate the LOF score for each observation
  LOF <- lof(X, k = mp)
  
  # compute distance matrix
  # This distance calculation is not overhead since the first iteration of 
  # kmeans will calculate this anyway
  dist_matrix <- as.matrix(dist(X, method = distance))
  
  # line 1: randomly choose a point as the reference point (it is not a centroid)
  # we use r to find the first centroid
  r <- sample(nrow(X), 1)
  
  # line 2 (we take m = 1 because in R, the first index of vector is 1;
  # in Python, we can follow the algorithm, take m = 0)
  m <- 1
  
  # line 3: choose K centroids
  while (m<=K) {
    if (m == 1){
      # line 5: sort the points in decreasing order of distances from r
      # (we extract the index of observation instead of the observation itself)
      index_decreasing <- sort.int(dist_matrix[r,], decreasing = T, index.return = T)$ix
      
    }else if(m == 2){
      # line 7: sort the points in decreasing order from the minimum distance from 
      # the current centroids
      # (we extract the index of observation instead of the observation itself)
      
      index_decreasing <- sort.int(dist_matrix[centroids,], 
                                   decreasing = T, index.return = T)$ix
    }else{
      
      index_decreasing <- sort.int(apply(dist_matrix[centroids[1:(m-1)],], 2, min),
                                   decreasing = T, index.return = T)$ix
      
    }
    # lines 9 - 14
    # rearrange the LOF following index_decreasing
    sorted.LOF <- LOF[index_decreasing]
    
    # the first point that has LOF < threshold is chosen as the centroid, 
    # which is furthest away from the current centroids
    index <- which(c(sorted.LOF<threshold))[1]
    # store the index of observation that chosen to be centroid
    centroids[m] <- index_decreasing[index] 
    
    # line 15
    m = m + 1
    
  }
  
  # line 17
  return(list(centroids = centroids, LOF = LOF))
}

# Calculate between cluster sum of squares in jth variable
# X - dataset
# cluster - cluster labels 
# K - number of clusters
# varW - variable weights
# obsW - observation weights
BCSS_var <- function(X, cluster, K, varW, obsW){
  
  # matrix of between cluster sum of squares for each variable in each cluster
  B <- matrix(NA, nrow = K, ncol = ncol(X))
  
  # Follow Eq (10)
  # calculate the mean of all weighted observations
  mean_total <- apply(obsW * X, 2, sum) / sum(obsW)
  # weighted total sum of squares (1st term in Eq(10))
  term_total <- obsW * (scale(X, center = mean_total, scale = F))^2
  
  # update between cluster sum of squares for each variable by clusters
  # then take the column summations to get BCSS of each variable
  # since the summation is linear, the final result is BCSS for each variable
  for (i in 1:K){
    # no observation is assigned to this cluster, move to next cluster
    # rare case
    if (sum(cluster == i) == 0){
      next;
    }else{
      # if there are some observations assigned to this cluster,
      # calculate the corresponding within cluster sum of squares
      
      # get the observations assigned to this cluster
      X_cluster <- X[cluster==i, , drop = F] 
      # observation weights in this cluster
      obsW_cluster <- obsW[cluster==i]
      # mean of within cluster
      mean_within <- apply(obsW_cluster * X_cluster, 2, sum) / sum(obsW_cluster, na.rm = T) 
      # weighted within cluster sum of squares (2nd term in Eq (10))
      term_within <- obsW_cluster * scale(X_cluster, center = mean_within, scale = F)^2 
      
      # total sum of squares contributed by observations with cluster label i
      term_t_cluster <- term_total[cluster==i, , drop = F] 
      
      # between-cluster sum of squares contributed by each observation
      # with cluster label i for p variables
      diff <- term_t_cluster - term_within 
      
      # Follow Eq (11)
      # if term_within = 0, no updates in variable weights
      new_varW <- as.vector(sum(varW, na.rm = T) / ((!is.na(term_within)) %*% varW ))
      new_BCSS <- diff * new_varW 
      # between-cluster sum of squares contributed by observations
      # with cluster label i for p variables
      B[i,] <- colSums(new_BCSS, na.rm = T) 
    }
  }
  
  # between-cluster sum of squares for each variable
  B <- colSums(B, na.rm=T) # 1 by p
  return(B)
}

# Apply a weighting function on the clusters detected by ROBIN to reveal outliers
# Follow Eqs (8) and (9)
# Compute observation weights by clusters
# The value is close to 1 indicates this observation is part of the cluster
# X - detected cluster
# mp - number of neighbors
# c - constant parameter
obs_Weights_cluster <- function(X, mp, c = 2){
  # calculate LOF scores for each cluster then standardize, Eq (8) 
  LOF <- lof(X, k = mp)
  LOF_scale <- scale(LOF)
  
  # Eq (9), use the translated bi-weight function
  M <- median(LOF_scale) + mad(LOF_scale)
  v1 <- (1 - ((LOF_scale - M)/(c - M))^2)^2
  v1[LOF >= c] <- 0
  v1[LOF <= M] <- 1
  
  return(v1)
}

# To carry out Phase 1 - Steps 1 and 2
# To calculate observation weights on weighted data
# X - dataset
# K - number of clusters
# mp - number of neighbors
# n_iter - maximum iterations 
# threshold - for identification of outliers, if an observation
#             has weight < threshold, it is considered as an outlier (in noise variable)

wrsk_phase1 <- function(X, K, mp = 10, n_iter = 50, threshold = 0.5){
  
  # dimension of dataset
  n <- nrow(X)
  p <- ncol(X)
  
  # ROBIN - find initial centroids
  # indices
  result <- ROBIN(X, K, mp)
  # observations that chosen to be initial centroids
  centroids <- as.matrix(X[result$centroids, ]) 
  
  # arbitrary vector of initial cluster labels of all observations
  prev_cluster <- rep(1,n)
  
  Iteration <- TRUE
  iter <- 1
  criterion <- 1e+09
  
  while(Iteration){
    
    # arbitrary distance matrix
    dist.centroids <- matrix(NA, nrow = n, ncol = K)
    
    # Step 1 of Phase 1
    # calculate distance matrix of observations from centroids
    # and carry out clusters assignment step
    for (i in 1:K){
      dist.centroids[,i] <- apply(scale(X, center = centroids[i,], 
                                        scale = FALSE)^2,1,sum)
    }
    
    # find the minimum distance of observation from centroid,
    # assign observation to this cluster
    cluster <- apply(dist.centroids, 1, which.min)
    
    # Step 2 of Phase 1
    # compute within-cluster sum of squares
    wcss <- apply(dist.centroids, 1, min)
    
    # update observetion weights and centroids on weighted data
    W <- obs_Weights(X, K, cluster, centroids)
    centroids <- as.matrix(W$centroids) 
    
    # according to the weights, we identify the outliers
    # assign outlier label (0) to observations which have weights < threshold
    new_cluster <- cluster
    new_cluster[W$v <= threshold] <- 0
    
    # Reach the local optima if there is no updates in cluster assignment
    # or stop the algorithm when the max iteration is reached
    if(!identical(prev_cluster, new_cluster) & (iter < n_iter)){
      Iteration <- T
      if(sum(wcss, na.rm = T) < criterion){
        WCSS <- sum(wcss, na.rm = T)
        outliers <- new_cluster
        new_centroids <- W$centroids
        new_weights <- W$v
      }
      iter <- iter + 1
      prev_cluster <- new_cluster
    }else{
      Iteration <- F
      break
    }
    
  }
  
  # clusters - cluster labels of all observations
  # centroids - the list of centroids
  # obsweights - observation weights of all observations
  # outliers - cluster and outlier labels
  # WCSS - within cluster sum of squares
  return(list(clusters = cluster,
              centroids = new_centroids,
              obsweights = new_weights,
              outliers = outliers,
              WCSS = WCSS))
}

# Carry out Phase 1 - Step 3
# Compute the observation weights of all observations
# on weighted or unweighted data
# X - dataset
# K - number of clusters
# cluster - cluster label of each observation
# centroids - centroids of K clusters
obs_Weights <- function(X, K, cluster, centroids){
  # sample size
  n <- nrow(X)
  # arbitrary vector of observation weights 
  v <- c()
  
  # calculate the observation weights of each cluster
  for(i in 1:K){
    # get the observations assigned to cluster i
    X.cluster <- X[cluster == i, , drop = F]
    
    if(nrow(X.cluster) <= 2){
      # only 1 or 2 observations in this cluster,
      # the point itself or the mean is the centroid 
      # observation weights must be 1 to be a part of this cluster
      # so the obs_Weights function is not necessary here
      centroids[i,] <- apply(X.cluster, 2, mean)
      vi <- rep(1,nrow(X.cluster))
    }else{
      # we take mp = 10, however if the number of observations in this cluster
      # is less than 11, we take the number of observations 
      # except the point itself as mp
      # calculate the observation weights by cluster
      n_neighbor <- ifelse(nrow(X.cluster) < 11, nrow(X.cluster) - 1, 10)
      vi <- as.vector(obs_Weights_cluster(X = X.cluster, mp = n_neighbor))
      # compute the next centroids using the observation weights, see Eq (10) 2nd term
      centroids[i,] <- apply(vi * X.cluster, 2, sum) / sum(vi)
    }
    
    v[cluster==i] <- vi
  }
  # v - observation weights of all observations
  # centroids - new centroids
  return(list(v = v, centroids = centroids))
}

# Phase 2 - Step 1
# update variable weights and BCSS
# X - dataset
# K - number of clusters
# s - sparsity parameter
# cluster - cluster labels
# varW - variable weights
# obsW - observation weights
var_Weights <- function(X, K, s, cluster, varW, obsW){
  # BCSS of each variable
  result <- BCSS_var(X, cluster, K, varW, obsW)
  # if BCSS is negative, replace it by 0
  result <- pmax(result,0)
  
  # follow algorithm for Tibshirani's sparse K-means clustering Step 2
  # to update variable weights
  delta <- bisection(result, s)
  # BCSS of p variables
  BCSS <- pmax(result - delta, 0) 
  varW <- BCSS/L2norm(BCSS)
  # update BCSS
  WBCSS <- sum(varW * result)
  
  # varW - updated variable weights
  # WBCSS - updated weighted BCSS
  return(list(varW = varW, WBCSS = WBCSS))
}

# Weighted robust sparse K-means without selection of parameters, K and s
# X - dataset
# K - number of clusters
# s - sparsity parameter
# mp - number of neighbors
# n_iter - maximum number of iteration
# threshold - for identification of outliers, if an observation
#             has weight < threshold, it is considered as an outlier (in noise variable)


WRSK <- function(X, K, s, mp = 10, n_iter = 10, threshold = 0.5){
  # dimension of dataset
  p <- ncol(X)
  n <- nrow(X)
  
  # initial BCSS, choose a very small number 
  prev_ss <- (-1e+09)
  
  # initial variable weights, all = 1/sqrt(p)
  varW <- rep(1/sqrt(p), p) 
  
  # Repeat Phase 1 and Phase 2 
  # Phase 1, Step 1: choose the K centroids using the weighted data
  for(i in 1:n_iter){
    # variable with variable weight = 0 is noise variable
    # exclude the noise variables - dimension reduction
    X.weighted <- t(t(X[ , varW != 0, drop = F]) * sqrt(varW[varW != 0]))
    
    # if i = 1, no updates
    if(i>1){
      # arbitrary matrix to store observations that chosen as centroids
      # K rows - K centroids
      centroids <- matrix(NA, ncol = ncol(X.weighted), nrow = K)
      
      # Find K centroids
      for (j in 1:K){
        # get the observations assigned to cluster i
        X.cluster <- X.weighted[cluster == j, , drop = F] 
        
        # take the mean to get the centroid
        centroids[j,] <- apply(obsW[cluster == j] * X.cluster, 2, sum) / sum(obsW[cluster == j])
        
      }
    }
    
    # Phase 1, Step 2
    # Calculate observation weights on weighted data
    V1 <- wrsk_phase1(X = X.weighted, K = K, threshold = threshold)
    
    # corresponding cluster labels
    cluster <- V1$clusters
    
    # Phase 1, Step 3
    # Calculate observation weights on unweighted data
    V2 <- obs_Weights(X = X, K = K, cluster = cluster,
                      centroids = matrix(NA, ncol = ncol(X), nrow = K))
    
    # follow Eq (12), choose the minimum observation weights to be the 
    # updated observation weights
    obsW <-apply(cbind(V1$obsweights, V2$v), 1, min)
    
    # Phase 2: update variable weights
    B <- var_Weights(X = X, K = K, s = s, 
                     cluster = cluster, varW = varW, obsW = obsW)
    
    # updated variable weights and weighted between-cluster sum of square
    varW <- B$varW
    new_ss <- B$WBCSS
    
    # The first iteration is completed
    # Phase 1 and Phase 2 are repeated until convergence is achieved
    # or the maximum number of iteration is reached
    if(((new_ss - prev_ss)/new_ss) < 1e-8 | i == n_iter){
      
      # convergence, last update
      # also exclude noise variables - dimension reduction
      X.weighted <- t(t(X[ , varW != 0, drop = F]) * sqrt(varW[varW != 0])) 
      
      # arbitrary matrix of centroids 
      centroids <- matrix(NA, ncol = ncol(X.weighted), nrow = K)
      
      # choose K centroids
      for (j in 1:K){
        # observations that assigned to cluster j
        X.cluster <- X.weighted[cluster == j, , drop = F]
        
        # take the mean to get the centroid
        centroids[j,] <- apply(obsW[cluster == j] * X.cluster, 2, sum) /sum(obsW[cluster == j])
      }
      
      # Calculate observation weights on weighted data
      V1 <- wrsk_phase1(X = X.weighted, K = K, threshold = threshold)
      # cluster labels
      cluster <- V1$clusters
      
      # Calculate observation weights on unweighted data
      V2 <- obs_Weights (X = X, K = K, cluster = cluster,
                         centroids = matrix(NA, ncol = ncol(X), nrow = K))
      
      # follow Eq (12), choose the minimum， update observation weights
      obsW <- apply(cbind(V1$obsweights, V2$v), 1 , min)
      
      # identify outliers if the observations have final weights <= threshold,
      # they are declared as outliers
      # outlier label is 0
      cluster_final <- cluster     
      cluster_final[obsW <= threshold] <- 0
      
      break
    }else{
      # repeat the process
      prev_ss <- new_ss
      i <- i + 1
    }
  } 
  
  # clusters - cluster labels of all observations
  # outliers - cluster and outlier labels
  # obsweights - final observation weights
  # varweights - final variable weights
  # WBCSS - final weighted between-cluster sum of squares
  return(list(clusters = cluster,
              outliers = cluster_final,
              obsweights = obsW,
              varweights = varW,
              WBCSS = B$WBCSS))
}

# Ball index (actually this function is same as the Hartigan)
# X - dataset
# K - number of clusters
# s - choices of sparsity parameter, in increasing order, (1,sqrt(q))

wrskB <- function(X, K, s){
  rownames(X) <- 1:nrow(X)
  colnames(X) <- 1:ncol(X)
  n <- nrow(X)

  B_sk <- W_sk <- T_sk <- varW <- c()
  H_sk <- c()
  # arbitrary list storing the results
  result <- list()
  
  for(i in 1:length(s)){
    # use the WRSK function to calculte WBCSS
    result_candidate <- WRSK(X = X, K = K, s = s[i])
    
    # store the results of candidate s
    result[[i]] <- result_candidate
    
    # WBSCC for a given s and K
    B_sk <- result_candidate$WBCSS
    
    # get the number of non-zero variable weights
    varW <- c(varW, length(which(result_candidate$varweights != 0)))
    
    # total sum of squares
    #x_final <- result_candidate$obsweights * X %*% result_candidate$varweights 
    #T_sk <- sum((scale(x_final, scale = F))^2)
    
    # within cluster sum of squares
    #W_sk <- c(W_sk, T_sk - B_sk)
    cluster <- result_candidate$clusters 
    obsW <- result_candidate$obsweights
    term_within <- c()
    for (j in 1:K){
      # no observation is assigned to this cluster, move to next cluster
      # rare case
      if (sum(cluster == j) == 0){
        next;
      }else{
        # if there are some observations assigned to this cluster,
        # calculate the corresponding within cluster sum of squares
        
        # get the observations assigned to this cluster
        X_cluster <- X[cluster==j, , drop = F] 
        # observation weights in this cluster
        obsW_cluster <- obsW[cluster==j]
        # mean of within cluster
        mean_within <- apply(obsW_cluster * X_cluster, 2, sum) / sum(obsW_cluster, na.rm = T) 
        # weighted within cluster sum of squares (2nd term in Eq (10))
        term_within_k <- (obsW_cluster * scale(X_cluster, center = mean_within, scale = F)^2)%*%result_candidate$varweights  
        term_within <- c(term_within,term_within_k) 
        
      }
    }
    W_s <- sum(term_within)
    W_sk <- c(W_sk, W_s)
    
    # stopping criteria:
    # 1. all variable weights are non-zero, i.e. no sparsity
    # the remaining choices of s will return the same results
    # or
    # 2. try all candidates s
    if(varW[i] == ncol(X) | s[i] == tail(s,1)){
      break
      return(list(W_sk = W_sk,
                  varW = varW, s = s[1:i], result = result))
    }
    
  }
  
  # W_sk = WCSS
  # varW - non-zero variable weights
  # s - vector of possible values of s, may not include the remaining choices 
  #     if we reach no sparsity beforehand
  # result - results of WRSK function on given s
  
  return(list(W_sk = W_sk, 
              varW = varW, s = s[1:i], result = result))
  
}

# To get the range of K
# Functions
ball_diff<-function(x){
  n<-length(x)
  Ball <- rep(0,(n-1))
  for(i in 1:(n-1)){
    Ball[i] <- abs(x[i] - x[i+1])
  }
  return(Ball)
}
ball_list <- function(x){
  K <- nrow(x)
  s <- ncol(x)
  diff.matrix <- matrix(0, ncol = s,nrow = (K-1))
  
  for(j in 1:s){
    for(i in 1:(K-1)){
      if(x[i,j]!=0 & x[i+1,j]!=0){
        diff.matrix[i,j] <- abs(x[i,j]-x[i+1,j])
        
      }else{
        diff.matrix[i,j] <-NA
      }
    }
  }
  for(i in 1:(K-1)){
    for(j in 1:s){
      if(is.na(diff.matrix[i,j])){
        diff.matrix[,j] <- NA
        
      }
    }
  }
  
  # assign column names
  colnames(diff.matrix) <- paste0("s", rep(1:s))
  # assign row names
  rownames(diff.matrix) <- paste0("d", rep(1:(K-1)))
  return(diff.matrix)
}
ball_index<-function(x,k){
  K <- nrow(x)
  s <- ncol(x)
  result<-c()
  for(i in 1:s){
    if(is.na(x[1,i])==F){
      index<-which.max(x[,i])
      result[i] <- k[index+1]
    }
  }
  as.numeric(majorityVote(result)$majority)
}

# CH index

# X - dataset
# K - number of clusters
# s - choices of sparsity parameter, in increasing order, (1,sqrt(q))


wrskCH <- function(X, K, s){
  rownames(X) <- 1:nrow(X)
  colnames(X) <- 1:ncol(X)
  n <- nrow(X)

  B_sk <- W_sk <- T_sk <- varW <- c()
  CH_sk <- c()
  # arbitrary list storing the results
  result <- list()
  
  for(i in 1:length(s)){
    # use the WRSK function to calculte WBCSS
    result_candidate <- WRSK(X = X, K = K, s = s[i])
    
    # store the results of candidate s
    result[[i]] <- result_candidate
    
    # WBSCC for a given s and K
    B_sk <- result_candidate$WBCSS
    
    # get the number of non-zero variable weights
    varW <- c(varW, length(which(result_candidate$varweights != 0)))
    
    # within cluster sum of squares
    #W_sk <- c(W_sk, T_sk - B_sk)
    cluster <- result_candidate$clusters 
    obsW <- result_candidate$obsweights
    term_within <- c()
    for (j in 1:K){
      # no observation is assigned to this cluster, move to next cluster
      # rare case
      if (sum(cluster == j) == 0){
        next;
      }else{
        # if there are some observations assigned to this cluster,
        # calculate the corresponding within cluster sum of squares
        
        # get the observations assigned to this cluster
        X_cluster <- X[cluster==j, , drop = F] 
        # observation weights in this cluster
        obsW_cluster <- obsW[cluster==j]
        # mean of within cluster
        mean_within <- apply(obsW_cluster * X_cluster, 2, sum) / sum(obsW_cluster, na.rm = T) 
        # weighted within cluster sum of squares (2nd term in Eq (10))
        term_within_k <- (obsW_cluster * scale(X_cluster, center = mean_within, scale = F)^2)%*%result_candidate$varweights  
        term_within <- c(term_within,term_within_k) 
        
      }
    }
    W_sk <- sum(term_within)
    # calculate the modified CH index
    CH_sk <- (c(CH_sk, (B_sk/(K-1))/(W_sk/(n-K))))
    
    # stopping criteria:
    # 1. all variable weights are non-zero, i.e. no sparsity
    # the remaining choices of s will return the same results
    # or
    # 2. try all candidates s
    if(varW[i] == ncol(X) | s[i] == tail(s,1)){
      break
      return(list(CH = CH_sk,
                  varW = varW, s = s[1:i], result = result))
    }
    
  }
  
  # CH - CH index
  # varW - non-zero variable weights
  # s - vector of possible values of s, may not include the remaining choices 
  #     if we reach no sparsity beforehand
  # result - results of WRSK function on given s
  
  return(list(CH = CH_sk, 
              varW = varW, s = s[1:i], result = result))
  
}

# Gap statistic
# cores - number of cores used

wrskGap_par <- function(X, K, s, n_per = 10,cores=6){
  # for programming convenience
  rownames(X) <- 1:nrow(X)
  colnames(X) <- 1:ncol(X)
  
  # Follow Eq (16) of modified Gap statistic
  # arbitrary vectors of WBCSS,
  # standard error of log(WBCSS_permuted), modified Gap statistic,
  # number of non-zero variable weights 
  B_sk <- se_Bskp <- gap_sk <- varW <- c()
  
  # arbitrary list storing the results
  result <- list()
  
  for(i in 1:length(s)){
    # use the WRSK function to calculte WBCSS
    result_candidate <- WRSK(X = X, K = K, s = s[i])
    
    # store the results of candidate s
    result[[i]] <- result_candidate
    
    # WBSCC for a given s and K
    B_sk <- result_candidate$WBCSS
    
    # get the number of non-zero variable weights 
    varW <- c(varW, length(which(result_candidate$varweights != 0)))
    
    # parallel computing
    result_permuted <- mclapply(1:n_per,
                                function(x, k, s, y){
                                  set.seed(y)
                                  x_permuted <- apply(x, 2, function(xx){sample(xx, length(xx))})
                                  r <- WRSK(X = x_permuted, K = k, s = s)
                                }, x = X, k = K, s = s[i],mc.cores = cores)
    
    
    # arbitrary vector of WBCSS of permuted dataset
    B_sk_p <- c()
    
    for(j in 1:n_per){
      
      # get the WBCSS of the permuted dataset
      B_sk_p[j] <- result_permuted[[j]]$WBCSS
      
    }
    
    # calculate the modified Gap statistic, follow Eq (16)
    gap_sk <- c(gap_sk, log(B_sk) - mean(log(B_sk_p)))
    
    # calculate the standard error
    # see Tibshirani' Gap statistic paper, Step 3
    # se = sd * sqrt(1+1/n_per) 
    n_Bskp  <- length(B_sk_p)
    se_Bskp <- c(se_Bskp, (sqrt(1 - 1/((n_Bskp)^2)) * sqrt(var(log(B_sk_p)))))
    
    # stopping criteria:
    # 1. all variable weights are non-zero, i.e. no sparsity
    # the remaining choices of s will return the same results
    # or
    # 2. try all candidates s
    if(varW[i] == ncol(X) | s[i] == tail(s,1)){
      break
      return(list(Gap = gap_sk, se = se_Bskp,
                  varW = varW, s = s[1:i], result = result))
    }
  }
  
  # Gap - modified Gap statistic
  # se - standard error of WBSCC 
  # varW - non-zero variable weights
  # s - vector of possible values of s, may not include the remaining choices 
  #     if we reach no sparsity beforehand
  # result - results of WRSK function on given s
  
  return(list(Gap = gap_sk, se = se_Bskp,
              varW = varW, s = s[1:i], result = result))
  
}
 # Hartigan
wrskH <- function(X, K, s){
  rownames(X) <- 1:nrow(X)
  colnames(X) <- 1:ncol(X)
  n <- nrow(X)
  
  B_sk <- W_sk <- T_sk <- varW <- c()
  H_sk <- c()
  # arbitrary list storing the results
  result <- list()
  
  for(i in 1:length(s)){
    # use the WRSK function to calculte WBCSS
    result_candidate <- WRSK(X = X, K = K, s = s[i])
    
    # store the results of candidate s
    result[[i]] <- result_candidate
    
    # WBSCC for a given s and K
    B_sk <- result_candidate$WBCSS
    
    # get the number of non-zero variable weights
    varW <- c(varW, length(which(result_candidate$varweights != 0)))
    
    # within cluster sum of squares
    #W_sk <- c(W_sk, T_sk - B_sk)
    cluster <- result_candidate$clusters 
    obsW <- result_candidate$obsweights
    term_within <- c()
    for (j in 1:K){
      # no observation is assigned to this cluster, move to next cluster
      # rare case
      if (sum(cluster == j) == 0){
        next;
      }else{
        # if there are some observations assigned to this cluster,
        # calculate the corresponding within cluster sum of squares
        
        # get the observations assigned to this cluster
        X_cluster <- X[cluster==j, , drop = F] 
        # observation weights in this cluster
        obsW_cluster <- obsW[cluster==j]
        # mean of within cluster
        mean_within <- apply(obsW_cluster * X_cluster, 2, sum) / sum(obsW_cluster, na.rm = T) 
        # weighted within cluster sum of squares (2nd term in Eq (10))
        term_within_k <- (obsW_cluster * scale(X_cluster, center = mean_within, scale = F)^2)%*%result_candidate$varweights  
        term_within <- c(term_within,term_within_k) 
        
      }
    }
    W_s <- sum(term_within)
    W_sk <- c(W_sk, W_s)
    
    # stopping criteria:
    # 1. all variable weights are non-zero, i.e. no sparsity
    # the remaining choices of s will return the same results
    # or
    # 2. try all candidates s
    if(varW[i] == ncol(X) | s[i] == tail(s,1)){
      break
      return(list(W_sk = W_sk,
                  varW = varW, s = s[1:i], result = result))
    }
    
  }
  
  # W_sk -  WCSS
  # varW - non-zero variable weights
  # s - vector of possible values of s, may not include the remaining choices 
  #     if we reach no sparsity beforehand
  # result - results of WRSK function on given s
  
  return(list(W_sk = W_sk, 
              varW = varW, s = s[1:i], result = result))
  
}

Hart.sk<-function(x,k,n = 30){
  lk<-length(k)
  Hartigan<-matrix(NA, ncol = ncol(x), nrow = lk - 1)
  for(i in 1:(lk-1)){
    for(j in 1:ncol(x)){
      if(x[i,j]!=0 & x[i+1,j]!=0){
        Hartigan[i,j] <-(n-k[i]-1)*abs(x[i,j]/x[i+1,j] - 1)
        
      }else{
        Hartigan[i,j] <-NA
      }
    }
  }
  for(i in 1:(lk-1)){
    for(j in 1:ncol(x)){
      if(is.na(Hartigan[i,j])){
        Hartigan[,j] <- NA
        
      }
    }
  }
  Hartigan
}

Hart.sk.p2<-function(x,k){
  n<-ncol(x)
  index<-c()
  hart<-c()
  for(i in 1:n){
    if(is.na(x[1,i])==F){
      index[i] <- which.min(x[,i])
      hart[i]<-x[index[i],i]
    }
  }
  K<-as.numeric(majorityVote(index)$majority)
  s.opt<-which(x[K,]<=10)[1]
  return(list(K=K,s=s.opt))
}

majority.K<-function(gap,ch,h){
  K<-as.matrix(rbind(gap,ch,h))
  n<-ncol(K)
  opt.K<-c()
  for(i in 1:n){
    opt.K[i]<-as.numeric(majorityVote(K[,i])$majority)
  }
  for(i in 1:n){
    if(opt.K[i]==2){
      opt.K[i]<-round(mean(K[,i]))
    }
  }
  opt.K
}
