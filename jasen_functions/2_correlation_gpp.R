# -------------------------------
#
# in the second part of the code, we want to get the correlation operator
# 
# we need the following helper functions:
#
# - cross_prod
#
#----------------------------------

res <- Rmat_diag_full$res
rho_diag <- Rmat_diag_full$rho_diag
NN <- length(patient_sel)
#dmax <- dmax

get_cor_gpp <- function(res=NULL, rho_diag=NULL,NN=NULL,dmax=NULL){
  
  # ------------------------------------------------------------------
  # 
  # goal: return a correlation square matrix that is (p x d times p x d)
  #
  # input:
  #
  # res (list of length equal to number of marks --> each list has 5 items)
  #  - data_i_count
  #  - intensity
  #  - mu
  #  - diag_mat
  #  - gamma
  # 
  # rho_diag
  #
  # NN (integer): number of replicates in our dataset
  #
  # dmax (integer): number of eigenfunctions before we truncate
  #
  # ------------------------------------------------------------------
  
  
  p = length(res) ## number of processes
  cov_est = matrix(0, dmax*p, dmax*p)
  
  for(i in 1:p){
    idx_i = 1:dmax + (i-1)*dmax # iterate in triplets 
    for(j in i:p){
      idx_j = 1:dmax + (j-1)*dmax
      
      if(i==j){
        cov_est[idx_i, idx_j] = rho_diag[[i]]$cov # if i=j, set the 3x3 block equal to the value generated previously
      }else{
        sigma_ij = cross_prod(res[[i]], res[[j]], NN=NN) # G_{i, j}(s, t) in equation 6
        
        # recall that rho_diag was calculated previously
        tmp = t(rho_diag[[i]]$eigenV) %*% sigma_ij %*% rho_diag[[j]]$eigenV
        cov_est[idx_i, idx_j] = tmp
      }
    }
  }
  cov_est = cov_est  + t(cov_est) # save time by using transpose
  diag(cov_est) = diag(cov_est)/2 
  #cor_est = cov2cor(cov_est)
  
  cov_est
}

X <- res[[1]]
Y <- res[[2]]
remove_diag <- T

cross_prod <- function(X=NULL, Y=NULL, NN=NULL,remove_diag = TRUE){
  
  # --------------------------------------------------------------------------
  #
  # goal: calculate G_{i, j}(t,s) in equation 6
  # 
  # input: 
  #  - X             (res[[i]] object): information pertaining to the ith mark
  #  - Y             (res[[j]] object): information pertaining to the jth mark
  #  - NN                    (integer): number of replicates in our simulation
  #  - remove_diag           (boolean): ???
  #
  #
  # --------------------------------------------------------------------------
  
  Sigma = matrix(0, nrow(X$intensity), nrow(X$intensity))  # square matrix for all the time instances
  
  if(is.null(Y)){ #if Y is null, which happens when i = j
    rho_i = X$mu
    if(remove_diag){
      rho_ii = (X$intensity %*% t(X$intensity) - X$diag_mat)/NN  
    }else{
      rho_ii = (X$intensity %*% t(X$intensity))/NN
    }
    
    Sigma_ii = log(rho_ii/(rho_i%*%t(rho_i)))
    Sigma = Sigma_ii
    #print("X")
  }else{
    rho_i = X$mu    # recall mu is a vector, each element corresponds with a time instance
    rho_j = Y$mu
    
    # for marks i and j, find the replicates for which both of them have at least one event
    idx = X$data_i_count$subject_num
    idy = Y$data_i_count$subject_num
    id = intersect(idx, idy)
    
    if(length(id)!=0){ # there's gotta be some replicate with both marks...
      
      # indices 
      matid_x = match(id, X$data_i_count$subject_num)
      matid_y = match(id, Y$data_i_count$subject_num)
      
      
      # rho_ij = 19 x 19 matrix
      rho_ij =  X$intensity[,matid_x,drop=FALSE] %*%  
        t(Y$intensity[,matid_y,drop=FALSE]) / NN 
      
      
      # G_{i, j}(t,s) in equation 6
      Sigma_ij = log(rho_ij/ (rho_i%*%t(rho_j)) )
      Sigma = Sigma_ij
    }
  }
  Sigma 
}