# ------------------------------------------
#
#
# finally, we use thresholding to determine the graph
#
# main function:
#
# - run_GPP_HT_BIC
#
#
# helper functions:
# 
# - get_lamseq
#   - get_groupNorm
# - Cov_hardT
# - find_min_eigen
# - calc_BIC_1
#   - get_nonzero
#
#
# ------------------------------------------

get_nonzero <- function(S=NULL, grpind=NULL){
  
  # ------------------------------------------------------------------
  #
  # goal:
  # - calculate q_{\tau_c, \tau_p} in equation (20) for GIC
  # - how many entries that are not on the main block diagonal, are non-zero?
  #
  # input:
  # - Rinv (80 x 80 matrix): candidate \hat{\theta} precision matrix
  # - grpind (matrix): matrix to indicate dimensions of block matrix
  #
  # --------------------------------------------------------------------  
  
  # for each diagonal block, set it equal to zero 
  for(i in 1:nrow(grpind) ){
    S[grpind[i,1]:grpind[i,2],grpind[i,1]:grpind[i,2]] = 0
  }
  
  sum(abs(S[lower.tri(S)]) !=0)
}

calc_BIC_1 <- function(Rinv=NULL, Rmat=NULL, n=NULL,grpind=NULL,factor=NULL){
  
  # ------------------------------------------------------------------
  #
  # goal:
  # - calculate the GIC 
  # - determinant = product of eigenvalues
  # - logdet = sum of logs of eigenvalues
  #
  # input:
  # - Rinv (80 x 80 matrix):   inverse matrix candidate = \hat{\theta}^{d, \tau_c, \tau_p}
  # - Rmat (80 x 80 matrix):   original Rmat to compare to
  # - grpind       (matrix):   start and end indices for all p marks
  # - factor       (number):   usually square root of ntrain
  #
  # --------------------------------------------------------------------
  
  # all 80 eigenvalues
  ss = eigen(Rinv)$values  
  
  # logdet(candidate) = sum of log of eigenvalues that are sufficiently positive 
  val = sum(log(ss[ss>pinv_eps]))
  
  # n * [ tr(AB) - logdet(\hat_theta) ]
  loglik = n*(sum(diag((Rinv%*%Rmat))) - val)
  df = get_nonzero(Rinv, grpind)
  loglik + factor*df
}

find_min_eigen <- function(S=NULL){
  
  # ----------------------------------
  #
  # input is a list of matrices
  #
  # return a list of min eigenvalues for each matrix
  #
  # ---------------------------
  eig = sapply(S, function(s){
    min(eigen(s)$val)
  })
  eig
}

get_groupNorm <- function(S=NULL, grpind=NULL){
  
  # --------------------------------------------------------------------------
  #
  # goal:
  # - Recall that the Rmat matrix is still augmented 
  # - We want to create a pxp matrix where each entry is the norm of the augmented elements
  #
  # example:
  # - [A_11 A_12        becomes     [norm(A_11) norm(A_12)
  #   [A_21 A_22]                    norm(A_21) norm(A_22)] 
  #
  # input: 
  # - S (matrix): the Rmat matrix
  # - grpind (matrix): the grpind matrix 
  #
  # output:
  # - a (p x p) matrix that seems to be taking the norm of each block matrix from Rmat
  # - it seems like we are taking the rmse: sqrt(sum(elements^2) / num(elements))
  #
  # --------------------------------------------------------------------------
  
  res = matrix(0, nrow(grpind), nrow(grpind)) # p x p
  
  for(i in 1:nrow(grpind) ){
    for(j in i:nrow(grpind) ){
      tmp = S[grpind[i,1]:grpind[i,2],
              grpind[j,1]:grpind[j,2], drop=FALSE] # obtain the relevant slice in the Rmat matrix
      
      res[i,j] = sum(tmp^2)/ length(tmp) # average of squared entries??? why???
    }
  }
  
  res = res + t(res)  # this is quicker, but we double counted the diagonals 
  diag(res) = diag(res)/2 # so we just divide by 2
  res = sqrt(res)         # square root everything
  res
}

get_lamseq <- function(Rmat=NULL, grpind=NULL, len=100){
  
  # ------------------------------------------------------------------
  #
  # goal:
  # - Helper function to extract off-diagonal entries of the F-normed Rmat
  #
  # input:
  # - Rmat (80 x 80 matrix): truncated from the original p x d matrix if enough eigenfunctions could explain 
  # - grpind (matrix): start and end indices for all p marks
  # - ntrain (integer): number of replicates
  # - factor (number): usually square root of ntrain
  #
  # --------------------------------------------------------------------
  
  # 1) take the norm
  Rmat_norm = get_groupNorm(Rmat,grpind)
  
  # 2) find all unique off-diagonal terms
  lam1seq  = sort(unique(Rmat_norm[upper.tri(Rmat_norm)]), decreasing = T)
  
  # 3) if there are more than 100 unique ones, take the quantile to get 100 pseudo-off-diagonal terms
  if(length(lam1seq)>100){
    lam1seq = quantile(lam1seq,seq(1,0,length=100))
  }
  lam1seq
}

Cov_hardT <- function(S=NULL, grpind=NULL, lamseq=NULL){
  
  # ------------------------------------------------------------------
  #
  # goal:
  # - Generate Rmat matrices of differing sparsity edits. 
  #
  # input:
  # - S (80 x 80 matrix): It is just Rmat. It is truncated from the original p x d square matrix if enough eigenfunctions could explain 
  # - grpind (matrix): start and end indices for all p marks
  # - lamseq (vector): obtained from the get_lamseq() function. 
  #
  # --------------------------------------------------------------------  
  
  Grp_norm = get_groupNorm(S, grpind)
  omega_arr = array(0, c(nrow(S), ncol(S), length(lamseq))) # 80 x 80 x 100
  
  for(i in seq_along(lamseq)){ # for i = 1, 2, ... 100
    lam1 = lamseq[i]
    X = S
    
    for(j in 1:(nrow(Grp_norm)-1)){ # j = 1, 2, ... 39
      for( k in (j+1):nrow(Grp_norm)){ # k = j+1, j+2, ... 40
        if(Grp_norm[j,k] <= lam1){   # if the triu elements <= lam1:
          X[grpind[j,1]:grpind[j,2], 
            grpind[k,1]:grpind[k,2]] = 0
          X[grpind[k,1]:grpind[k,2],
            grpind[j,1]:grpind[j,2]] = 0
        }
      }
    }
    #diag(X) = diag(X)/2
    omega_arr[,,i] = X
  }
  
  omega_arr
}

Rmat <- Rmat_IC

run_GPP_HT_BIC_troubleshoot <- function(Rmat=NULL, grpind=NULL,ntrain=NULL, factor=NULL, tol = 1e-8){
  
  # ------------------------------------------------------------------
  #
  # only use for loops to figure out where the issue lies
  #
  
  # 1) get the off-diagonals of the normed RMat
  lamseq = get_lamseq(Rmat, grpind)
  
  # 2) construct the omega_arr matrix of increasing sparsity
  Rmat_thre_1 = Cov_hardT(S=Rmat, grpind=grpind, lamseq=lamseq)
  
  print("Are all matrices symmetric in Rmat_thre_1?") # check if all of them are symmetric
  print(apply(Rmat_thre_1, 3, isSymmetric) %>% all())

  BIC_list_troubleshoot <- list()
  
  for(i in 1:dim(Rmat_thre_1)[3]){
    Rmat_1st  = Rmat_thre_1[,,i] # the i-th 80 x 80 matrix
    
    min_eig_val = min(find_min_eigen(c(list(Rmat_1st))))
    
    if(min_eig_val < pinv_eps){ # if the min eigenvalu is negative (aka we need to fix it), do the adjusting in part 3
      min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
      Rmat_1st = adjust_S(Rmat_1st,min_eig_val)
    }
    
    # 4) then take the moore-penrose pseudo-inverse
    Rinv = pinv(Rmat_1st,pinv_eps)
    lamseq_2 = get_lamseq(Rinv, grpind) # lamseq the inv
    Rinv_thre_2 = Cov_hardT(S=Rinv, grpind=grpind, lamseq=lamseq_2) # construct omega_arr matrix of increasing sparsity again
    
    # troubleshoot: let all values below tolerance be equal to 0
    
    Rinv_thre_2[abs(Rinv_thre_2) < tol] <- 0
    
    print(i)
    print("Are all matrices symmetric in Rmat_thre_2?") # check if all of them are symmetric
    print(apply(Rinv_thre_2, 3, isSymmetric) %>% all())    
    
    # 5) for each varying level of sparsity of the pinv:
    
    BIC_val_troubleshoot <- c()
    
    for(j in 1:dim(Rinv_thre_2)[3]){
      
      x <- Rinv_thre_2[,,j]
      min_eig_val = min(find_min_eigen(c(list(x))))
      if(min_eig_val < pinv_eps){
        min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
        diag(x) = diag(x) + min_eig_val
      }
      
      
      BIC <- calc_BIC_1(x, Rmat, ntrain, grpind, factor) 
      BIC_val_troubleshoot <- c(BIC_val_troubleshoot, BIC)
    }
    
    
    BIC_list_troubleshoot[[i]] <- BIC_val_troubleshoot    
  }
  
  
  
  # 6) we have a list of vectors. find the minimum BIC for each sparsity level
  BIC_list_1 = sapply(BIC_list_troubleshoot, function(x){min(x)})
  min_ind_1 = which.min(BIC_list_1) # which sparsity level got us the least BIC
  mid_ind_2 = which.min(BIC_list_troubleshoot[[min_ind_1]]) # which pinv level got us the least BIC
  
  
  # 7) after obtaining the indices that got us the overall least BIC, we recreate the matrix:
  
  # 7a) extract sparsity i
  Rmat_1st  = Rmat_thre_1[,,min_ind_1]
  
  # 7b) make pd if needed
  min_eig_val = min(find_min_eigen(c(list(Rmat_1st))))
  if(min_eig_val < pinv_eps){
    min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
    Rmat_1st = adjust_S(Rmat_1st,min_eig_val)
  }
  
  # 7c) get pseudo-inverse 
  Rinv = pinv(Rmat_1st,pinv_eps)
  
  # 7d) keep the j-th pseudo-inverse sparsity level
  lamseq_2 = get_lamseq(Rinv, grpind)
  Rinv_thre_2 = Cov_hardT(S=Rinv, grpind=grpind, lamseq=lamseq_2[mid_ind_2]) # 80 x 80 x 1
  
  res = Rinv_thre_2[,,1] # 80 x 80
  
  # 7e) make psd
  min_eig_val = min(find_min_eigen(c(list(res))))
  if(min_eig_val < pinv_eps){
    min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
    diag(res) = diag(res) + min_eig_val
  }
  res
}

run_GPP_HT_BIC <- function(Rmat=NULL, grpind=NULL,ntrain=NULL, factor=NULL){
  
  # ------------------------------------------------------------------
  #
  # goal:
  # - vary the level of sparsity of correlation matrix via lamseq
  # - take the pseudoinverse, then vary its sparsity
  # - sparsify --> pseudoinverse --> sparsify --> compare GIC wrt original Rmat
  #
  # input:
  # - Rmat (80 x 80 matrix): truncated from the original p x d matrix if enough eigenfunctions could explain 
  # - grpind (matrix): start and end indices for all p marks
  # - ntrain (integer): number of replicates
  # - factor (number): usually square root of ntrain
  #
  #
  
  # 1) get the off-diagonals of the normed RMat
  lamseq = get_lamseq(Rmat, grpind)
  
  # 2) construct the omega_arr matrix of increasing sparsity
  Rmat_thre_1 = Cov_hardT(S=Rmat, grpind=grpind, lamseq=lamseq)
  
  
  BIC_list = lapply(1:dim(Rmat_thre_1)[3], function(i){
    
    # 3) for each 80 x 80 matrix, first ensure that it is properly positive definite
    Rmat_1st  = Rmat_thre_1[,,i] # the i-th 80 x 80 matrix
    
    min_eig_val = min(find_min_eigen(c(list(Rmat_1st))))
    
    if(min_eig_val < pinv_eps){ # if the min eigenvalu is negative (aka we need to fix it), do the adjusting in part 3
      min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
      Rmat_1st = adjust_S(Rmat_1st,min_eig_val)
    }
    
    # 4) then take the moore-penrose pseudo-inverse
    Rinv = pinv(Rmat_1st,pinv_eps)
    lamseq_2 = get_lamseq(Rinv, grpind) # lamseq the inv
    Rinv_thre_2 = Cov_hardT(S=Rinv, grpind=grpind, lamseq=lamseq_2) # construct omega_arr matrix of increasing sparsity again
    
    # added the tolerance factor
    Rinv_thre_2[abs(Rinv_thre_2) < tol] <- 0
    
    # 5) for each varying level of sparsity of the pinv:
    BIC_val = apply(Rinv_thre_2, 3, function(x){
      
      # 5a) ensure its sufficiently pd
      min_eig_val = min(find_min_eigen(c(list(x))))
      if(min_eig_val < pinv_eps){
        min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
        diag(x) = diag(x) + min_eig_val
      }
      
      # 5b) calculate the BIC wrt the original 80 x 80 Rmat and the 80 x 80 current pinv
      calc_BIC_1(x, Rmat, ntrain, grpind, factor)
    })
    BIC_val
  })
  
  # 6) we have a list of vectors. find the minimum BIC for each sparsity level
  BIC_list_1 = sapply(BIC_list, function(x){min(x)})
  min_ind_1 = which.min(BIC_list_1) # which sparsity level got us the least BIC
  mid_ind_2 = which.min(BIC_list[[min_ind_1]]) # which pinv level got us the least BIC
  
  
  # 7) after obtaining the indices that got us the overall least BIC, we recreate the matrix:
  
  # 7a) extract sparsity i
  Rmat_1st  = Rmat_thre_1[,,min_ind_1]
  
  # 7b) make pd if needed
  min_eig_val = min(find_min_eigen(c(list(Rmat_1st))))
  if(min_eig_val < pinv_eps){
    min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
    Rmat_1st = adjust_S(Rmat_1st,min_eig_val)
  }
  
  # 7c) get pseudo-inverse 
  Rinv = pinv(Rmat_1st,pinv_eps)
  
  # 7d) keep the j-th pseudo-inverse sparsity level
  lamseq_2 = get_lamseq(Rinv, grpind)
  Rinv_thre_2 = Cov_hardT(S=Rinv, grpind=grpind, lamseq=lamseq_2[mid_ind_2]) # 80 x 80 x 1
  
  res = Rinv_thre_2[,,1] # 80 x 80
  
  # 7e) make psd
  min_eig_val = min(find_min_eigen(c(list(res))))
  if(min_eig_val < pinv_eps){
    min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
    diag(res) = diag(res) + min_eig_val
  }
  res
}

factor = sqrt(ntrain)
res_GPP_BIC = run_GPP_HT_BIC(Rmat_IC, grpind, ntrain,factor)
graph_all[["GPP_BIC"]] = as.matrix((get_groupNorm(res_GPP_BIC, grpind)!=0)+0)
res_all[["GPP_BIC"]] = res_GPP_BIC


