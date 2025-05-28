#
#
# part 1: construct the intensity functions
# 
# parent function:
# - get_rho_diag_pp
#
# helper functions:
# - get_gamma
# - cross_prod
# - get_diag_mat
# - get_intensity
#
# 

#### functions to get correlation matrix ####

pinv_eps = 1e-2

get_gamma <- function(tseq=NULL){
  
  # -----------------------------------------------
  # GOAL: get bandwidth estimate from equation (10)
  #
  # inputs:
  # - tseq (vector):   list of timestamps across all subjects
  #
  # outputs:
  # - res (number): 1 / h_i^2 in equation (10)
  # -----------------------------------------------
  
  n = length(tseq)
  if(n>2& n<10000){ # if num timestamps is under 10K
    res = 1/(sum(dist(tseq))*2/n/(n-1))^2
  }else if(n >= 10000){ # if num timestamps is over 10K, sample 10K of them
    tseq = sample(tseq, 10000)
    n = length(tseq)
    res = 1/(sum(dist(tseq))*2/n/(n-1))^2
  }else{
    res = NaN
  }
  return(res)
}

cross_prod <- function(X=NULL, Y=NULL, NN=NULL,remove_diag = TRUE){
  
  # ---------------------------------------------------------------------------- 
  #
  # goal: estimate \hat{G}_{i,j}(s,t) in equation (13)
  #
  # inputs:
  # - X (list): includes relevant information to the i-th mark
  #
  #   - data_i_count
  #   - intensity (19 x 375 matrix): remember that each replicate has a 19-dim vector
  #   - mu
  #   - diag_mat
  #   - gamma
  # 
  # - NN (integer): number of replicates
  #
  #
  
  Sigma = matrix(0, nrow(X$intensity), nrow(X$intensity)) # matrix of zeros
  
  if(is.null(Y)){ # if we only supply one list 
    
    # 1a) rho_i is already created
    rho_i = X$mu
    
    if(remove_diag){
      
      # rho_ii could be zero when each replicate has at most one event!!
      
      if(max(X$data_i_count$count) <= 1){
        rho_ii <- matrix(0, length(rho_i), length(rho_i))
      } else{
        rho_ii = (X$intensity %*% t(X$intensity) - X$diag_mat)/NN  # sometimes there are errors here
      }
      
    }else{
      rho_ii = (X$intensity %*% t(X$intensity))/NN
    }
    
    Sigma_ii = log(rho_ii/(rho_i%*%t(rho_i)))
    Sigma = Sigma_ii

  }else{
    rho_i = X$mu
    rho_j = Y$mu
    
    idx = X$data_i_count$subject_num
    idy = Y$data_i_count$subject_num
    id = intersect(idx, idy)
    
    if(length(id)!=0){
      
      matid_x = match(id, X$data_i_count$subject_num)
      matid_y = match(id, Y$data_i_count$subject_num)
      
      rho_ij =  X$intensity[,matid_x,drop=FALSE] %*%  
        t(Y$intensity[,matid_y,drop=FALSE]) / NN 
      
      Sigma_ij = log(rho_ij/ (rho_i%*%t(rho_j)) )
      Sigma = Sigma_ij
    }
  }
  Sigma 
}

get_diag_mat <- function(xx=NULL, Tseq=NULL, gamma=NULL,denom=NULL){
  
  # ----------------------------------------------------------------------------
  # GOAL: obtain Gamma_ij^k(s, t) bivariate density estimate  in equation (11)
  #
  # xx      (50-dim vector)  all timestamps for mark i for subject k 
  # Tseq    (19-dim vector)  evenly spaced out times between  0 and 1
  # gamma   (number)         bandwidth parameter
  # denom   (19-dim vector)  w_h(t) denominator for each time in Tseq
  #
  # output:
  # - matrix   (19x19 matrix)  
  # ----------------------------------------------------------------------------
  
  diff = outer(Tseq,xx,'-')            # 19 x 50
  diff = exp(-gamma* (diff^2))/denom   # 19 x 50
  diff%*%t(diff)                       # 19 x 19
}

get_intensity <- function(xx=NULL, Tseq=NULL, gamma=NULL,denom=NULL){
  
  # -----------------------------------------------
  # GOAL: obtain rho_i(t) intensity estimate  in equation (14)
  #
  # the output is w_h(t_i) * sum_{k=1}^n \sum_{l = 1}^{\xi_i^k} kappa()
  # it includes everything in equation (14) except normalizing by 1/n
  #
  #
  # xx      (3805-dim vector)         all timestamps for mark i for subject k 
  # Tseq    (19-dim vector)           evenly spaced out times between  0 and 1
  # gamma   (number)                  bandwidth parameter, 1/h^2
  # denom   (19-dim vector)           w_h(t) denominator for each time in Tseq
  #
  # output:
  # - res   (19-dim vector)           intensity estimate in equation (14) across Tseq times for subject k 
  # -----------------------------------------------
  
  diff = outer(Tseq,xx,'-')       # 19 x 3805
  diff = exp(-gamma* (diff^2))    # 19 x 3805
  res = apply(diff,1,sum)/denom   # 19-dim vec / 19-dim vec
  names(res) = paste("V",1:length(Tseq),sep="")
  res
}

get_rho_diag_pp <- function(data_all=NULL,patient_sel=NULL,feature_sel=NULL,
                            Tseq=NULL, dmax=4, gamma_ratio = 1, 
                            gamma_max = 100, remove_diag = TRUE, ncores=3){

  # ------------------------
  #
  # goal: 
  #
  #
  # inputs:
  # - data_all    (data.frame)   dataframe of 3 columns (subject_num, feature_id, time)
  # - patient_sel (vector)       vector of 1:n to denote patient ID
  # - feature_sel (vector)       vector of 1:p to denote mark ID
  # - Tseq        (vector)       evenly spaced vector from 0 to 1 to denote time
  # - dmax        (integer)
  # - gamma_ratio (integer)
  # - gamma_max   (integer)
  # - remove_diag (bool)
  # - ncores      (integer)
  #  
  #
  #
  # outputs:
  #
  # - 
  #
    

  p = length(feature_sel) ## number of processes
  
  #### estimate intensity ###
  NN = length(patient_sel)
  data_all = data_all[is.element(feature_id, feature_sel),]
  data_all = data_all[is.element(subject_num, patient_sel),]
  
  res_troubleshoot <- list()
  
  for(i in feature_sel){
    
    # 1a) extract info for mark i
    data_i = data_all[feature_id==i,]
    data_i_count = data_i[,.(count=.N),by="subject_num"]
    
    
    # 1b) bandwidth h
    gamma = get_gamma(data_i$time)*gamma_ratio
    gamma[is.na(gamma)] = gamma_max
    
    # 1c) w_h(t) denominator for each Tseq
    denom = sapply(Tseq, function(x){
      truncNorm_denom(x,gamma,0,1)
    })
    
    # 1d) obtain n * rho_i(t) in equation (14) at all Tseq times by subject num
    # thus, this is a (573 x 19)-dim vector
    intensity = data_i[,get_intensity(time,Tseq,gamma,denom),by=c("subject_num")]
    
    
    # 1e) rho_i(t) marginal intensity in equation (12)
    patient_now = unique(intensity$subject_num)                          # patient IDs for mark i that have at least one event (under 400)
    intensity = matrix(intensity$V1,nrow=length(Tseq))                   # rearrange to be a 19 x 380 matrix
    data_i_count = data_i_count[match(patient_now, subject_num),]        # filter out subjects that are not in "patient_now". Rather trivial.
    mu = apply(intensity, 1, sum) / NN                                   # rho_i(t) for mark i 
    
    
    # 1f) rho_{i,j}(s,t) bivariate intensity in equation (12)
    diag_mat = data_i[,get_diag_mat(time,Tseq,gamma,denom),by=c("subject_num")]   # for all k subjects, compute the 19 x 19 matrix entries
    diag_mat = matrix(diag_mat$V1,nrow=length(Tseq)*length(Tseq))                 # rearrange into a 361 x k matrix
    diag_mat = apply(diag_mat, 1, sum)                                            # take the sum across k subjects to return a 361-lengthed vector
    diag_mat = matrix(diag_mat, length(Tseq), length(Tseq))                       # rearrange it to be a 19 x 19 matrix
    
    # 1g) return the following
    # - data_i_count = \xi_i^k counts
    # - intensity    = 19 x n matrix of intensities for all replicates across 19 prespecified timepoints
    # - mu           = average intensity across all replicates for each of the 19 pre-specified timepoints
    # - diag_mat     = 19 x 19 matrix of extraneous terms we will delete later
    # - gamma        = h_i bandwidth parameter
    
    res_troubleshoot[[i]] <- list(data_i_count=data_i_count, intensity=intensity, 
         mu = mu, diag_mat=diag_mat, gamma=gamma)    
  }
  
  res = mclapply(feature_sel, function(i){ # paralellize for each mark
    
    
    # 1a) extract info for mark i
    data_i = data_all[feature_id==i,]
    data_i_count = data_i[,.(count=.N),by="subject_num"]

    
    # 1b) bandwidth h
    gamma = get_gamma(data_i$time)*gamma_ratio
    gamma[is.na(gamma)] = gamma_max

    # 1c) w_h(t) denominator for each Tseq
    denom = sapply(Tseq, function(x){
      truncNorm_denom(x,gamma,0,1)
    })
    
    # 1d) obtain rho_i(t) in equation (14) at all Tseq times by subject num
    # thus, this is a (573 x 19)-dim vector
    intensity = data_i[,get_intensity(time,Tseq,gamma,denom),by=c("subject_num")]
    
  
    # 1e) rho_i(t) marginal intensity
    patient_now = unique(intensity$subject_num)                          # patient IDs for mark i that have at least one event (under 400)
    intensity = matrix(intensity$V1,nrow=length(Tseq))                   # rearrange to be a 19 x 380 matrix
    data_i_count = data_i_count[match(patient_now, subject_num),]       # count table for number of events 
    mu = apply(intensity, 1, sum) / NN                                   # rho_i(t) for mark i 
    

    # 1f) rho_{i,j}(s,t) bivariate intensity in equation (12)
    diag_mat = data_i[,get_diag_mat(time,Tseq,gamma,denom),by=c("subject_num")]   # for all k subjects, compute the 19 x 19 matrix entries
    diag_mat = matrix(diag_mat$V1,nrow=length(Tseq)*length(Tseq))                 # rearrange into a 361 x k matrix
    diag_mat = apply(diag_mat, 1, sum)                                            # take the sum across k subjects to return a 361-lengthed vector
    diag_mat = matrix(diag_mat, length(Tseq), length(Tseq))                       # rearrange it to be a 19 x 19 matrix

    # 1g) return the following
    # - data_i_count = \xi_i^k counts
    # - intensity    = 19 x n matrix of intensities for all replicates across 19 prespecified timepoints
    # - mu           = average intensity across all replicates for each of the 19 pre-specified timepoints
    # - diag_mat     = 19 x 19 matrix 
    # - gamma        = h_i bandwidth parameter
    
    list(data_i_count=data_i_count, intensity=intensity, 
         mu = mu, diag_mat=diag_mat, gamma=gamma)
  }, mc.cores = ncores)
  
  gamma_vec = sapply(res, function(x){
    x$gamma
  })
  
  # 2) approximate C_{X_i, X_j} with Mercer's theorem
  
  rho_diag_troubleshoot <- list()
  
  for(i in 1:length(res)){
    
    xx <- res[[i]]
    
    # 2a) obtain G_{ii}(s,t) 
    Sigma_ii = cross_prod(X=xx,NN=NN,remove_diag=remove_diag)
    
    # 2b) eigendecomposition, and pad the negative eigenvalues as 0
    eigen_res = eigen(Sigma_ii)
    eigen_res$values[ eigen_res$values<0] = 0 
    
    # 2c) keep non-negative eigenvalues (trivially all of them)
    positiveInd = eigen_res$values >= 0
    d = eigen_res$values[positiveInd]
    
    # 2d) keep the eigenvectors associated with non-negative eigenvalues
    eigenV = eigen_res$vectors[, positiveInd, drop=FALSE]
    
    # 2e) cumulative explanatory power of eigenvalues
    FVE = cumsum(d) / sum(d)
    
    # 2f) keep top 3 eigenvalue-eigenvector pairs
    d = d[1:dmax]
    if(any(d < 0)){stop("negative d values!")}
    eigenV = eigenV[, 1:dmax, drop=FALSE]
    
    
    # 2g) variance of KL coefficients???
    cov_diag = t(eigenV) %*% Sigma_ii %*% eigenV
    
    rho_diag_troubleshoot[[i]] <- list(d = d, 
         eigenV = eigenV, 
         cumFVE = FVE, 
         cov = diag(diag(cov_diag),length(d),length(d)), # 2h) only keep the diagonal ones???
         Sigma_ii=Sigma_ii)
    
  }
  
  rho_diag = mclapply(res, function(xx){
    
    Sigma_ii = cross_prod(X=xx,NN=NN,remove_diag=remove_diag)
    ## pca with svd
    eigen_res = eigen(Sigma_ii)
    eigen_res$values[ eigen_res$values<0] = 0 
    positiveInd = eigen_res$values >= 0
    d = eigen_res$values[positiveInd]
    eigenV = eigen_res$vectors[, positiveInd, drop=FALSE]
    FVE = cumsum(d) / sum(d)
    
    d = d[1:dmax]
    if(any(d < 0)){stop("negative d values!")}
    eigenV = eigenV[, 1:dmax, drop=FALSE]
    
    cov_diag = t(eigenV) %*% Sigma_ii %*% eigenV
    list(d = d, 
         eigenV = eigenV, 
         cumFVE = FVE, 
         cov = diag(diag(cov_diag),length(d),length(d)),
         Sigma_ii=Sigma_ii)
  }, mc.cores = ncores)
  
  list(res=res,rho_diag=rho_diag, gamma_vec=gamma_vec)
}

