
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

#### truncated normal denominator for [0,T] ####
truncNorm_denom <- function(t=NULL, gamma=1 ,a=0, b=1){
  
  # ----------------------------------------------------------------------------
  # GOAL: obtain w_h(t) parameter in equation (9)
  #
  # input:
  # - t       (number)  time between 0 and 1 to evalulate w_h(t) at 
  # - gamma   (number)  bandwidth parameter, 1/h_i^2
  # - a       (number)  lower integral limit 
  # - b       (number)  upper integral limit
  #
  # output:
  # - integral to obtain 1/w_h(t)
  # ----------------------------------------------------------------------------
  
  integrate(function(x) {exp(-gamma*(t-x)^2)}, a, b)$value
}


truncNorm_denom_v2 <- function(t=NULL, gamma=1 ,a=0, b=1){
  return(pnorm((b - t)/gamma) - pnorm((a - t)/gamma))
}

truncNorm_denom_v3 <- function(t=NULL, gamma=1 ,a=0, b=1){
  (1 / sqrt(2 * pi * gamma^2)) * integrate(function(x) {exp(-(t-x)^2 / (2 * gamma^2))}, a, b)$value
}

#### function to estimate intensity
get_intensity <- function(xx=NULL, Tseq=NULL, gamma=NULL,denom=NULL){
  
  # -----------------------------------------------
  # GOAL: obtain rho_i(t) intensity estimate  in equation (12)
  #
  # xx      (vector)  all timestamps for mark i for subject k 
  # Tseq    (vector)  evenly spaced out times between  0 and 1
  # gamma   (number)  bandwidth parameter, 1/h^2
  # denom   (vector)  w_h(t) denominator for each time in Tseq
  #
  # output:
  # - res   (vector)  intensity estimate in equation (12) across Tseq times for subject k 
  # -----------------------------------------------
  
  diff = outer(Tseq,xx,'-')
  diff = exp(-gamma* (diff^2))
  res = apply(diff,1,sum)/denom
  names(res) = paste("V",1:length(Tseq),sep="")
  #data.table(t(res))
  res
}

get_intensity_v2 <- function(xx=NULL, Tseq=NULL, gamma=NULL,denom=NULL){
  
  diff = outer(Tseq,xx,'-') # 19 x n matrix
  diff = (1 / sqrt(2 * pi * gamma^2)) * exp(- (diff^2) / (2 * gamma^2)) # kernel density 
  
  res = apply(diff,1,mean)/denom # take the mean not the sum!! for the xi_i^k inverse term.
  names(res) = paste("V",1:length(Tseq),sep="")
  #data.table(t(res))
  res
}

### function to diagomal elements matrix
get_diag_mat <- function(xx=NULL, Tseq=NULL, gamma=NULL,denom=NULL){
  
  # ----------------------------------------------------------------------------
  # GOAL: obtain Gamma_ij^k(s, t) bivariate density estimate  in equation (11)
  #
  # xx      (vector)  all timestamps for mark i for subject k 
  # Tseq    (vector)  evenly spaced out times between  0 and 1
  # gamma   (number)  bandwidth parameter
  # denom   (vector)  w_h(t) denominator for each time in Tseq
  #
  # output:
  # - res   (vector)  density estimation in equation (9) across Tseq times
  # ----------------------------------------------------------------------------
  
  diff = outer(Tseq,xx,'-')
  diff = exp(-gamma* (diff^2))/denom
  diff%*%t(diff)
}

# 
# data_all <- data_df
# dmax <- 4
# gamma_ratio <- 1
# gamma_max <- 100
# remove_diag <- T
# ncores <- 3

get_rho_diag_pp <- function(data_all=NULL,patient_sel=NULL,feature_sel=NULL,
                            Tseq=NULL, dmax=4, gamma_ratio = 1, 
                            gamma_max = 100, remove_diag = TRUE, ncores=3){
  
  # ------------------------
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
  
  # Tmax = 1 ## define domain [0,T]
  # M = 50 ## t_1,...,t_M
  # Tseq = seq(0.1,Tmax-0.1,length=M)
  p = length(feature_sel) ## number of processes
  
  #### estimate intensity ###
  NN = length(patient_sel)
  data_all = data_all[is.element(feature_id, feature_sel),]
  data_all = data_all[is.element(subject_num, patient_sel),]
  
  res = mclapply(feature_sel, function(i){
    
    # for each mark, parallelize this computation
    
    # extract info for mark i
    data_i = data_all[feature_id==i,]
    data_i_count = data_i[,.(count=.N),by="subject_num"]

    
    # get bandwidth h
    gamma = get_gamma(data_i$time)*gamma_ratio
    gamma[is.na(gamma)] = gamma_max # if the get_gamma() function returned NA, fill it in with gamma_max

    
    # get w_h(t) denominator term for each Tseq
    denom = sapply(Tseq, function(x){
      truncNorm_denom(x,gamma,0,1)
    })
    denom_v2 = sapply(Tseq, function(x){
      truncNorm_denom_v2(x,gamma,0,1)
    })    
    denom_v3 = sapply(Tseq, function(x){
      truncNorm_denom_v3(x,gamma,0,1)
    })    
    
    # obtain intensities for all k replicates at all Tseq times
    intensity    = data_i[,get_intensity(time,Tseq,gamma,denom),by=c("subject_num")] # for each subject with non-zero amount of events of mark i, compute the intensity at all Tseq (19) points
    intensity_v2 = data_i[,get_intensity_v2(time,Tseq,gamma,denom),by=c("subject_num")]
    
    
    
    patient_now = unique(intensity$subject_num)                          # patient IDs for mark i that have at least one event (under 400)
    intensity = matrix(intensity$V1,nrow=length(Tseq))                   # rearrange to be a 19 x 380 matrix
    data_i_count = data_i_count[match(patient_now, subject_num),]        # count table for number of events 
    mu = apply(intensity, 1, sum) / NN                                   # rho_i(t) for mark i 
    

    # now, obtain bivariate intensity in equation (12)
    diag_mat = data_i[,get_diag_mat(time,Tseq,gamma,denom),by=c("subject_num")]   # for all k subjects, compute the 19 x 19 matrix entries
    diag_mat = matrix(diag_mat$V1,nrow=length(Tseq)*length(Tseq))                 # rearrange into a 361 x k matrix
    diag_mat = apply(diag_mat, 1, sum)                                            # take the sum across k subjects to return a 361-lengthed vector
    diag_mat = matrix(diag_mat, length(Tseq), length(Tseq))                       # rearrange it to be a 19 x 19 matrix
    # we now have a matrix that is 

    
    # return the following:
    # data_i_count = xi_i^k counts
    # intensity = 19 x n matrix of intensities for all replicates across 19 prespecified timepoints
    # mu = average intensity across all replicates for each of the 19 pre-specified timepoints
    # diag_mat = 19 x 19 matrix 
    # gamma = h_i bandwidth parameter
    list(data_i_count=data_i_count, intensity=intensity, 
         mu = mu, diag_mat=diag_mat, gamma=gamma)
  }, mc.cores = ncores)
  
  gamma_vec = sapply(res, function(x){
    x$gamma
  })
  
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
    list(d = d, eigenV = eigenV, cumFVE = FVE, 
         cov = diag(diag(cov_diag),length(d),length(d)),Sigma_ii=Sigma_ii)
  }, mc.cores = ncores)
  
  list(res=res,rho_diag=rho_diag, gamma_vec=gamma_vec)
}


#### functions to get correlation matrix ####

cross_prod <- function(X=NULL, Y=NULL, NN=NULL,remove_diag = TRUE){
  
  # ---------------------------------------------------------------------------- 
  # inputs:
  # - X (list)
  #
  #
  #
  
  Sigma = matrix(0, nrow(X$intensity), nrow(X$intensity))
  
  if(is.null(Y)){
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