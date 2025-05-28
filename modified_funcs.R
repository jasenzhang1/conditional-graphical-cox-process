
# i want this to return a graph with at least q edges with the minimum GIC

library(pbmcapply)

# min dges
# tolerance
# mclapply

pinv_eps = 1e-2

run_GPP_HT_BIC_v2 <- function(Rmat=NULL, grpind=NULL,ntrain=NULL, factor=NULL, num_edges=0, tol = 1e-8, ncores = 1){
  
  lamseq = get_lamseq(Rmat, grpind)
  Rmat_thre_1 = Cov_hardT(S=Rmat, grpind=grpind, lamseq=lamseq)

  BIC_list = pbmclapply(1:5, function(i){  
  # BIC_list = pbmclapply(1:dim(Rmat_thre_1)[3], function(i){
    Rmat_1st  = Rmat_thre_1[,,i]
    min_eig_val = min(find_min_eigen(c(list(Rmat_1st))))
    if(min_eig_val < pinv_eps){
      min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
      Rmat_1st = adjust_S(Rmat_1st,min_eig_val)
    }
    
    Rinv = pinv(Rmat_1st,pinv_eps)
    
    if(! isSymmetric(Rinv)){
      Rinv = 0.5 * (Rinv + t(Rinv))
    }    
    
    lamseq_2 = get_lamseq(Rinv, grpind)
    Rinv_thre_2 = Cov_hardT(S=Rinv, grpind=grpind, lamseq=lamseq_2)
    
    # thresholding
    Rinv_thre_2[abs(Rinv_thre_2) < tol] <- 0
    
    BIC_val = apply(Rinv_thre_2, 3, function(x){
      
      min_eig_val = min(find_min_eigen(c(list(x))))
      if(min_eig_val < pinv_eps){
        min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
        diag(x) = diag(x) + min_eig_val
      }
      
      # modify this
      BIC_i <- calc_BIC_1(x, Rmat, ntrain, grpind, factor)
      graph_i <- as.matrix((get_groupNorm(x, grpind)!=0)+0)
      num_edges_i <- sum( graph_i[upper.tri(graph_i)] )
        
      c(BIC_i, num_edges_i)
    })
    BIC_val
  }, mc.cores = ncores)
  
  # filter the BIC's with at least 10 edges and then remove # edges 
  

  BIC_list_1 = sapply(BIC_list, function(x){withCallingHandlers({
    min(x[1, x[2,] >= num_edges])
  },
  warning = function(w) {
    if (grepl("no non-missing arguments to min; returning Inf", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  }
  )}
  )
  
  min_ind_1 = which.min(BIC_list_1)
  
  x2 <- BIC_list[[min_ind_1]]
  
  cond <- x2[2, ] >= num_edges
  idx <- which(cond)            
  mid_ind_2 <- idx[which.min(x2[1,][cond])] 
  
  
  Rmat_1st  = Rmat_thre_1[,,min_ind_1]
  min_eig_val = min(find_min_eigen(c(list(Rmat_1st))))
  if(min_eig_val < pinv_eps){
    min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
    Rmat_1st = adjust_S(Rmat_1st,min_eig_val)
  }
  
  Rinv = pinv(Rmat_1st,pinv_eps)
  
  if(! isSymmetric(Rinv)){
    Rinv = 0.5 * (Rinv + t(Rinv))
    print('symmetrizing...')
    print(isSymmetric(Rinv))
  }  
  
  lamseq_2 = get_lamseq(Rinv, grpind)
  Rinv_thre_2 = Cov_hardT(S=Rinv, grpind=grpind, lamseq=lamseq_2[mid_ind_2])
  
  res = Rinv_thre_2[,,1]
  min_eig_val = min(find_min_eigen(c(list(res))))
  if(min_eig_val < pinv_eps){
    min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
    diag(res) = diag(res) + min_eig_val
  }
  res
}

# modify to add tolerance

run_GPP_HT_BIC_tol <- function(Rmat=NULL, grpind=NULL,ntrain=NULL, factor=NULL, tol = 1e-8){
  
  lamseq = get_lamseq(Rmat, grpind)
  Rmat_thre_1 = Cov_hardT(S=Rmat, grpind=grpind, lamseq=lamseq)

  print("Are all matrices symmetric in Rmat_thre_1?") # check if all of them are symmetric
  print(apply(Rmat_thre_1, 3, isSymmetric) %>% all())    
  
  BIC_list = lapply(1:dim(Rmat_thre_1)[3], function(i){
    print(i)
    Rmat_1st  = Rmat_thre_1[,,i]
    min_eig_val = min(find_min_eigen(c(list(Rmat_1st))))
    if(min_eig_val < pinv_eps){
      min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
      Rmat_1st = adjust_S(Rmat_1st,min_eig_val)
    }
    
    Rinv = pinv(Rmat_1st,pinv_eps)
    lamseq_2 = get_lamseq(Rinv, grpind)
    Rinv_thre_2 = Cov_hardT(S=Rinv, grpind=grpind, lamseq=lamseq_2)
    
    # thresholding
    Rinv_thre_2[abs(Rinv_thre_2) < tol] <- 0 
    
    BIC_val = apply(Rinv_thre_2, 3, function(x){
      
      min_eig_val = min(find_min_eigen(c(list(x))))
      if(min_eig_val < pinv_eps){
        min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
        diag(x) = diag(x) + min_eig_val
      }
      
      calc_BIC_1(x, Rmat, ntrain, grpind, factor)
    })
    BIC_val
  })
  
  BIC_list_1 = sapply(BIC_list, function(x){min(x)})
  min_ind_1 = which.min(BIC_list_1)
  mid_ind_2 = which.min(BIC_list[[min_ind_1]])
  
  
  Rmat_1st  = Rmat_thre_1[,,min_ind_1]
  min_eig_val = min(find_min_eigen(c(list(Rmat_1st))))
  if(min_eig_val < pinv_eps){
    min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
    Rmat_1st = adjust_S(Rmat_1st,min_eig_val)
  }
  
  Rinv = pinv(Rmat_1st,pinv_eps)
  lamseq_2 = get_lamseq(Rinv, grpind)
  Rinv_thre_2 = Cov_hardT(S=Rinv, grpind=grpind, lamseq=lamseq_2[mid_ind_2])
  
  res = Rinv_thre_2[,,1]
  min_eig_val = min(find_min_eigen(c(list(res))))
  if(min_eig_val < pinv_eps){
    min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
    diag(res) = diag(res) + min_eig_val
  }
  res
}


run_GPP_HT_BIC_troubleshoot <- function(Rmat=NULL, grpind=NULL,ntrain=NULL, factor=NULL, num_edges = 0, tol = 1e-8){
  
  # ------------------------------------------------------------------
  #
  # only use for loops to figure out where the issue lies
  #
  
  print('Is Rmat symmetric?')
  print(isSymmetric(Rmat))
  
  # 1) get the off-diagonals of the normed RMat
  lamseq = get_lamseq(Rmat, grpind)
  
  # 2) construct the omega_arr matrix of increasing sparsity
  Rmat_thre_1 = Cov_hardT(S=Rmat, grpind=grpind, lamseq=lamseq)
  
  print("Are all matrices symmetric in Rmat_thre_1?") # check if all of them are symmetric
  print(apply(Rmat_thre_1, 3, isSymmetric) %>% all())
  
  BIC_list_troubleshoot <- list()
  
  #for(i in 1:dim(Rmat_thre_1)[3]){
  for(i in 1:5){
    print(i)
    
    Rmat_1st  = Rmat_thre_1[,,i] # the i-th 80 x 80 matrix
    
    min_eig_val = min(find_min_eigen(c(list(Rmat_1st))))
    
    if(min_eig_val < pinv_eps){ # if the min eigenvalu is negative (aka we need to fix it), do the adjusting in part 3
      min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
      Rmat_1st = adjust_S(Rmat_1st,min_eig_val)
    }
    
    # 4) then take the moore-penrose pseudo-inverse
    Rinv = pinv(Rmat_1st,pinv_eps)
    
    print('Performed inversion, checking if pinv is symmetric')
    print(isSymmetric(Rinv))
    
    if(! isSymmetric(Rinv)){
      Rinv = 0.5 * (Rinv + t(Rinv))
      print('symmetrizing...')
      print(isSymmetric(Rinv))
    }
    
    lamseq_2 = get_lamseq(Rinv, grpind) # lamseq the inv
    Rinv_thre_2 = Cov_hardT(S=Rinv, grpind=grpind, lamseq=lamseq_2) # construct omega_arr matrix of increasing sparsity again
    
    # troubleshoot: let all values below tolerance be equal to 0
    
    print("Are all matrices symmetric in Rmat_thre_2 before applying tolerance?") # check if all of them are symmetric
    print(apply(Rinv_thre_2, 3, isSymmetric) %>% all())       
    
    # thresholding
    Rinv_thre_2[abs(Rinv_thre_2) < tol] <- 0
    

    print("Are all matrices symmetric in Rmat_thre_2 after applying tolerance?") # check if all of them are symmetric
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
      
      # obtain BIC + number of edges if we map it 
      BIC_i <- calc_BIC_1(x, Rmat, ntrain, grpind, factor)
      graph_i <- as.matrix((get_groupNorm(x, grpind)!=0)+0)
      num_edges_i <- sum( graph_i[upper.tri(graph_i)] )      
      

      BIC_val_troubleshoot <- cbind(BIC_val_troubleshoot, c(BIC_i, num_edges_i))
    }
    
    
    BIC_list_troubleshoot[[i]] <- BIC_val_troubleshoot    
  }
  

  
  # 6) we have a list of vectors. find the minimum BIC for each sparsity level
  
  print('looking at BIC_list_troubleshoot')
  print(class(BIC_list_troubleshoot))
  print(length(BIC_list_troubleshoot))
  print(dim(BIC_list_troubleshoot[[2]]))
  print(BIC_list_troubleshoot[[2]])

  
  BIC_list_1 = sapply(BIC_list_troubleshoot, function(x){withCallingHandlers({
    min(x[1, x[2,] >= num_edges])
  },
  warning = function(w) {
    if (grepl("no non-missing arguments to min; returning Inf", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  }
  )}
  )
  
  min_ind_1 = which.min(BIC_list_1) # which sparsity level got us the least BIC
  
  x2 <- BIC_list_troubleshoot[[min_ind_1]]
  
  cond <- x2[2, ] >= num_edges
  idx <- which(cond)            
  mid_ind_2 <- idx[which.min(x2[1,][cond])]    # which pinv level got us the least BIC
  
  
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
  
  if(! isSymmetric(Rinv)){
    Rinv = 0.5 * (Rinv + t(Rinv))
  }  
  
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


# print out a lot of statements
# - Error in rho_diag[[j]]$eigenV : $ operator is invalid for atomic vectors

get_cor_gpp_troubleshoot <- function(res=NULL, rho_diag=NULL,NN=NULL,dmax=NULL){
  
  
  p = length(res) ## number of processes
  
  print(class(p))
  print(p)
  print(class(dmax))
  print(dmax)
  
  cov_est = matrix(0, dmax*p, dmax*p)
  
  for(i in 1:p){
    idx_i = 1:dmax + (i-1)*dmax
    for(j in i:p){
      idx_j = 1:dmax + (j-1)*dmax
      if(i==j){
        cov_est[idx_i, idx_j] = rho_diag[[i]]$cov
      }else{
        
        print('printing raw values')
        print(res[[j]])
        print(rho_diag[[j]])
        
        sigma_ij = cross_prod(res[[i]], res[[j]], NN=NN)
        print('got sigma')
        
        tmp = t(rho_diag[[i]]$eigenV) %*% sigma_ij %*% rho_diag[[j]]$eigenV
        print('got tmp')
        cov_est[idx_i, idx_j] = tmp
        print('stored tmp')
        print(paste('success ', as.character(j), sep = as.character(i)))
      }
    }
  }
  print('done with cor_gpp compuation')
  cov_est = cov_est  + t(cov_est) 
  diag(cov_est) = diag(cov_est)/2
  #cor_est = cov2cor(cov_est)
  
  cov_est
}

get_cor_gpp_troubleshoot_2 <- function(i, j, res=NULL, rho_diag=NULL,NN=NULL,dmax=NULL){
  
  

  idx_i = 1:dmax + (i-1)*dmax

  idx_j = 1:dmax + (j-1)*dmax
  
  if(i==j){
    wow <- rho_diag[[i]]$cov
    print('success')
  }else{
    sigma_ij = cross_prod(res[[i]], res[[j]], NN=NN)
    
    print(class(rho_diag[[i]]$eigenV))
    print(class(sigma_ij))
    print(class(rho_diag[[j]]$eigenV))
    
    print(dim(rho_diag[[i]]$eigenV))
    print(dim(sigma_ij))        
    print(dim(rho_diag[[j]]$eigenV))
    
    
    tmp = t(rho_diag[[i]]$eigenV) %*% sigma_ij %*% rho_diag[[j]]$eigenV

    print('success')
  }


  
  return(0)
}

get_rho_diag_pp_reproduce <- function(data_all=NULL,patient_sel=NULL,feature_sel=NULL,
                            Tseq=NULL, dmax=4, gamma_ratio = 1, 
                            gamma_max = 100, remove_diag = TRUE, ncores=3){
  
  # Tmax = 1 ## define domain [0,T]
  # M = 50 ## t_1,...,t_M
  # Tseq = seq(0.1,Tmax-0.1,length=M)
  p = length(feature_sel) ## number of processes
  
  #### estimate intensity ###
  NN = length(patient_sel)
  data_all = data_all[is.element(feature_id, feature_sel),]
  data_all = data_all[is.element(subject_num, patient_sel),]
  
  res = mclapply(feature_sel, function(i){
    
    data_i = data_all[feature_id==i,]
    data_i_count = data_i[,.(count=.N),by="subject_num"]

    
    gamma = get_gamma_reproduce(data_i$time)*gamma_ratio
    gamma[is.na(gamma)] = gamma_max
    # gamma = 1
    # gamma = 100
    denom = sapply(Tseq, function(x){
      truncNorm_denom(x,gamma,0,1)
    })
    
    intensity = data_i[,get_intensity(time,Tseq,gamma,denom),by=c("subject_num")]
    
    patient_now = unique(intensity$subject_num)
    intensity = matrix(intensity$V1,nrow=length(Tseq))
    data_i_count = data_i_count[match(patient_now, subject_num),]
    mu = apply(intensity, 1, sum) / NN
    

    
    diag_mat = data_i[,get_diag_mat(time,Tseq,gamma,denom),by=c("subject_num")]
    diag_mat = matrix(diag_mat$V1,nrow=length(Tseq)*length(Tseq))
    diag_mat = apply(diag_mat, 1, sum)
    diag_mat = matrix(diag_mat, length(Tseq), length(Tseq))
    #diag_mat1 = diag_mat
    
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

get_rho_diag_pp_troubleshoot <- function(data_all=NULL,patient_sel=NULL,feature_sel=NULL,
                            Tseq=NULL, dmax=4, gamma_ratio = 1, 
                            gamma_max = 100, remove_diag = TRUE, ncores=1){
  
  # Tmax = 1 ## define domain [0,T]
  # M = 50 ## t_1,...,t_M
  # Tseq = seq(0.1,Tmax-0.1,length=M)
  p = length(feature_sel) ## number of processes
  
  #### estimate intensity ###
  NN = length(patient_sel)
  data_all = data_all[is.element(feature_id, feature_sel),]
  data_all = data_all[is.element(subject_num, patient_sel),]
  
  print('at res')
  
  # res = mclapply method
  # res2 = lapply
  # res3 = for loop 
  
  res = mclapply(feature_sel, function(i){
    
    data_i = data_all[feature_id==i,]
    data_i_count = data_i[,.(count=.N),by="subject_num"]

    gamma = get_gamma_reproduce(data_i$time)*gamma_ratio
    gamma[is.na(gamma)] = gamma_max

    denom = sapply(Tseq, function(x){
      truncNorm_denom(x,gamma,0,1)
    })
    
    intensity = data_i[,get_intensity(time,Tseq,gamma,denom),by=c("subject_num")]
    
    patient_now = unique(intensity$subject_num)
    intensity = matrix(intensity$V1,nrow=length(Tseq))
    data_i_count = data_i_count[match(patient_now, subject_num),]
    mu = apply(intensity, 1, sum) / NN
    
    # diag_mat = lapply(patient_now, function(id){
    #     time_id = data_i$time[data_i$subject_num==id]
    #     get_diag_mat(time_id, Tseq,gamma,denom)
    # })
    # diag_mat = Reduce("+", diag_mat)
    
    diag_mat = data_i[,get_diag_mat(time,Tseq,gamma,denom),by=c("subject_num")]
    diag_mat = matrix(diag_mat$V1,nrow=length(Tseq)*length(Tseq))
    diag_mat = apply(diag_mat, 1, sum)
    diag_mat = matrix(diag_mat, length(Tseq), length(Tseq))
    #diag_mat1 = diag_mat
    
    list(data_i_count=data_i_count, intensity=intensity, 
         mu = mu, diag_mat=diag_mat, gamma=gamma)
  }, mc.cores = ncores)  
  
  res2 <- lapply(feature_sel, function(i){
    
    
    data_i = data_all[feature_id==i,]
    data_i_count = data_i[,.(count=.N),by="subject_num"]

    
    gamma = get_gamma_reproduce(data_i$time)*gamma_ratio
    gamma[is.na(gamma)] = gamma_max
    # gamma = 1
    # gamma = 100
    denom = sapply(Tseq, function(x){
      truncNorm_denom(x,gamma,0,1)
    })
    
    intensity = data_i[,get_intensity(time,Tseq,gamma,denom),by=c("subject_num")]
    
    patient_now = unique(intensity$subject_num)
    intensity = matrix(intensity$V1,nrow=length(Tseq))
    data_i_count = data_i_count[match(patient_now, subject_num),]
    mu = apply(intensity, 1, sum) / NN
    
    # diag_mat = lapply(patient_now, function(id){
    #     time_id = data_i$time[data_i$subject_num==id]
    #     get_diag_mat(time_id, Tseq,gamma,denom)
    # })
    # diag_mat = Reduce("+", diag_mat)
    
    diag_mat = data_i[,get_diag_mat(time,Tseq,gamma,denom),by=c("subject_num")]
    diag_mat = matrix(diag_mat$V1,nrow=length(Tseq)*length(Tseq))
    diag_mat = apply(diag_mat, 1, sum)
    diag_mat = matrix(diag_mat, length(Tseq), length(Tseq))
    #diag_mat1 = diag_mat
    
    list(data_i_count=data_i_count, intensity=intensity, 
                      mu = mu, diag_mat=diag_mat, gamma=gamma)
  })  
  
  res3 <- list()
  
  for(j in 1:length(feature_sel)){
    
    print(j)
    i <- feature_sel[j]
    
    data_i = data_all[feature_id==i,]
    data_i_count = data_i[,.(count=.N),by="subject_num"]
    
    gamma = get_gamma_reproduce(data_i$time)*gamma_ratio
    gamma[is.na(gamma)] = gamma_max
    # gamma = 1
    # gamma = 100
    denom = sapply(Tseq, function(x){
      truncNorm_denom(x,gamma,0,1)
    })
    
    intensity = data_i[,get_intensity(time,Tseq,gamma,denom),by=c("subject_num")]
    
    patient_now = unique(intensity$subject_num)
    intensity = matrix(intensity$V1,nrow=length(Tseq))
    data_i_count = data_i_count[match(patient_now, subject_num),]
    mu = apply(intensity, 1, sum) / NN
    
    # diag_mat = lapply(patient_now, function(id){
    #     time_id = data_i$time[data_i$subject_num==id]
    #     get_diag_mat(time_id, Tseq,gamma,denom)
    # })
    # diag_mat = Reduce("+", diag_mat)
    
    diag_mat = data_i[,get_diag_mat(time,Tseq,gamma,denom),by=c("subject_num")]
    diag_mat = matrix(diag_mat$V1,nrow=length(Tseq)*length(Tseq))
    diag_mat = apply(diag_mat, 1, sum)
    diag_mat = matrix(diag_mat, length(Tseq), length(Tseq))
    #diag_mat1 = diag_mat
    
    res3[[j]] <- list(data_i_count=data_i_count, intensity=intensity, 
         mu = mu, diag_mat=diag_mat, gamma=gamma)
  }
  
  gamma_vec = sapply(res, function(x){
    x$gamma
  })
  gamma_vec2 = sapply(res2, function(x){
    x$gamma
  })
  gamma_vec3 = sapply(res3, function(x){
    x$gamma
  })
  print('at rho_diag')
  
  # rho_diag = mclapply
  # rho_diag2 <- lapply
  # rho_diag3 = for loop
  
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

  rho_diag2 = lapply(res, function(xx){
    
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
  })    
    
  rho_diag3 <- list()
  
  for(i in 1:length(res)){
    
    print(i)
    
    xx <- res[[i]]

    
    Sigma_ii = cross_prod(X=xx,NN=NN,remove_diag=remove_diag)
    
    ## pca with svd
    eigen_res = eigen(Sigma_ii)
    eigen_res$values[ eigen_res$values<0] = 0 
    positiveInd = eigen_res$values >= 0
    d = eigen_res$values[positiveInd]
    
    print(d)
    
    eigenV = eigen_res$vectors[, positiveInd, drop=FALSE]
    FVE = cumsum(d) / sum(d)
    
    d = d[1:dmax]
    if(any(d < 0)){stop("negative d values!")}
    eigenV = eigenV[, 1:dmax, drop=FALSE]
    
    cov_diag = t(eigenV) %*% Sigma_ii %*% eigenV
    
    rho_diag3[[i]] <- list(d = d, eigenV = eigenV, cumFVE = FVE, 
                          cov = diag(diag(cov_diag),length(d),length(d)),Sigma_ii=Sigma_ii)
  }
  
  list(res       = res,
       res2      = res2,
       res3      = res3,
       rho_diag  = rho_diag,
       rho_diag2 = rho_diag2,
       rho_diag3 = rho_diag3,
       gamma_vec = gamma_vec,
       gamma_vec2 = gamma_vec2,
       gamma_vec3 = gamma_vec3)
}

get_gamma_reproduce <- function(tseq=NULL){
  n = length(tseq)
  if(n>2){
    res = 1/(sum(dist(tseq))*2/n/(n-1))^2
  }else{
    res = NaN
  }
  res
}

cov2cor_manual <- function(Rmat){
  Rmat / sqrt(diag(Rmat) %*% t(diag(Rmat)))
}

adjust_R_troubleshoot <- function(Rmat=NULL){
  
  print('is Rmat symmetric?')
  print(isSymmetric(Rmat))
  print('Rmat dim')
  print(dim(Rmat))
  
  print('checking for NA and Inf')
  print(any(is.na(Rmat)))        # Check for NA
  print(any(is.infinite(Rmat)))  # Check for Inf  
  
  print('Any variances equal to 0?')
  print(any(diag(Rmat) == 0))
  print(which(diag(Rmat) == 0))
  
  Rmat_v2 = cov2cor(Rmat)
  Rmat <- cov2cor_manual(Rmat)
  
  print('are the two cov2cor methods the same?')
  print(table(Rmat_v2 == Rmat))
  
  print(Rmat[45:49, 45:49])
  
  print('done with cov2cor')
  Rmat[Rmat >= 0.99] = 0.99
  Rmat[Rmat <= -0.99] = -0.99
  diag(Rmat) = 1
  
  print('checking for NA and Inf')
  print(any(is.na(Rmat)))        # Check for NA
  print(any(is.infinite(Rmat)))  # Check for Inf  
  
  
  
  print(which(is.na(Rmat), arr.ind = TRUE)) # print NA's
  
  min_eig_val = min(find_min_eigen(c(list(Rmat))))
  if(min_eig_val < pinv_eps){
    min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
    Rmat = adjust_S(Rmat,min_eig_val)
  }
  Rmat
}
