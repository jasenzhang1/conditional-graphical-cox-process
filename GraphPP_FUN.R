## Graph-PP fun ##
require(Matrix)
require(Rcpp)
require(methods)
require(parallel)
require(data.table)
require(pracma)
require(PtProcess)
require(MASS)

pinv_eps = 1e-2
#### get bandwidth estimate  ####
get_gamma <- function(tseq=NULL){
    n = length(tseq)
    if(n>2& n<10000){
        res = 1/(sum(dist(tseq))*2/n/(n-1))^2
    }else if(n >= 10000){
        tseq = sample(tseq, 10000)
        n = length(tseq)
        res = 1/(sum(dist(tseq))*2/n/(n-1))^2
    }else{
        res = NaN
    }
    res
}

# get_gamma <- function(tseq=NULL){
#     n = length(tseq)
#     if(n>2){
#         res = 1/(sum(dist(tseq))*2/n/(n-1))^2
#     }else{
#         res = NaN
#     }
#     res
# }


#### truncated normal denominator for [0,T] ####
truncNorm_denom <- function(t=NULL, gamma=1 ,a=0, b=1){
    integrate(function(x) {exp(-gamma*(t-x)^2)}, a, b)$value
}

#### function to estimate intensity
get_intensity <- function(xx=NULL, Tseq=NULL, gamma=NULL,denom=NULL){
    diff = outer(Tseq,xx,'-')
    diff = exp(-gamma* (diff^2))
    res = apply(diff,1,sum)/denom
    names(res) = paste("V",1:length(Tseq),sep="")
    #data.table(t(res))
    res
}

### function to diagomal elements matrix
get_diag_mat <- function(xx=NULL, Tseq=NULL, gamma=NULL,denom=NULL){
    diff = outer(Tseq,xx,'-')
    diff = exp(-gamma* (diff^2))/denom
    diff%*%t(diff)
}

#### functions to get correlation matrix ####

cross_prod <- function(X=NULL, Y=NULL, NN=NULL,remove_diag = TRUE){
    
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


get_rho_diag_pp <- function(data_all=NULL,patient_sel=NULL,feature_sel=NULL,
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
        # gamma_tmp = data_i[,get_gamma(time),
        #                    by=c("subject_num")]
        # gamma_tmp = gamma_tmp[!is.na(V1),]
        # gamma = median(gamma_tmp$V1)
        # gamma[is.na(gamma)] = gamma_max
        # gamma[gamma > gamma_max] = gamma_max
        
        gamma = get_gamma(data_i$time)*gamma_ratio
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

get_rho_off_diag_pp <- function(data_all=NULL,patient_sel=NULL,feature_sel=NULL,
                                Tseq=NULL, dmax=4,  gamma_vec = NULL, ncores=3){
    
    # Tmax = 1 ## define domain [0,T]
    # M = 50 ## t_1,...,t_M
    # Tseq = seq(0.1,Tmax-0.1,length=M)
    p = length(feature_sel) ## number of processes
    
    #### estimate intensity ###
    NN = length(patient_sel)
    data_all = data_all[is.element(feature_id, feature_sel),]
    data_all = data_all[is.element(subject_num, patient_sel),]
    
    res = mclapply(seq_along(feature_sel), function(i){
        
        fid = feature_sel[i]
        data_i = data_all[feature_id==fid,]
        data_i_count = data_i[,.(count=.N),by="subject_num"]
        gamma = gamma_vec[i]
        
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
        
        list(data_i_count=data_i_count, intensity=intensity, 
             mu = mu, diag_mat=diag_mat, gamma=gamma)
    }, mc.cores = ncores)
    
    list(res=res)
}

get_cor_gpp <- function(res=NULL, rho_diag=NULL,NN=NULL,dmax=NULL){
    
    
    p = length(res) ## number of processes
    cov_est = matrix(0, dmax*p, dmax*p)
    
    for(i in 1:p){
        idx_i = 1:dmax + (i-1)*dmax
        for(j in i:p){
            idx_j = 1:dmax + (j-1)*dmax
            if(i==j){
                cov_est[idx_i, idx_j] = rho_diag[[i]]$cov
            }else{
                sigma_ij = cross_prod(res[[i]], res[[j]], NN=NN)
                tmp = t(rho_diag[[i]]$eigenV) %*% sigma_ij %*% rho_diag[[j]]$eigenV
                cov_est[idx_i, idx_j] = tmp
            }
        }
    }
    cov_est = cov_est  + t(cov_est) 
    diag(cov_est) = diag(cov_est)/2
    #cor_est = cov2cor(cov_est)
    
    cov_est
}


get_cor_gpp_sub <- function(iseq=NULL, res=NULL, rho_diag=NULL,NN=NULL,dmax=NULL){
    
    
    p = length(res) ## number of processes
    cov_est = matrix(0, dmax*p, dmax*p)
    
    for(i in iseq){
        idx_i = 1:dmax + (i-1)*dmax
        for(j in i:p){
            idx_j = 1:dmax + (j-1)*dmax
            if(i==j){
                cov_est[idx_i, idx_j] = rho_diag[[i]]$cov
            }else{
                sigma_ij = cross_prod(res[[i]], res[[j]], NN=NN)
                tmp = t(rho_diag[[i]]$eigenV) %*% sigma_ij %*% rho_diag[[j]]$eigenV
                cov_est[idx_i, idx_j] = tmp
            }
        }
    }
    
    #    cov_est = cov_est  + t(cov_est) 
    #    diag(cov_est) = diag(cov_est)/2
    #    cor_est = cov2cor(cov_est)
    
    cov_est
}

create_rho_mat <- function(data_all=NULL,patient_sel=NULL,feature_sel=NULL,
                           Tseq=NULL, dmax=4,  
                           gamma_max = 100, remove_diag = TRUE, ncores=3){
    
    Rmat_diag = get_rho_diag_pp(data_all=data_all,patient_sel=patient_sel,
                                feature_sel=feature_sel,
                                Tseq=Tseq, dmax=dmax,  
                                gamma_max = gamma_max, 
                                remove_diag = remove_diag, ncores=ncores)
    Rmat = get_cor_gpp(res=Rmat_diag$res, rho_diag=Rmat_diag$rho_diag,
                       NN=length(patient_sel),dmax=dmax)
    Rmat
}


create_rho_mat_V2 <- function(data_all=NULL,patient_sel=NULL,feature_sel=NULL,
                           Tseq=NULL, dmax=4,   gamma_max = 100, 
                           nfold=0, remove_diag = TRUE, ncores=3){
    
    Rmat_diag = get_rho_diag_pp(data_all=data_all,patient_sel=patient_sel,
                                feature_sel=feature_sel,
                                Tseq=Tseq, dmax=dmax,  
                                gamma_max = gamma_max, 
                                remove_diag = remove_diag, ncores=ncores)
    Rmat = get_cor_gpp(res=Rmat_diag$res, rho_diag=Rmat_diag$rho_diag,
                       NN=length(patient_sel),dmax=dmax)
    
    Rmat_train = list()
    Rmat_test = list()
    
    if(nfold>0){
        test_list = split(patient_sel, rep(1:nfold,length=length(patient_sel)))
        train_list = lapply(seq_along(test_list), function(i){
            setdiff(patient_sel, test_list[[i]])
        })
        
        
        Rmat_train = lapply(train_list, function(id){
            Rmat_diag_tmp = get_rho_diag_pp(data_all=data_all,patient_sel=id,
                                        feature_sel=feature_sel,
                                        Tseq=Tseq, dmax=dmax,  
                                        gamma_max = gamma_max, 
                                        remove_diag = remove_diag, ncores=ncores)
            get_cor_gpp(res=Rmat_diag_tmp$res, rho_diag=Rmat_diag$rho_diag,
                        NN=length(id),dmax=dmax)
        })
        
        Rmat_test = lapply(test_list, function(id){
            Rmat_diag_tmp = get_rho_diag_pp(data_all=data_all,patient_sel=id,
                                            feature_sel=feature_sel,
                                            Tseq=Tseq, dmax=dmax,  
                                            gamma_max = gamma_max, 
                                            remove_diag = remove_diag, ncores=ncores)
            get_cor_gpp(res=Rmat_diag_tmp$res, rho_diag=Rmat_diag$rho_diag,
                        NN=length(id),dmax=dmax)
        })
        
    }
    
    list(Rmat=Rmat, Rmat_diag=Rmat_diag,
         Rmat_train=Rmat_train, Rmat_test=Rmat_test)
}

find_min_eigen <- function(S=NULL){
    eig = sapply(S, function(s){
        min(eigen(s)$val)
    })
    eig
}

adjust_S <- function(S=NULL, val=NULL){
    diag(S) = diag(S)+val
    cov2cor(S)
}

adjust_R <- function(Rmat=NULL){
    Rmat = cov2cor(Rmat)
    Rmat[Rmat >= 0.99] = 0.99
    Rmat[Rmat <= -0.99] = -0.99
    diag(Rmat) = 1
    
    min_eig_val = min(find_min_eigen(c(list(Rmat))))
    if(min_eig_val < pinv_eps){
        min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
        Rmat = adjust_S(Rmat,min_eig_val)
    }
    Rmat
}

#### functions to get sparse correlation matrix ####

get_groupNorm <- function(S=NULL, grpind=NULL){
    
    res = matrix(0, nrow(grpind), nrow(grpind))
    for(i in 1:nrow(grpind) ){
        for(j in i:nrow(grpind) ){
            tmp = S[grpind[i,1]:grpind[i,2],
                    grpind[j,1]:grpind[j,2], drop=FALSE]
            res[i,j] = sum(tmp^2)/ length(tmp)
        }
    }
    
    res = res + t(res)
    diag(res) = diag(res)/2
    res = sqrt(res)
    res
}


Cov_hardT <- function(S=NULL, grpind=NULL, lamseq=NULL){
    
    Grp_norm = get_groupNorm(S, grpind)
    omega_arr = array(0, c(nrow(S), ncol(S), length(lamseq)))
    
    for(i in seq_along(lamseq)){
        lam1 = lamseq[i]
        X = S
        
        for(j in 1:(nrow(Grp_norm)-1)){
            for( k in (j+1):nrow(Grp_norm)){
                if(Grp_norm[j,k] <= lam1){
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


calc_error_loglik <- function(omega=NULL, S=NULL){
    val = determinant(omega)$modulus
    sum(diag((omega%*%S))) - val
}

calc_error_F <- function(R=NULL, S=NULL){
    sum((R-S)^2)
}


get_lamseq <- function(Rmat=NULL, grpind=NULL, len=100){
    
    Rmat_norm = get_groupNorm(Rmat,grpind)
    lam1seq  = sort(unique(Rmat_norm[upper.tri(Rmat_norm)]), decreasing = T)
    
    if(length(lam1seq)>100){
        lam1seq = quantile(lam1seq,seq(1,0,length=100))
    }
    lam1seq
}

tune_HT_1st <- function(Rmat=NULL, Rmat_train_list=NULL, Rmat_test_list=NULL,
                        grpind=NULL){
    
    lamseq = get_lamseq(Rmat, grpind)
    cv_err = NULL
    for(i in seq_along(Rmat_train_list)){
        Rmat_train_thre_1 = Cov_hardT(S=Rmat_train_list[[i]], grpind=grpind, lamseq=lamseq)
        
        cv_err_tmp = sapply(seq_along(lamseq), function(ii){
            sum((Rmat_train_thre_1[,,ii] - Rmat_test_list[[i]])^2)
        })
        cv_err = rbind(cv_err, cv_err_tmp)
    }
    
    cv_se = apply(cv_err, 2, sd)/ sqrt(length(Rmat_train_list))
    cv_err = apply(cv_err, 2, mean)
    
    
    min_ind = which.min(cv_err)
    lam_min = lamseq[min_ind]
    lam_1se = lamseq[ which(cv_err <= (cv_err+cv_se)[min_ind])[1] ]
    
    Rmat_min = Cov_hardT(S=Rmat, grpind=grpind, lamseq=lam_min)
    Rmat_min = Rmat_min[,,1]
    Rmat_1se = Cov_hardT(S=Rmat, grpind=grpind, lamseq=lam_1se)
    Rmat_1se = Rmat_1se[,,1]
    
    list(cv_err=cv_err, cv_se=cv_se,lamseq=lamseq,
         lam_min=lam_min, lam_1se=lam_1se,
         Rmat_min=Rmat_min, Rmat_1se=Rmat_1se)
}

tune_HT_2nd <- function(Rinv=NULL, Rmat_train_list=NULL, Rmat_test_list=NULL,
                        grpind=NULL, lam1=NULL){
    
    lamseq = get_lamseq(Rinv, grpind)
    cv_err = NULL
    for(i in seq_along(Rmat_train_list)){
        
        Rinv_train = Cov_hardT(S=Rmat_train_list[[i]], grpind=grpind, lamseq=lam1)
        Rinv_train = Rinv_train[,,1]
        #Rinv_train = adjust_S(Rinv_train)
        #Rinv_train = solve(Rinv_train)
        Rinv_train = pinv(Rinv_train,pinv_eps)
        
        Rmat_train_thre_2 = Cov_hardT(S=Rinv_train, grpind=grpind, lamseq=lamseq)
        
        cv_err_tmp = sapply(seq_along(lamseq), function(ii){
            omega = Rmat_train_thre_2[,,ii]
            S = Rmat_test_list[[i]]
            val = determinant(omega)$modulus
            sum(diag((omega%*%S))) - val
        })
        cv_err = rbind(cv_err, cv_err_tmp)
    }
    
    #cv_err = apply(cv_err, 2, mean)
    
    # lam2 = lam2seq[which.min(cv_err)]
    # 
    # 
    # Rinv = Cov_hardT(S=Rinv, grpind=grpind, lam1seq=lam2)
    # Rinv = Rinv[,,1]
    # 
    cv_se = apply(cv_err, 2, sd)/ sqrt(length(Rmat_train_list))
    cv_err = apply(cv_err, 2, mean)
    
    
    min_ind = which.min(cv_err)
    lam_min = lamseq[min_ind]
    lam_1se = lamseq[ which(cv_err <= (cv_err+cv_se)[min_ind])[1] ]
    
    Rinv_min = Cov_hardT(S=Rinv, grpind=grpind, lamseq=lam_min)
    Rinv_min = Rinv_min[,,1]
    Rinv_1se = Cov_hardT(S=Rinv, grpind=grpind, lamseq=lam_1se)
    Rinv_1se = Rinv_1se[,,1]
    
    list(cv_err=cv_err, cv_se=cv_se,lamseq=lamseq,
         lam_min=lam_min, lam_1se=lam_1se,
         Rinv_min=Rinv_min, Rinv_1se=Rinv_1se)
}

run_GPP_HT <- function(Rmat=NULL, Rmat_train_list=NULL, 
                       Rmat_test_list=NULL, grpind=NULL,
                       type = c("min","1se","mixed")[1]) {
    
    # 1st hard-thresholding 
    res_tune_1st = tune_HT_1st(Rmat=Rmat, Rmat_train_list=Rmat_train_list,
                               Rmat_test_list=Rmat_test_list, grpind=grpind)
    
    if(type=="min" | type=="mixed"){
        Rmat_1st = res_tune_1st$Rmat_min
        Rinv = pinv(Rmat_1st,pinv_eps)
        lam1 = res_tune_1st$lam_min
    }else{
        Rmat_1st = res_tune_1st$Rmat_1se
        Rinv = pinv(Rmat_1st,pinv_eps)
        lam1 = res_tune_1st$lam_1se
    }
   
    # 2nd hard-thresholding 
    res_tune_2nd = tune_HT_2nd(Rinv=Rinv, Rmat_train_list=Rmat_train_list,
                               Rmat_test_list=Rmat_test_list,
                               grpind=grpind, lam1 = lam1)
    if(type=="min"){
        Rinv = res_tune_2nd$Rinv_min
    }else{
        Rinv = res_tune_2nd$Rinv_1se
    }
    
    Rinv
}




tune_HT_2nd_2 <- function(Rinv=NULL, Rmat_train_list=NULL, Rmat_test_list=NULL,
                        grpind=NULL, lam1=NULL){
    
    lamseq = get_lamseq(Rinv, grpind)
    cv_err = NULL
    for(i in seq_along(Rmat_train_list)){
        
        Rinv_train = Cov_hardT(S=Rmat_train_list[[i]], grpind=grpind, lamseq=lam1)
        Rinv_train = Rinv_train[,,1]
        #Rinv_train = adjust_S(Rinv_train)
        #Rinv_train = solve(Rinv_train)
        Rinv_train = pinv(Rinv_train,pinv_eps)
        
        Rmat_train_thre_2 = Cov_hardT(S=Rinv_train, grpind=grpind, lamseq=lamseq)
        
        cv_err_tmp = sapply(seq_along(lamseq), function(ii){
            omega = Rmat_train_thre_2[,,ii]
            S = Rmat_test_list[[i]]
            val = determinant(omega)$modulus
            #val = eigen(omega)$value
            #val = val[val>pinv_eps]
            #val = sum(log(val))
            sum(diag((omega%*%S))) - val
        })
        cv_err = rbind(cv_err, cv_err_tmp)
    }
    
    #cv_err = apply(cv_err, 2, mean)
    
    # lam2 = lam2seq[which.min(cv_err)]
    # 
    # 
    # Rinv = Cov_hardT(S=Rinv, grpind=grpind, lam1seq=lam2)
    # Rinv = Rinv[,,1]
    # 
    cv_se = apply(cv_err, 2, sd)/ sqrt(length(Rmat_train_list))
    cv_err = apply(cv_err, 2, mean)
    
    
    min_ind = which.min(cv_err)
    lam_min = lamseq[min_ind]
    lam_1se = lamseq[ which(cv_err < (cv_err+cv_se)[min_ind])[1] ]
    
    Rinv_min = Cov_hardT(S=Rinv, grpind=grpind, lamseq=lam_min)
    Rinv_min = Rinv_min[,,1]
    Rinv_1se = Cov_hardT(S=Rinv, grpind=grpind, lamseq=lam_1se)
    Rinv_1se = Rinv_1se[,,1]
    
    list(cv_err=cv_err, cv_se=cv_se,lamseq=lamseq,
         lam_min=lam_min, lam_1se=lam_1se,
         Rinv_min=Rinv_min, Rinv_1se=Rinv_1se)
}

run_GPP_HT_2 <- function(Rmat=NULL, Rmat_train_list=NULL, 
                       Rmat_test_list=NULL, grpind=NULL,
                       type = c("min","1se")[1]) {
    
    # 1st hard-thresholding 
    res_tune_1st = tune_HT_1st(Rmat=Rmat, Rmat_train_list=Rmat_train_list,
                               Rmat_test_list=Rmat_test_list, grpind=grpind)
    
    if(type=="min"){
        Rmat_1st = res_tune_1st$Rmat_min
        Rinv = pinv(Rmat_1st,pinv_eps)
        lam1 = res_tune_1st$lam_min
    }else{
        Rmat_1st = res_tune_1st$Rmat_1se
        Rinv = pinv(Rmat_1st,pinv_eps)
        lam1 = res_tune_1st$lam_1se
    }
    
    # 2nd hard-thresholding 
    res_tune_2nd = tune_HT_2nd_2(Rinv=Rinv, Rmat_train_list=Rmat_train_list,
                               Rmat_test_list=Rmat_test_list,
                               grpind=grpind, lam1 = lam1)
    if(type=="min"){
        Rinv = res_tune_2nd$Rinv_min
    }else{
        Rinv = res_tune_2nd$Rinv_1se
    }
    
    Rinv
}


get_cor_count <- function(data_all=NULL, patient_sel=NULL, feature_sel=NULL){
    
    N = length(patient_sel)
    p = length(feature_sel)
    
    data_all = data_all[is.element(feature_id, feature_sel),]
    data_all = data_all[is.element(subject_num, patient_sel),]
    
    setorder(data_all, subject_num, feature_id, time)
    
    data.cov = data_all[,.(count = .N),by=c("subject_num","feature_id")]
    data.cov[, i:= match(subject_num, patient_sel) ]
    data.cov[, j:= match(feature_id , feature_sel)]
    
    mat_dx = sparseMatrix(i=data.cov[,i], j=data.cov[,j], x=data.cov[,count],
                          dims = c(N,p))
    mat_dx@x = log(mat_dx@x+1)
    
    cor(as.matrix(mat_dx))
}

get_accuracy <- function(EstG=NULL, TrueG = NULL){
    est = EstG[upper.tri(EstG)]
    true = TrueG[upper.tri(TrueG)]
    accu = sum(est==true)/length(est)
    TP = sum(est==1 & true==1)
    FP = sum(est==1 & true==0)
    FN = sum(est==0 & true==1)
    TN = sum(est==0 & true==0) 
    TPR = TP / sum(true==1)
    TNR = TN / sum(true==0)
    FPR = FP / sum(true==0)
    FNR = FN / sum(true==1)
    PPV = TP / sum(est==1)
    NPV = TN / sum(est==0)
    F1 = 2*TP/(2*TP+FP+FN)
    res = c(accu,TPR,TNR,FPR,FNR, PPV, NPV,F1,TP,FP,FN,TN)
    names(res) = c('Accuracy','TPR','TNR', 'FPR', 'FNR','PPV','NPV','F1',
                   'TP','FP','FN','TN')
    res
}



get_nonzero <- function(S=NULL, grpind=NULL){
    
    for(i in 1:nrow(grpind) ){
        S[grpind[i,1]:grpind[i,2],grpind[i,1]:grpind[i,2]] = 0
    }
    sum(abs(S[lower.tri(S)]) !=0)
}

calc_BIC_1 <- function(Rinv=NULL, Rmat=NULL, n=NULL,grpind=NULL,factor=NULL){
    #val = determinant(Rinv)$modulus
    ss = eigen(Rinv)$values
    val = sum(log(ss[ss>pinv_eps]))
    loglik = n*(sum(diag((Rinv%*%Rmat))) - val)
    df = get_nonzero(Rinv, grpind)
    loglik + factor*df
}


calc_AIC_1 <- function(Rinv=NULL, Rmat=NULL, n=NULL,grpind=NULL){
    #val = determinant(Rinv)$modulus
    ss = eigen(Rinv)$values
    val = sum(log(ss[ss>pinv_eps]))
    loglik = n*(sum(diag((Rinv%*%Rmat))) - val)
    df = get_nonzero(Rinv, grpind)
    loglik + 2*df
}



run_GPP_HT_BIC <- function(Rmat=NULL, grpind=NULL,ntrain=NULL, factor=NULL){
    
    lamseq = get_lamseq(Rmat, grpind)
    Rmat_thre_1 = Cov_hardT(S=Rmat, grpind=grpind, lamseq=lamseq)
    
    BIC_list = lapply(1:dim(Rmat_thre_1)[3], function(i){
        Rmat_1st  = Rmat_thre_1[,,i]
        min_eig_val = min(find_min_eigen(c(list(Rmat_1st))))
        if(min_eig_val < pinv_eps){
            min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
            Rmat_1st = adjust_S(Rmat_1st,min_eig_val)
        }
        
        Rinv = pinv(Rmat_1st,pinv_eps)
        lamseq_2 = get_lamseq(Rinv, grpind)
        Rinv_thre_2 = Cov_hardT(S=Rinv, grpind=grpind, lamseq=lamseq_2)
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


run_GPP_HT_AIC <- function(Rmat=NULL, grpind=NULL, ntrain=NULL){
    
    lamseq = get_lamseq(Rmat, grpind)
    Rmat_thre_1 = Cov_hardT(S=Rmat, grpind=grpind, lamseq=lamseq)
    
    AIC_list = lapply(1:dim(Rmat_thre_1)[3], function(i){
        Rmat_1st  = Rmat_thre_1[,,i]
        min_eig_val = min(find_min_eigen(c(list(Rmat_1st))))
        if(min_eig_val < pinv_eps){
            min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
            Rmat_1st = adjust_S(Rmat_1st,min_eig_val)
        }
        
        Rinv = pinv(Rmat_1st,pinv_eps)
        lamseq_2 = get_lamseq(Rinv, grpind)
        Rinv_thre_2 = Cov_hardT(S=Rinv, grpind=grpind, lamseq=lamseq_2)
        AIC_val = apply(Rinv_thre_2, 3, function(x){
            #calc_BIC_test(x, Rmat_1st, ntrain, grpind)
            
            min_eig_val = min(find_min_eigen(c(list(x))))
            if(min_eig_val < pinv_eps){
                min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
                diag(x) = diag(x) + min_eig_val
            }
            
            calc_AIC_1(x, Rmat, ntrain, grpind)
        })
        AIC_val
    })
    
    AIC_list_1 = sapply(AIC_list, function(x){min(x)})
    min_ind_1 = which.min(AIC_list_1)
    mid_ind_2 = which.min(AIC_list[[min_ind_1]])
    
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

