
# i want this to return a graph with at least q edges with the minimum GIC

run_GPP_HT_BIC_v2 <- function(Rmat=NULL, grpind=NULL,ntrain=NULL, factor=NULL, num_edges=NULL){
  
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
      
      # modify this
      BIC_i <- calc_BIC_1(x, Rmat, ntrain, grpind, factor)
      graph_i <- as.matrix((get_groupNorm(x, grpind)!=0)+0)
      num_edges_i <- sum( graph_i[upper.tri(graph_i)] )
        
      c(BIC_i, num_edges_i)
    })
    BIC_val
  })
  
  # filter the BIC's with at least 10 edges and then remove # edges 
  

  BIC_list_1 = sapply(BIC_list, function(x){ min(x[1, x[2,] >= num_edges]) })
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
