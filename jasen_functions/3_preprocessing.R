# -------------------------------------------------------------------
# 
# Goal: in this part, we have our correlation operator
# 
# recall that FVE_thre = 0.9
#
#
#

adjust_R <- function(Rmat=NULL){
  
  # -------------------------------------------------------------------------
  # 
  # goal: when we truncate Rmat because we don't need all those eigenfunctions,
  #       we need to adjust the shrunken matrix 
  #  
  #
  # input: 
  # - Rmat (80 x 80 matrix), that used to be (120 x 120) but we truncated
  #
  # -------------------------------------------------------------------------
  
  Rmat = cov2cor(Rmat)
  Rmat[Rmat >= 0.99] = 0.99
  Rmat[Rmat <= -0.99] = -0.99
  diag(Rmat) = 1
  
  min_eig_val = min(find_min_eigen(c(list(Rmat)))) # find min eigenvalue
  
  # remember that covariance matrices must have nonnegative eigenvalues
  if(min_eig_val < pinv_eps){
    min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
    Rmat = adjust_S(Rmat,min_eig_val)
  }
  Rmat
}

adjust_S <- function(S=NULL, val=NULL){
  
  # remember that adding c * identity increases all the eigenvalues by c
  diag(S) = diag(S)+val
  cov2cor(S)
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

### choose d based on FVE
FVE_list = lapply(Rmat_diag_full$rho_diag, function(x){x$cumFVE})
d_seq = sapply(FVE_list, function(x){
  which(x>FVE_thre)[1] # which eigenfunction is the first one to cross 90% explanatory power
})

grpind = cumsum(d_seq)
grpind = cbind(c(1, grpind[1:(p-1)]+1),
               grpind[1:p] ) # starting index and ending index


# for each mark, do the following:
# - let d = number of eigenfunctions needed
# - if d exceeds dmax, we only keep dmax eigenfunctions
# - return a list of boolean vectors, 
col_keep =  lapply(FVE_list, function(x){
  d = which(x>FVE_thre)[1]
  col_ind = rep(FALSE,dmax)
  col_ind[1:d] = TRUE
  col_ind
})

# vectorize
col_keep = do.call(c, col_keep)

# only keep the important rows and columns from Rmat (120 x 120)
Rmat = Rmat[col_keep,col_keep]

### deal with small/negative eigenvalues of the corr matrix 
Rmat_IC = adjust_R(Rmat)

###
graph_all = list()
res_all = list()