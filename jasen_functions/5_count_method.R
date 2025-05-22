#
#
# Use the default count method
#
#
#

data_all <- data_df

get_cor_count <- function(data_all=NULL, patient_sel=NULL, feature_sel=NULL){
  
  # --------------------------------------------------------------------------
  #
  # goal:
  # - 
  #
  #
  # input: 
  # - data_all      (data.frame):    3 column dataframe of the raw data
  #                                    - subject_num: which replicate? 1:n
  #                                    - feature_id:  which mark? 1:p
  #                                    - time:        which timestamp? any time in [0, TT]
  #
  #
  # - patient_sel       (vector):    1:ntrain
  # - feature_sel       (vector):    1:p
  #
  # output:
  # - a (p x p) matrix that seems to be taking the norm of each block matrix from Rmat
  # - it seems like we are taking the rmse: sqrt(sum(elements^2) / num(elements))
  #
  # --------------------------------------------------------------------------  
  
  N = length(patient_sel)
  p = length(feature_sel)
  
  # 1) ensure the dataset is proper + sorted
  data_all = data_all[is.element(feature_id, feature_sel),]
  data_all = data_all[is.element(subject_num, patient_sel),]
  
  setorder(data_all, subject_num, feature_id, time)
  
  # 2) create a table of counts
  data.cov = data_all[,.(count = .N),by=c("subject_num","feature_id")]
  
  # 2a) for each subject_num and feature_id, find its respective index in patient_sel or feature_sel
  #     but if patient_sel and feature_sel are just 1:n, 1:p, this is trivially just subject_num and feature_id again
  data.cov[, i:= match(subject_num, patient_sel) ]
  data.cov[, j:= match(feature_id , feature_sel)]
  
  
  # 3) create a sparse matrix object (400 x 40)
  mat_dx = sparseMatrix(i=data.cov[,i], j=data.cov[,j], x=data.cov[,count],
                        dims = c(N,p))
  
  # 3a) log the counts
  mat_dx@x = log(mat_dx@x+1)
  
  # 4) return the correlation matrix (40 x 40)
  cor(as.matrix(mat_dx))
  
  
  
}


# 40 x 40 sample correlation
S_count = get_cor_count(data_all=data_df,patient_sel=patient_sel,
                        feature_sel=feature_sel)


# obtain the trivial indexing scheme. Each class only has width 1 
grpind = cbind(seq_along(feature_sel),seq_along(feature_sel))

# BIC and AIC
factor = sqrt(ntrain)
res_count_BIC = run_GPP_HT_BIC(S_count, grpind, ntrain, factor)
graph_all[["count_BIC"]] = as.matrix((get_groupNorm(res_count_BIC, grpind)!=0)+0)
res_all[["count_BIC"]] = res_count_BIC
