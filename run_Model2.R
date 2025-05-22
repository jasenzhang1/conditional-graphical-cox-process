
### Model 2 ###
rm(list=ls())

ffs = 1
set.seed(ffs)
source("GraphPP_FUN.R")

## simulation setting

nrep = 100
#p = 100
ncores = 2
#ntrain = 500

pseq = 40
nseq = c(400,500,600)

np_seq = cbind(rep(pseq, each=length(nseq)), 
               rep(nseq, length(pseq)))

simu_now  = "simu2"
rep_num = ffs %% nrep
if(rep_num==0){rep_num=100}

np_ind = ceiling(ffs/nrep)
p = np_seq[np_ind,1]
ntrain = np_seq[np_ind,2]

simu_file = paste(simu_now, ".R",sep="")

#### simulate the data ####

source(file = simu_file)

### get corr matrix for full data
patient_sel = 1:ntrain
feature_sel = 1:p
Tseq = seq(0.05,0.95,length=19)
dmax = 3
FVE_thre = 0.9

quantile(data_df[,(.N),by=c("subject_num","feature_id")]$V1)

Rmat_diag_full = get_rho_diag_pp(data_all=data_df,patient_sel=patient_sel,
                                 feature_sel=feature_sel,Tseq=Tseq, dmax=dmax,
                                  ncores=ncores)

Rmat = get_cor_gpp(res=Rmat_diag_full$res, rho_diag=Rmat_diag_full$rho_diag,
                   NN=length(patient_sel),dmax=dmax)


### choose d based on FVE
FVE_list = lapply(Rmat_diag_full$rho_diag, function(x){x$cumFVE})
d_seq = sapply(FVE_list, function(x){
    which(x>FVE_thre)[1]
})
grpind = cumsum(d_seq)
grpind = cbind(c(1, grpind[1:(p-1)]+1),
               grpind[1:p] )
col_keep =  lapply(FVE_list, function(x){
    d = which(x>FVE_thre)[1]
    col_ind = rep(FALSE,dmax)
    col_ind[1:d] = TRUE
    col_ind
})
col_keep = do.call(c, col_keep)

Rmat = Rmat[col_keep,col_keep]

### deal with small/negative eigenvalues of the corr matrix 
Rmat_IC = adjust_R(Rmat)

###
graph_all = list()
res_all = list()

### get sparse corr matrix using GPP method ####

# BIC 
factor = sqrt(ntrain)
res_GPP_BIC = run_GPP_HT_BIC(Rmat_IC, grpind, ntrain,factor)
graph_all[["GPP_BIC"]] = as.matrix((get_groupNorm(res_GPP_BIC, grpind)!=0)+0)
res_all[["GPP_BIC"]] = res_GPP_BIC


#### get sparse corr matrix using count method ###

S_count = get_cor_count(data_all=data_df,patient_sel=patient_sel,
                        feature_sel=feature_sel)

grpind = cbind(seq_along(feature_sel),seq_along(feature_sel))

# BIC and AIC
factor = sqrt(ntrain)
res_count_BIC = run_GPP_HT_BIC(S_count, grpind, ntrain, factor)
graph_all[["count_BIC"]] = as.matrix((get_groupNorm(res_count_BIC, grpind)!=0)+0)
res_all[["count_BIC"]] = res_count_BIC


### save results 
accu_all = lapply(graph_all, function(tmp){
    get_accuracy(EstG = tmp, TrueG = graph_true)
})



# filename = paste("result_simu/",simu_now,"_BIC_n=",ntrain,"_p=",p,"_",rep_num ,".rda", sep="" )
# 
# save(accu_all, graph_true, graph_all, res_all, file=filename)
# 
# quit(save="no")


