
# qrsh

ffs = 1
set.seed(ffs)
source("GraphPP_FUN.R")
source("modified_funcs.R")
library(dplyr)



#### simulate the data ####

# 1) Open 99_Create_Spiketrain_Dataset and load the data_df file

IDs <- c('346', '351', '366', '361', '362', '368')  # mouse ID
ID2 <- c('Tau1', 'Tau2', 'Tau3', 'WT1', 'WT2', 'WT3') # our name 

IDs <- c('346', '361')  # mouse ID
ID2 <- c('Tau1', 'WT1') # our name 

n <- 50    #number of replicates
movements <- c(0, 1, 2)
movements <- c(1)
weeks <- c(18, 22, 26)
VR <- 0
min_edges <- 1
ncores <- parallel::detectCores() - 1

print('Number of Cores')
print(ncores)
print(Sys.time()) 

for(week in weeks){ # for each week
    
    for(movement in movements){ # for each movement
        
        for(i in 1:length(IDs)){ # for each mouse
            
            print(paste(as.character(week), as.character(i), sep = as.character(movement)))
            
            data_root <- "spike_data/"
            
            setting_ID <- paste('w', '_m', sep = as.character(week))                # w17_m
            setting_ID <- paste(setting_ID, 'vr', sep = as.character(movement))     # w17_m2vr
            setting_ID <- paste(setting_ID, '_n', sep = as.character(VR))           # w17_m2vr0_n
            setting_ID <- paste(setting_ID, as.character(n), sep = '')              # w17_m2vr0_n400
            setting_ID2 <- paste(setting_ID, as.character(min_edges), sep = '_me')  # w17_m2vr0_n400_me1
            
            file_ID  <- paste(ID2[i], setting_ID, sep = '_')                        # Tau1_w17_m2vr0_n400
            file_ID2 <- paste(ID2[i], setting_ID2, sep = '_')                       # Tau1_w17_m2vr0_n400_me1
            
            data_root <- paste(data_root, '/', sep  = setting_ID)                   # ../spike_data/w17_m2vr0_n400/
            
            data_ID <- paste(data_root, '.RData', sep = file_ID)                    # ../spike_data/w17_m2vr0_n400/Tau1_e1_m2vr0_n400.RData
            
            load(data_ID) # data_df2
            
            # which neurons fire way too little and must be discarded (under 50 total firings = discard)
            
            data_df2 <- as.data.table(data_df2)
            
            data_df3 <- data_df2[, if (.N >= 50) .SD, by = feature_id]
            
            print(length(unique(data_df2$feature_id)))
            print(length(unique(data_df3$feature_id)))
            
            missing_neurons <- setdiff(unique(data_df2$feature_id),
                                       unique(data_df3$feature_id))
            
            # retroactively get parameters
            
            ntrain <- length(unique(data_df3$subject_num))
            p <- length(unique(data_df3$feature_id))
            
            
            ### get corr matrix for full data
            patient_sel = 1:ntrain
            feature_sel = unique(data_df3$feature_id) %>% sort()
            
            Tseq = seq(0.05,0.95,length=19)
            dmax = 3
            FVE_thre = 0.9 # originally 0.9
            
            quantile(data_df3[,(.N),by=c("subject_num","feature_id")]$V1)
            
            # Rmat_diag_full_ts = get_rho_diag_pp_troubleshoot(data_all=data_df3,
            #                                  patient_sel=patient_sel,
            #                                  feature_sel=feature_sel,
            #                                  Tseq=Tseq,
            #                                  dmax=dmax,
            #                                  ncores=ncores)
            
            Rmat_diag_full = get_rho_diag_pp_reproduce(data_all=data_df3,
                                                       patient_sel=patient_sel,
                                                       feature_sel=feature_sel,
                                                       Tseq=Tseq,
                                                       dmax=dmax,
                                                       ncores=ncores)
            
            print('checkpoint 1')
            
            # Rmat_ts = get_cor_gpp_troubleshoot(res=Rmat_diag_full_ts$res, rho_diag=Rmat_diag_full_ts$rho_diag,
            #                                 NN=length(patient_sel),dmax=dmax)

            Rmat = get_cor_gpp(res=Rmat_diag_full$res, rho_diag=Rmat_diag_full$rho_diag,
                               NN=length(patient_sel),dmax=dmax)
            
            print(min(diag(Rmat)))
                        
            print('checkpoint 2')
            
            ### choose d based on FVE
            FVE_list = lapply(Rmat_diag_full$rho_diag, function(x){x$cumFVE})
            

            

            
            d_seq = sapply(FVE_list, function(x){
                which(x>FVE_thre)[1]
            })
            
            if(max(d_seq) > dmax){
                print('needed to trigger dseq truncation!')
                d_seq <- pmin(d_seq, dmax)
            }

            
            grpind = cumsum(d_seq)
            

            
            grpind = cbind(c(1, grpind[1:(p-1)]+1),
                           grpind[1:p] )
            

            
            col_keep =  lapply(FVE_list, function(x){
                d = min(which(x>FVE_thre)[1], dmax) # modification to ensure we stay below dmax
                col_ind = rep(FALSE,dmax)
                col_ind[1:d] = TRUE
                col_ind
            })
            

            col_keep = do.call(c, col_keep)
            
            print(table(col_keep))
            

            
            Rmat = Rmat[col_keep,col_keep]
            
            
            
            ### deal with small/negative eigenvalues of the corr matrix 
            Rmat_IC_ts = adjust_R_troubleshoot(Rmat)
            Rmat_IC = adjust_R(Rmat)
            ###
            graph_all = list()
            
            ### get sparse corr matrix using GPP method ####
            
            # BIC 
            factor = sqrt(ntrain)
            res_GPP_BIC_v2 = run_GPP_HT_BIC_v2(Rmat=Rmat_IC,
                                               grpind=grpind,
                                               ntrain=ntrain,
                                               factor=factor,
                                               num_edges=min_edges,
                                               ncores=ncores)

            # res_GPP_BIC <- run_GPP_HT_BIC_troubleshoot(Rmat=Rmat_IC,
            #                                   grpind=grpind,
            #                                   ntrain=ntrain,
            #                                   factor=factor,
            #                                   num_edges=min_edges)
            # 
            # res_GPP_BIC_v3 <- run_GPP_HT_BIC_tol(Rmat=Rmat_IC,
            #                                      grpind=grpind,
            #                                      ntrain=ntrain,
            #                                      factor=factor)
            
            print('checkpoint 3')
            
            # print(table(res_GPP_BIC_v2 == res_GPP_BIC))
            # print(table(res_GPP_BIC_v2 == res_GPP_BIC_v3))
            # print(table(res_GPP_BIC == res_GPP_BIC_v3))
            
            graph_all[["GPP_BIC"]] = as.matrix((get_groupNorm(res_GPP_BIC_v2, grpind)!=0)+0)
            graph_all[['missing_neurons']] <- missing_neurons
            
            
            ### save results 
            
            save_dir <- paste('result_simu/', '/', sep = setting_ID2)                       #../result_simu/e1_m2vr0_n400_me1/
            
            # if the folder doesn't exist, create it 
            if (!dir.exists(save_dir)) {
                dir.create(save_dir)
            }    
            
            save_file <- paste(save_dir, '.rda', sep = file_ID2)                            #../result_simu/e1_m2vr0_n400/me1/Tau1_e1_m2vr0_n400_me1.rda
            
            save(graph_all, file=save_file)
            
            # print for job scheduler
            print(paste(ID2[i], ' just finished', sep = ''))
            print(Sys.time())    
        }
    
    }
}
