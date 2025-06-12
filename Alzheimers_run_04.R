

start_time <- Sys.time()
ffs = 1
set.seed(ffs)
source("GraphPP_FUN.R")
source("modified_funcs.R")
library(dplyr)



#### simulate the data ####

# 1) Open 99_Create_Spiketrain_Dataset and load the data_df file

# 1a) Parameters that don't change

IDs <- c('346', '351', '366', '361', '362', '368')  # mouse ID
ID2 <- c('Tau1', 'Tau2', 'Tau3', 'WT1', 'WT2', 'WT3') # our name 

num_neurons <- c(169, 250, 240, 249, 235, 294)
names(num_neurons) <- ID2

# 1b) Parameters that change, double check!!!

task_num <- '04'

movements <- c(0, 1, 2)
VR <- 0

ews <- 1:4
ew_symb <- 'e'

tn_symb <- 't'
tns <- 10  # replicates or timescale

min_edges <- 1
ncores <- parallel::detectCores() - 1

print('Number of Cores')
print(ncores)



for(ew in ews){ # for each epoch/week
    
    for(movement in movements){ # for each movement
        
        data_root <- paste0("spike_data/task_", task_num, '/') # spike_data/task_04/
        
        setting_ID <- paste0(ew_symb, ew, '_m', movement, 'vr', VR, '_', tn_symb, tns) # e1_m0vr0_t20
        setting_ID2 <- paste0(setting_ID, '_me', min_edges)                            # e1_m0vr0_t20_me1
        
        data_root <- paste0(data_root, setting_ID)                    # ../spike_data/w17_m2vr0_n400/
        
        file_names <- list.files(path = data_root, full.names = TRUE)
        
        m <- length(file_names)
        mice_names <-  sub(".*/([^/_]+)_.*", "\\1", file_names)  # get mice names, between last '/' and very next '_'        
        
        for(i in 1:m){ # for each mouse
            
            start_time_i <- Sys.time()
            
            print('--------------------------------------------------')
            print(paste('week or epoch: ', ew))
            print(paste('movement: ', movement))
            print(paste('mouse: ', mice_names[i]))
            print('\n')
            
            
            load(file_names[i]) # data_df2
            
            time_scale <- tns
            num_replicates <- unique(data_df2$subject_num) %>% length()
            print(paste0('timescale of each replicate: ', time_scale, ' secs'))
            print(paste0('number of replicates: ', num_replicates))
            
            # which neurons fire way too little and must be discarded (under 100 total firings = discard)
            
            
            data_df3 <- data_df2[, if (.N >= 100) .SD, by = feature_id]
            
            # there are some rare cases when gamma is so inflated that it makes calculations insane
            
            f_ids <- c()
            for(f_id in unique(data_df3$feature_id)){
                data_i = data_df3[feature_id==f_id,]
                data_i_count = data_i[,.(count=.N),by="subject_num"]
                
                gamma = get_gamma_reproduce(data_i$time)
                if(gamma > 100){
                    f_ids <- c(f_ids, f_id)
                }
            }
            
            if(length(f_ids) > 0){
                print('GAMMA IS INFLATED, WE REMOVED THESE NEURONS')
                print(f_ids)
            }
            
            data_df4 <- data_df3[!feature_id %in% f_ids]

            
            
            
            neuron_count <- num_neurons[mice_names[i]]
            
            silent_neurons <- setdiff(1:neuron_count,
                                      unique(data_df2$feature_id))
            
            missing_neurons <- setdiff(unique(data_df2$feature_id),
                                       unique(data_df4$feature_id))
            
            inactive_neurons <- c(silent_neurons, missing_neurons)
            
            # print number of neurons
            print(paste('total neurons: ', unname(neuron_count)))                     # total neurons
            print(paste('active neurons: ', length(unique(data_df4$feature_id))))     # active neurons
            print(paste('discarded neurons: ', length(missing_neurons)))              # discarded neurons
            print(paste('silent neurons: ', length(silent_neurons)))                  # silent neurons
            print(paste('do they add up? ', length(unique(data_df4$feature_id)) + length(missing_neurons) + length(silent_neurons) == unname(neuron_count)))
            
            
            # retroactively get parameters
            
            ntrain <- length(unique(data_df4$subject_num))
            p <- length(unique(data_df4$feature_id))
            
            
            ### get corr matrix for full data
            patient_sel = 1:ntrain
            feature_sel = unique(data_df4$feature_id) %>% sort()
            
            Tseq = seq(0.05,0.95,length=19)
            dmax = 4
            FVE_thre = 0.9 # originally 0.9
            
            quantile(data_df4[,(.N),by=c("subject_num","feature_id")]$V1)
            
            # Rmat_diag_full_ts = get_rho_diag_pp_troubleshoot(data_all=data_df4,
            #                                  patient_sel=patient_sel,
            #                                  feature_sel=feature_sel,
            #                                  Tseq=Tseq,
            #                                  dmax=dmax,
            #                                  ncores=ncores)
            
            Rmat_diag_full = get_rho_diag_pp_reproduce(data_all=data_df4,
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
            
            # print(min(diag(Rmat)))
            
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
            
            # print(table(col_keep))
            
            
            
            Rmat = Rmat[col_keep,col_keep]
            
            
            
            ### deal with small/negative eigenvalues of the corr matrix 
            # Rmat_IC_ts = adjust_R_troubleshoot(Rmat)
            Rmat_IC = adjust_R(Rmat)
            ###
            graph_all = list()
            
            ### get sparse corr matrix using GPP method ####
            
            # BIC 
            factor = sqrt(ntrain)
            res_GPP_BIC_v3 = run_GPP_HT_BIC_v3(Rmat=Rmat_IC,
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
            
            graph_all[["GPP_BIC"]] = as.matrix((get_groupNorm(res_GPP_BIC_v3[[1]], grpind)!=0)+0)
            graph_all[["weighted_GPP_BIC"]] = as.matrix((get_groupNorm(res_GPP_BIC_v3[[1]], grpind)))
            graph_all[['missing_neurons']] <- inactive_neurons
            graph_all[['tuning_parameters']] <- c(res_GPP_BIC_v3[[2]], res_GPP_BIC_v3[[3]])  # params 
            graph_all[['tuning_parameter_indices']] <- c(res_GPP_BIC_v3[[4]], res_GPP_BIC_v3[[5]])  # indices 
            graph_all[['tuning_parameter_max_indices']] <- c(res_GPP_BIC_v3[[6]], res_GPP_BIC_v3[[7]])  # max indices 
            

            graph_all[['time_scale']] <- time_scale
            graph_all[['num_replicates']] <- num_replicates
            
            ### save results 
            
            save_dir <- paste0('results_alzheimers/task_', task_num, '/')  # results_alzheimers/task_04/
            if (!dir.exists(save_dir)) {
                dir.create(save_dir)
            }               
            
            save_dir <- paste0(save_dir, setting_ID2, '/')                       #../task_04/e1_m2vr0_t10_me1/
            
            # if the folder doesn't exist, create it 
            if (!dir.exists(save_dir)) {
                dir.create(save_dir)
            }    
            
            file_ID2 <- paste(mice_names[i], setting_ID2, sep = '_')
            
            save_file <- paste(save_dir, '.rda', sep = file_ID2)                            #../result_simu/e1_m2vr0_n400/me1/Tau1_e1_m2vr0_n400_me1.rda
            
            save(graph_all, file=save_file)
            
            # print for job scheduler
            # print(paste(ID2[i], ' just finished', sep = ''))
            
            end_time_i <- Sys.time()
            elapsed <- as.numeric(end_time_i - start_time_i, units = 'mins') %>% round(2)
            print(paste('minutes taken: ', elapsed, sep = ''))
            
        }
    }
}
print('-----------------------------------------------------------')
end_time <- Sys.time()
total_time <- as.numeric(end_time - start_time, units = "mins") %>% round(2)
print('COMPLETE!!')
print(paste('total minutes taken: ', total_time, sep = ''))


