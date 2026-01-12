

# ------------------------------------------------------------------------------
#
# GOAL: script to create mice replicate (y_c, y_d) data for LGCP estimation
# 
#
#       and save the data as "LGCP_data" in the file: data/with_ts/ID2_t5_data.rda
#
# inputs:
# 
# - ID           (string)   mouse 3-digit ID
# - time_scales  (integer)  how long is each data replicate? 2, 5, or 10 seconds
# 
#
# variables:
#
# - ID2          (string)  nickname like "Tau1"
# ------------------------------------------------------------------------------



args <- commandArgs(trailingOnly = TRUE)


ID <- args[1]
time_scale <- as.numeric(args[2])


source('functions/40_mice_data_functions.R')
data_dir <- 'data/'



weekly_dataset <- load_realigned_dataset(data_dir)[[1]]
neuron_df <- load_brain_region_dataset(data_dir)[[1]]

ID2 <- as.character(neuron_df$ID2[which(neuron_df$Mouse == ID)][1])  

discrete_covariates <- c('movement', 'VR')
continuous_covariates <- c('week', 'timestamp')

save_dir <- "data/with_ts"
if(!dir.exists(save_dir)){
  dir.create(save_dir, recursive = TRUE)
}


    
# obtain data
LGCP_data <- get_spiketrain_dataset_conditional_LGCP(ID, time_scale, discrete_covariates, continuous_covariates, weekly_dataset, neuron_df)

# prepare to save
file_name <- paste0(save_dir, '/', ID2, '_t', time_scale, '_data.RData')

# save 
save(LGCP_data, file = file_name)


