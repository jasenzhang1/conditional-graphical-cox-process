
# script to create mice replicate (y_c, y_d) data for LGCP estimation

t0 <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)


ID <- args[1]
time_scale <- as.numeric(args[2])


library(dplyr)
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
t_start <- Sys.time()
LGCP_data <- get_spiketrain_dataset_conditional_LGCP(ID, time_scale, discrete_covariates, continuous_covariates, weekly_dataset, neuron_df)

# prepare to save
file_name <- paste0(save_dir, '/', ID2, '_t', time_scale, '_data.rda')

# save 
save(LGCP_data, file = file_name)

# report time taken
t_end <- Sys.time()
t_elapsed <- difftime(t_end, t_start, units = 'mins') %>% as.numeric() %>% round(2)
print(paste0(ID2, '_t', time_scale, ' finished in ', t_elapsed, ' mins'))

