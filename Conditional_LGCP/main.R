ffs = 1
set.seed(ffs)


source('functions/01_kernel_estimation.R')
source('functions/02_intensity_estimation.R')

library(dplyr)
library(data.table)

load('data/WT1_w19_demo_data.rda')

ncores <- parallel::detectCores() - 1
ncores <- 4
print(paste0('number of cores: ', ncores))

data_df <- LGCP_data[[1]]
y_d_df <- LGCP_data[[2]]
y_c_df <- LGCP_data[[3]]

# extract the binary levels
cov_df <- y_d_df %>% dplyr::select(- subject_num)

discrete_strata <- cov_df %>% unique()

min_edges <- 0

for(y_ind in c(1)){ # 1) for each discrete variable level
  
  yd_i <- discrete_strata[y_ind, ]
  s_yd <- y_d_df$subject_num[apply(cov_df, 1, function(x) all(x == yd_i))]  # s_yd = indices
  data_df2 <- data_df %>% filter(subject_num %in% s_yd)
  
  data_df4 <- data_df2
  
  neuron_count <- 249
  
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
  patient_sel = unique(data_df4$subject_num) %>% sort()
  feature_sel = unique(data_df4$feature_id) %>% sort()
  
  #Tseq = seq(0.01,0.99,length=99)
  Tseq = seq(0.05,0.95,length=19)
  
  dmax = 4
  FVE_thre = 0.9 # originally 0.9  
  
  # estimation of rho
  
  time_1 <- Sys.time()
  rho_list <- estimate_intensities_stratum_parallel(data_df4, patient_sel, feature_sel, Tseq, ncores)
  time_2 <- Sys.time()
  
  t_algo <- as.numeric(difftime(time_2, time_1, units = "mins")) %>% round(2)
  
  print(paste0('parallel algorithm took ', t_algo, ' minutes'))
  
}