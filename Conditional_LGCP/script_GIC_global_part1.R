suppressPackageStartupMessages(library(RhpcBLASctl))
blas_set_num_threads(1); omp_set_num_threads(1)
source('functions/00_function_wrapper.R')

args <- commandArgs(trailingOnly = TRUE)
model_type <- args[1] 

if(model_type == 'mice'){
  ID            <- args[2]
  y_c_structure <- args[3]
  time_scale    <- as.numeric(args[4])
  method        <- args[5]
  movement      <- as.numeric(args[6])
  VR            <- as.numeric(args[7])
  cont_inds     <- as.numeric(args[8]) # Total number of queries
  
  discrete_level <- paste0('m', movement, 'vr', VR)
  temp_file_dir  <- 'temp_data/mice'
  GIC_folder     <- paste0(temp_file_dir, '/GIC_joint_', ID, '_', discrete_level, '_t', time_scale)
  mouse <- TRUE
} else if(model_type == 'simu'){
  n_large   <- as.numeric(args[2])
  n         <- as.numeric(args[3])
  rep_i     <- as.numeric(args[4])
  adj_type  <- args[5]
  method    <- args[6]
  cont_inds <- as.numeric(args[7])
  
  temp_file_dir <- 'temp_data/simu'
  GIC_folder    <- paste0(temp_file_dir, '/GIC_joint_', adj_type, '_n_', n, '_rep', rep_i)
  mouse <- FALSE
}

if (!dir.exists(GIC_folder)) dir.create(GIC_folder, recursive = TRUE)

# --- 1. Load Data from all queries ---
if(mouse){
  step_3_paths <- paste0(temp_file_dir, '/part3_', ID, '_', discrete_level, '_t', time_scale, '_nquery', 1:cont_inds, '.rds')
  step_2_paths <- paste0(temp_file_dir, '/part2_', ID, '_', discrete_level, '_t', time_scale, '_nquery', 1:cont_inds, '.rds')
} else{
  step_3_paths <- paste0(temp_file_dir, '/part3_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '_rep_', rep_i, '.rds')
  step_2_paths <- paste0(temp_file_dir, '/part2_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '_rep_', rep_i, '.rds')
}

step_2_data <- lapply(step_2_paths, readRDS)
W_y_list    <- lapply(step_2_data, `[[`, "W_y")
p           <- step_2_data[[1]]$p

results_all_queries <- lapply(step_3_paths, readRDS)

# --- 2. Multi-ID Setup Logic (Task Mapping) ---
# Get names like "KL_cor_est_eig1", "KL_cor_est_eig2"
suffix_names <- step_00_grab_ID(names(results_all_queries[[1]]$step_5b), 'KL_cor')

task_map <- data.frame()
max_k <- 0

for(suffix in suffix_names){
  # Create a list where each element corresponds to one query, 
  # but only contains the matrix for this specific suffix
  C_cond_list_suffix <- lapply(results_all_queries, function(res) res$step_5b[[suffix]])
  
  # Run your existing setup logic, but specific to this suffix
  # Note: Modified GIC_joint_part1_setup to accept a suffix to avoid file overwrites
  num_k_suffix <- GIC_joint_part1_setup(C_cond_list_suffix, p, W_y_list, GIC_folder, suffix)
  
  task_map <- rbind(task_map, data.frame(id = suffix, num_k = num_k_suffix))
  
  if(num_k_suffix > max_k) max_k <- num_k_suffix
}

# --- 3. Output for Bash ---
num_suffixes <- nrow(task_map)
write.csv(task_map, paste0(GIC_folder, "/task_map.csv"), row.names = FALSE)

cat("max_k=", max_k, "\n")
cat("num_suffixes=", num_suffixes, "\n")