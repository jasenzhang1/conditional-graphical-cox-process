library(RhpcBLASctl)
blas_set_num_threads(1); omp_set_num_threads(1)
source('functions/00_function_wrapper.R')


args <- commandArgs(trailingOnly = TRUE)
model_type <- args[1] 

if(model_type == 'mice'){
  ID <- args[2]                       # ID <- 'Tau1'
  y_c_structure <- args[3]            # y_c_structure <- "week_only" or "time_and_week"
  time_scale <- as.numeric(args[4])   # time_scale <- 10   (each replicate is 5 seconds)
  method <- args[5]                   # method <- 'CPGM'
  movement <- as.numeric(args[6])     # movement <- 0
  VR <- as.numeric(args[7])           # VR <- 0
  cont_ind <- as.numeric(args[8])     # cont_ind <- 1  
  eigen_setting <- args[9]            # eigen_seting <- 'only_joint'
  
  X_truth <- F
  
  discrete_level <- paste0('m', movement, 'vr', VR)
  
  setting_info_list <- list(ID = ID,
                            y_c_structure = y_c_structure,
                            time_scale = time_scale,
                            method = method,
                            discrete_level = discrete_level)
  
  temp_file_dir <- 'temp_data/mice'
  mouse <- T
  
} else if(model_type == 'simu'){
  
  n_large <- as.numeric(args[2])
  n <- as.numeric(args[3])
  rep_i <- as.numeric(args[4])
  adj_type <- args[5]
  method <- args[6]
  X_truth <- as.logical(args[7])
  cont_ind <- as.numeric(args[8])
  eigen_setting <- args[9]
  
  setting_info_list <- list(n_large = n_large,
                            n = n,
                            rep_i = rep_i,
                            adj_type = adj_type,
                            method = method)
  
  temp_file_dir <- paste0('temp_data/simu_', adj_type, '_n_', n_large, '_rep_', rep_i)
  mouse <- F
  
} else{
  stop('model_type not supported')
}

# Construct file name to load the results from your previous estimation step

if(mouse){
  file_name <- paste0('part2b_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.rds')
  GIC_folder <- paste0(temp_file_dir, '/GIC_local_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind)
  if (!dir.exists(GIC_folder)) dir.create(GIC_folder)  
} else{
  file_name <- paste0("part2b_", adj_type, '_n_', n, '_nquery', cont_ind, '_rep_', rep_i, '.rds')
  GIC_folder <- paste0(temp_file_dir, '/GIC_local_', adj_type, '_n_', n, '_nquery', cont_ind, '_rep', rep_i)
  if (!dir.exists(GIC_folder)) dir.create(GIC_folder)  
}

# --- 1. Load the Estimation Data (.rds safe) ---

full_path <- paste0(temp_file_dir, "/", file_name)

if (!file.exists(full_path)) stop("File not found: ", full_path)

# Handle .rds vs .RData
if (grepl("\\.rds$", file_name, ignore.case = TRUE)) {
  estimated_graphs <- readRDS(full_path)
} else {
  load(full_path) # Assumes it contains 'estimated_graphs'
}

# Ensure variables (step_5b, p, W_y) are standalone
if (exists("estimated_graphs")) {
  list2env(estimated_graphs, envir = .GlobalEnv)
} else {
  stop("Object 'estimated_graphs' not found in loaded file.")
}

# --- 2. Multi-ID Setup Logic ---
# Identify which items in step_5b need GIC (e.g., 'est', 'X_truth')
core_names <- step_00_grab_ID(names(step_5b), 'KL_cor')
input_names <- names(step_5b)

task_map <- data.frame()
max_k <- 0

cat('checkpoint 1')

for(i in seq_along(core_names)){
  C_cond_i <- step_5b[[input_names[i]]]
  id_i <- core_names[i]
  
  # Calculate quantiles and save ID-specific initial data
  # This function saves 'GIC_local_initial_data_[id_i].RData'
  num_k_i <- GIC_step1_precompute(C_cond_i, p, W_y, GIC_folder, id_i)
  
  task_map <- rbind(task_map, data.frame(id = id_i, num_k = num_k_i))
  
  # Keep track of the largest k range for the Bash loop
  if(num_k_i > max_k) max_k <- num_k_i
}

# --- 3. Save Task Map and Return num_k ---


num_suffixes <- nrow(task_map)
  
write.csv(task_map, paste0(GIC_folder, "/task_map.csv"), row.names = FALSE)

# We return max_k so Bash knows the maximum range it might need to loop through

cat("max_k=", max_k, "\n")
cat("num_suffixes=", num_suffixes, "\n")
