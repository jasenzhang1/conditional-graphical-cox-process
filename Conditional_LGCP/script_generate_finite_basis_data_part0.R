
library(RhpcBLASctl)

# limit threads in BLAS/LAPACK
blas_set_num_threads(1)   # limit BLAS
omp_set_num_threads(1)    # limit OpenMP

args <- commandArgs(trailingOnly = TRUE)

n <- as.numeric(args[1])
n_query <- as.numeric(args[2])
adj_type <- args[3]
adj_params <- as.numeric(args[4:length(args)])


source('functions/00_function_wrapper.R')
source('functions/20_simulation_function_wrapper.R')


# n <- 100
# adj_type <- 'block_banded_v2'
# adj_params <- c(0, 1, 0.4, 0.8, 2)
# adj_type <- 'block_banded_c0'
# adj_params <- c(0.5, 0.5, 2)

d <- 2
p <- 12
m <- 30
m_est <- 30
T_max <- 1


beta_0 <- 5
time_grid <- make_time_grid(m)
time_grid_est <- make_time_grid(m_est)
time_grid_both <- sort(union(time_grid, time_grid_est))



y_c_query <- make_time_grid(n_query) %>% matrix(nrow = n_query)

 

# 1) system parameters
seed = NULL
# ncores <- parallel::detectCores() - 1
ncores = 1

# 0) preprocessing

Y_c <- generate_y_c_adj_type(n, adj_type, adj_params, seed = NULL)
basis_list <- trig_basis(d)
mean_vec <- rep(0, p*d)

# 1) generate dataset ----------------------------------------------------------

param_list <- list(
  Y_c = Y_c,
  basis_list = basis_list,
  n = n,
  d = d,
  p = p, 
  adj_type = adj_type,
  adj_params = adj_params, 
  beta_0 = beta_0, 
  time_grid = time_grid, 
  time_grid_est = time_grid_est,
  time_grid_both = time_grid_both,
  T_max = T_max, 
  y_c_query = y_c_query,
  seed = seed,
  mean_vec = mean_vec
)

temp_file_dir <- 'temp_data'
if (!dir.exists(temp_file_dir)) dir.create(temp_file_dir)  # temp_data

temp_file_dir <- paste0(temp_file_dir, '/simu_data')
if (!dir.exists(temp_file_dir)) dir.create(temp_file_dir)  # temp_data/simu_data

# save data

datafile_name <- paste0("part0_", adj_type, '_n_', n, '.rds')

saveRDS(param_list, file = file.path(temp_file_dir, datafile_name)) 
                                          





