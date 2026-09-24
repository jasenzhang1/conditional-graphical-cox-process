
library(RhpcBLASctl)
source('functions/00_function_wrapper.R')
source('functions/20_simulation_function_wrapper.R')

# limit threads in BLAS/LAPACK
blas_set_num_threads(1)   # limit BLAS
omp_set_num_threads(1)    # limit OpenMP

args <- commandArgs(trailingOnly = TRUE)

p        <- as.numeric(args[1])
d        <- as.numeric(args[2])
n_large  <- as.numeric(args[3])
rep_i    <- as.numeric(args[4])
n_query  <- as.numeric(args[5])
beta_0   <- as.numeric(args[6])
adj_type <- args[7]
adj_params <- as.numeric(args[8:length(args)])




# p <- 12
# d <- 2
# n_large <- 100
# adj_type <- 'block_banded_v2'
# adj_params <- c(0, 1, 0.4, 0.8, 2)
# adj_type <- 'block_banded_c0'
# adj_params <- c(0.5, 0.5, 2)

# 1) intermediate parameters

m <- 30
m_est <- 30
T_max <- 1


time_grid <- make_time_grid(m)
time_grid_est <- make_time_grid(m_est)
time_grid_both <- sort(union(time_grid, time_grid_est))
m_both <- length(time_grid_both)


y_c_query <- make_y_c_grid(n_query) %>% matrix(nrow = n_query)



# 1) system parameters
seed <- simulation_seed(adj_type, rep_i, stage = 0)
cat('seed (part 0) =', seed, '\n')
# ncores <- parallel::detectCores() - 1
ncores = 1

# 0) preprocessing

Y_c <- generate_y_c_adj_type(n_large, adj_type, adj_params, seed = seed)
basis_list <- trig_basis(d)
mean_vec <- rep(0, p*d) # mean for generating beta



mu_t        <- rep(beta_0, m)    # mu(t) = beta_0 (constant)
mu_t_both   <- rep(beta_0, m_both)
mu_t_coarse <- rep(beta_0, m_est)

# 1) generate dataset ----------------------------------------------------------

param_list <- list(
  Y_c = Y_c,                     # (n x q_c matrix)
  basis_list = basis_list,       # (d-dim list)
  n = n_large,                   # integer
  d = d,                         # integer
  p = p,                         # integer
  m = m,                         # integer
  m_est = m_est,                 # integer
  m_both = m_both,               # integer
  adj_type = adj_type,           # string
  adj_params = adj_params,       # vector of parameters
  beta_0 = beta_0,               # scalar
  time_grid = time_grid,           # (m-dim vector)
  time_grid_est = time_grid_est,   # (m_est dim vector)
  time_grid_both = time_grid_both,  # (m_both-dim vector)
  T_max = T_max,                    # integer
  y_c_query = y_c_query,            # (n_query x q_c matrix)
  seed = seed,                      # integer
  mean_vec = mean_vec,              # (pd-dim vector)
  mu_t = mu_t,                      # (m-dim vector)
  mu_t_both = mu_t_both,            # (m_both-dim vector)
  mu_t_coarse = mu_t_coarse         # (m_est-dim vector)
)

temp_file_dir <- 'temp_data'
if (!dir.exists(temp_file_dir)) dir.create(temp_file_dir)  # temp_data

temp_file_dir <- paste0(temp_file_dir, '/simu_data_', adj_type, '_n_', n_large, '_rep_', rep_i)
if (!dir.exists(temp_file_dir)) dir.create(temp_file_dir)  # temp_data/simu_data

# save data

datafile_name <- paste0("part0_", adj_type, '_n_', n_large, '_rep_', rep_i, '.rds')

saveRDS(param_list, file = file.path(temp_file_dir, datafile_name)) 
                                          





