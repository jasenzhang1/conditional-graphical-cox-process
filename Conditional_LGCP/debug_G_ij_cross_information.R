# incorporate G_ij into estimation


t0 <- Sys.time()

source('functions/20_simulation_function_wrapper.R')
source('functions/00_function_wrapper.R')

n_large <- 100


p = 10                      # Number of processes (p)  
T_max = 1                   # Time horizon (T)
q_c = 1                     # Continuous conditioning dimension (q_c)
y_c_borders = list(1:4)     # Border values
K = 4                       # Discrete combinations (K)
sparsity = 0.2              # Graph sparsity (s)
theta = 1.0                 # Signal strength (theta)
adj_type = "banded_v1"        # Graph topology
dependence_type = "constant"  # Conditional dependence type
time_grid_size = 50         # Time discretization (m)
seed = 1
ncores = parallel::detectCores() - 1
terse = TRUE


# m x m GP kernel parameters

m <- time_grid_size

base_kernel_params <- list(base_gamma = 20,      # won't be pd, but we can massage it 
                           base_kernel = 'rbf',
                           base_variance = 1,
                           base_GP_mean = 5)

time_grid <- seq(0, T_max, length.out = time_grid_size)
time_grid_est <- 1:19/20

time_grid_both <- sort(union(time_grid, time_grid_est))

# 1) generate dataset ----------------------------------------------------------

dataset <- simulate_conditional_cox_data_v4(n_large, p, T_max, q_c, y_c_borders,
                                            sparsity,
                                            theta,
                                            dependence_type,
                                            adj_type,
                                            time_grid,
                                            time_grid_est,
                                            base_kernel_params,
                                            seed,
                                            ncores)

# 2) prepare for estimation

time_grid <- dataset$time_grid
time_grid_est <- dataset$time_grid_est
time_grid_both <- dataset$time_grid_both
true_graphs <- dataset$true_graphs
true_graph_indices <- sort(names(true_graphs)) 

n <- dim(dataset$Y_continuous)[1]


# create dataset for inference
data_df4 <- convert_data_for_estimation(dataset$subject_data) %>% as.data.table()
data_df4$time <- data_df4$time / T_max # normalize to [0, 1]

query_yd <- rep(1, n) %>% as.data.frame()  # assume they are all from the same discrete strata
y_c_strata <- dataset$Y_continuous
query_y_cs <- dataset$Y_continuous %>% unique()
discrete_strata <- query_yd %>% unique()
Tseq_est <- time_grid_est




# size of dataset
print(paste0('number of subjects: ', n))
print(paste0('number of spikes: ', nrow(data_df4)))


# 3.3) retroactively get parameters for subjects (ntrain, patient_sel) and neurons (p, feature_sel)

patient_sel = unique(data_df4$subject_num) %>% sort()
feature_sel = unique(data_df4$feature_id) %>% sort()


# 4) estimation --------------------------------------------------------------


rho_list_1 <- estimate_intensities_stratum_parallel_v3(data_df4, patient_sel, feature_sel, Tseq_est, ncores)
rho_list_2 <- estimate_intensities_stratum_parallel_v4(data_df4, patient_sel, feature_sel, Tseq_est, T, ncores)

g_ij_1   <- estimate_covariance_functions_ii(rho_list_1)
g_ij_2   <- estimate_covariance_functions_ii_cross_informed(rho_list_2)

wow1 <- g_ij_1[,,1]
wow2 <- g_ij_2[,,1]



grid.arrange(visualize_nonneg_matrix_heatmap(g_ij_1[,,1]),
             visualize_nonneg_matrix_heatmap(g_ij_2[,,1]),
             visualize_nonneg_matrix_heatmap(dataset$true_graphs$`1`$P_block_kronecker$base_cov), nrow = 1)
