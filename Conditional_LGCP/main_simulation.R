# 1) generate simulated data
# 2) estimate ground truth
# 3) calculate accuracy metrics

source('functions/20_simulation_function_wrapper.R')
source('functions/00_function_wrapper.R')

n = 100                     # Sample size (n)
p = 10                      # Number of processes (p)  
T_max = 1                   # Time horizon (T)
q_c = 1                     # Continuous conditioning dimension (q_c)
y_c_borders = list(1:4)     # Border values
K = 4                       # Discrete combinations (K)
sparsity = 0.2              # Graph sparsity (s)
theta = 1.0                 # Signal strength (theta)
graph_type = "random"       # Graph topology
dependence_type = "constant"  # Conditional dependence type
time_grid_size = 50         # Time discretization (m)
seed = 1
ncores = parallel::detectCores() - 1
  
# m x m GP kernel parameters

m <- time_grid_size

base_kernel_params <- list(base_gamma = 0.2 * m^2,  # chatgpt said this would ensure rbf is pd
                           base_kernel = 'rbf',
                           base_variance = 1,
                           base_GP_mean = 5)

time_grid <- seq(0, T_max, length.out = time_grid_size)


dataset <- simulate_conditional_cox_data_v4(n, p, T_max, q_c, y_c_borders,
                                            sparsity,
                                            theta,
                                            dependence_type,
                                            time_grid,
                                            base_kernel_params,
                                            seed,
                                            ncores)


# visualize intensities --------------------------------------------------------

# mean of exp(log_intensities)
mat_list <- lapply(dataset$subject_data, function(x) exp(x$X_functions))
rho_simu <- Reduce("+", mat_list) / length(mat_list)

g_data_2 <- visualize_log_intensity(rho_simu, time_grid) # average intensities across all subjects

# but we also have ground truth of GP mean

m <- time_grid_size
GP_means_m <- dataset$true_graphs[[1]]$P_block_kronecker$GP_simu_mean
GP_means_pm <- rep(GP_means_m, p)

GP_kernel_pm <- dataset$true_graphs[[1]]$P_block_kronecker$GP_simu_var

GP_means_exp_pm <- exp(GP_means_pm + 0.5 * diag(GP_kernel_pm))

GP_mean_ground_truth <- matrix(GP_means_exp_pm, nrow = p, ncol = m)

g_GP_baseline <- visualize_log_intensity(GP_mean_ground_truth, time_grid)

grid.arrange(g_data_2, g_GP_baseline, nrow = 1)



# pre-processing ---------------------------------------------------------------

simu_threshold <- dataset$simulation_params$threshold_p

# first, convert dataset into a format that can be used for estimation

df_estimate <- convert_data_for_estimation(dataset$subject_data) %>% as.data.table()
df_estimate$time <- df_estimate$time / T_max # normalize to [0, 1]

query_yd <- rep(1, n) %>% as.data.frame()  # assume they are all from the same discrete strata
y_c_strata <- dataset$Y_continuous
query_y_cs <- dataset$Y_continuous %>% unique()
discrete_strata <- query_yd %>% unique()
Tseq <- 1:19/20

ncores <- parallel::detectCores() - 1

# estimation -------------------------------------------------------------------


estimated_results <- full_conditional_estimation(df_estimate, y_c_strata, query_y_cs, Tseq, simu_threshold, ncores)


# visualize rho_ij -------------------------------------------------------------

rho_ij_simu <- ground_truth_rho_ij(dataset$subject_data)
rho_11_simu <- rho_ij_simu$`1_1`
rho_11_est <- estimated_results$rho_list[[1]]$rho_ii_mat

df_rho_11_simu <- reshape2::melt(rho_11_simu)
df_rho_11_est <- reshape2::melt(rho_11_est)

df_rho_11_simu$heatmap <- 'Simulation'
df_rho_11_est$heatmap <- 'Estimation'

# df_rho_11_simu$Var1 <- time_grid[df_rho_11_simu$Var1]
# df_rho_11_est$Var1 <- Tseq[df_rho_11_est$Var1]
# df_rho_11_simu$Var2 <- time_grid[df_rho_11_simu$Var2]
# df_rho_11_est$Var2 <- Tseq[df_rho_11_est$Var2]


df_rho_11 <- rbind(df_rho_11_simu,  df_rho_11_est)

# Plot with facet_wrap
ggplot(df_rho_11, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  facet_wrap(~heatmap) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(strip.text = element_text(size = 14))


# visualize intensities (estimate) ---------------------------------------------

rho_i_est <- estimated_results$intensities

g_rho <- visualize_log_intensity(rho_i_est, Tseq)

grid.arrange(g_data_2, g_rho)


# compare Precision matrices (yc = 1) ------------------------------------------

simu_mats <- dataset$true_graphs[[1]]$P_block_kronecker
#simu_mats <- dataset$true_graphs

# P_simu <- kronecker(simu_mats$H_mat_prec_thresh, simu_mats$base_precision)

P1 <- solve(simu_mats$base_cov)
P1 <- (P1 + t(P1))/2

P2 <- (simu_mats$prec_mat + t(simu_mats$prec_mat))/2

P_simu <- kronecker(P2, P1) # ground truth pm x pm matrix (500 x 500)

P_est <- estimated_results$estimates[[1]]$P_mat #estimated pm x pm matrix (190 x 190)

delta_t_simu <- time_grid[2] - time_grid[1]
delta_t_est <- Tseq[2] - Tseq[1]

#P_simu_HS <- sqrt(sum(simu_mats$base_precision^2)) * delta_t_simu * abs(simu_mats$H_mat_prec_thresh)
P_simu_HS <- block_matrix_HS(P_simu, 10) * delta_t_simu
P_est_HS <- block_matrix_HS(P_est, 10) * delta_t_est
  


# finally, calculate accuracy metrics

ground_truth <- dataset$true_graphs
estimated_graphs <- estimated_results$estimates


graph_id_1 <- names(dataset$true_graphs)
graph_id_2 <- names(estimated_results$estimates)


id_both <- intersect(graph_id_1, graph_id_2)


for(id in id_both){
  adj_i <- ground_truth[[id]]$adj_mat
  
  adj_j <- list_to_adj_mat(estimated_graphs[[id]]$edge_strengths, p)
  
  print(mean(adj_j[adj_i == 1]))
  print(mean(adj_j[adj_i == 0]))
  
  
}
