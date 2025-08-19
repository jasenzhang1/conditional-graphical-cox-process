# 3) pre-processing before estimation ------------------------------------------

simu_threshold <- dataset$simulation_params$threshold_p

# 3.1) first, convert dataset into a format that can be used for estimation

df_estimate <- convert_data_for_estimation(dataset$subject_data) %>% as.data.table()
df_estimate$time <- df_estimate$time / T_max # normalize to [0, 1]

query_yd <- rep(1, n) %>% as.data.frame()  # assume they are all from the same discrete strata
y_c_strata <- dataset$Y_continuous
query_y_cs <- dataset$Y_continuous %>% unique()
discrete_strata <- query_yd %>% unique()
Tseq_est <- 1:19/20
#Tseq <- time_grid

ncores <- parallel::detectCores() - 1

# 4) estimation ----------------------------------------------------------------


estimated_results <- full_conditional_estimation(df_estimate, y_c_strata, query_y_cs, Tseq_est, simu_threshold, ncores)


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