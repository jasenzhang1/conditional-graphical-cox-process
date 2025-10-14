# unpack
library(grid)
library(gridExtra)
library(tidyverse)
source('functions/28_Simulation_Visualization.R')
source('functions/28b_Simulation_Visualization_2.R')
source('functions/00d_debugging.R')

# parameters

T_max = 1                                           # Time horizon (T)
time_grid_size = 50                                 # Time discretization (m)

time_grid <- seq(0, T_max, length.out = time_grid_size)
time_grid_est <- 1:19/20
time_grid_both <- sort(union(time_grid, time_grid_est))
p <- 10

n <- 500
query_id <- 1

# only look at one dataset
unpacking_pipeline('simu_results/banded_c2/OG', '112', time_grid_est, time_grid, n, p, query_id)

# convergence of intermediate estimators

metrics <- c('rho_i_dist', 'rho_ij_dist', 'g_ij_dist', 'C_HS', 'P_HS', 'auc')
grid.newpage()
metrics_summary <- visualize_metrics('simu_results/banded_trig2/OG', metrics, 1, 2)
View(metrics_summary$metric_table)
print(grid.arrange(metrics_summary$metric_graph))



# 2) ground truth vs estimated p x p matrix




visualize_truths_from_est(step_list_temp, '44', time_grid_est, time_grid, time_grid_both, p)

visualize_truths_from_est(step_list_OG, '113', time_grid_est, time_grid, time_grid_both, p)
visualize_truths_from_est(step_list_OG, '112', time_grid_est, time_grid, time_grid_both, p)

visualize_truths_from_est(step_list_CPGM, '113', time_grid_est, time_grid, time_grid_both, p)
visualize_truths_from_est(step_list_CPGM, '112', time_grid_est, time_grid, time_grid_both, p)

# troubleshoot KL coeffs

eigen_OG <- graph_results_OG$estimated_graphs_part_1$step_4$eigen_decomp_est



eigen_CPGM <- graph_results_CPGM$estimated_graphs_part_2[[1]]$step_4$eigen_decomp_est

KL_OG <- graph_results_OG$estimated_graphs_part_1$step_5$kl_coeffs_est
var_OG <- apply(KL_OG, c(2,3), var)

cov_CPGM <- graph_results_CPGM$estimated_graphs_part_2[[1]]$step_5$KL_cov_est
var_CPGM_1 <- diag(cov_CPGM[['1_1']])





# visualize ||V - V_Hat||_HS error
grid.newpage()
visualize_V_cond_convergence('simu_results_banded_c0_5', 'C', 1, 2) %>% print()

# visualize AUC versus n and y_c
grid.newpage()
visualize_AUC_across_n('simu_results_banded_c1_5') %>% print()



