# unpack
library(grid)
library(gridExtra)
library(tidyverse)
source('functions/28_Simulation_Visualization.R')

# 1) graphs that don't change over time 
# 2) graphs that change over time, but values below adjacency threshold are 0
# 3) ...

# 10) continuous covariates = (week, velocity, VR)

# simu_results = constant

# unpack and then graph (28b)
step_list_OG <- list(step_1 = graph_results_OG$estimated_graphs_part_1$step_1,
                  step_2 = graph_results_OG$estimated_graphs_part_1$step_2,
                  step_3 = graph_results_OG$estimated_graphs_part_1$step_3,
                  step_4 = graph_results_OG$estimated_graphs_part_1$step_4,
                  step_5 = graph_results_OG$estimated_graphs_part_1$step_5,
                  step_6 = NA,
                  step_7 = NA,
                  step_8 = graph_results_OG$estimated_graphs_part_2[[1]]$step_8,
                  step_9 = graph_results_OG$estimated_graphs_part_2[[1]]$step_9,
                  step_10 = graph_results_OG$estimated_graphs_part_2[[1]]$step_10,
                  step_11 = graph_results_OG$estimated_graphs_part_2[[1]]$step_11,
                  step_12 = graph_results_OG$estimated_graphs_part_2[[1]]$step_12,
                  step_2b = graph_results_OG$estimated_graphs_part_1$step_2b
                  )

step_list_CPGM <- list(step_1 = graph_results_CPGM$estimated_graphs_part_1$step_1,
                      step_2 = graph_results_CPGM$estimated_graphs_part_2[[1]]$step_2,
                      step_3 = graph_results_CPGM$estimated_graphs_part_2[[1]]$step_3,
                      step_4 = graph_results_CPGM$estimated_graphs_part_2[[1]]$step_4,
                      step_5 = graph_results_CPGM$estimated_graphs_part_2[[1]]$step_5,
                      step_6 = NA,
                      step_7 = NA,
                      step_8 = NA,
                      step_9 = graph_results_CPGM$estimated_graphs_part_2[[1]]$step_9,
                      step_10 = graph_results_CPGM$estimated_graphs_part_2[[1]]$step_10,
                      step_11 = graph_results_CPGM$estimated_graphs_part_2[[1]]$step_11,
                      step_12 = graph_results_CPGM$estimated_graphs_part_2[[1]]$step_12,
                      step_2b = graph_results_CPGM$estimated_graphs_part_2[[1]]$step_2b
)

step_list <- list(step_1 = graph_results_JASA$estimated_graphs_part_1$step_1,
                  step_2 = graph_results_JASA$estimated_graphs_part_1$step_2,
                  step_3 = graph_results_JASA$estimated_graphs_part_1$step_3,
                  step_4 = graph_results_JASA$estimated_graphs_part_1$step_4,
                  step_5 = graph_results_JASA$estimated_graphs_part_1$step_5,
                  step_6 = NA,
                  step_7 = NA,
                  step_8 = graph_results_JASA$estimated_graphs_part_2[[1]]$step_8,
                  step_9 = graph_results_JASA$estimated_graphs_part_2[[1]]$step_9,
                  step_10 = graph_results_JASA$estimated_graphs_part_2[[1]]$step_10,
                  step_11 = graph_results_JASA$estimated_graphs_part_2[[1]]$step_11,
                  step_12 = graph_results_JASA$estimated_graphs_part_2[[1]]$step_12,
                  step_2b = graph_results_JASA$estimated_graphs_part_1$step_2b
)

step_list_temp <- list(step_1 = step_1,
                  step_2 = step_2,
                  step_3 = step_3,
                  step_4 = step_4,
                  step_5 = NA,
                  step_6 = NA,
                  step_7 = NA,
                  step_8 = NA,
                  step_9 = NA,
                  step_10 = NA,
                  step_11 = NA,
                  step_12 = NA,
                  step_2b = step_2b
)


T_max = 1                                           # Time horizon (T)
time_grid_size = 50                                 # Time discretization (m)

time_grid <- seq(0, T_max, length.out = time_grid_size)
time_grid_est <- 1:19/20
time_grid_both <- sort(union(time_grid, time_grid_est))
p <- 10

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

# unpack the graphs 
load('simu_results_banded_c1_5/n_500.RData')

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_1$g_01) # ground truths

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_1$g_11) # X_hat_k_i

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2[[1]]$g_22) # rho_i

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2[[1]]$g_25) # rho_ij

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2[[1]]$g_94) # C_cond pm x pm

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2[[1]]$g_104) # P_cond pm x pm

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2[[1]]$g_111b) #HS norms vs adj

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2[[1]]$g_112) #HS norms

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2[[2]]$g_112b) #HS norms

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2[[2]]$g_113) #ROC curves





# visualize ||V - V_Hat||_HS error
grid.newpage()
visualize_V_cond_convergence('simu_results_banded_c0_5', 'C', 1, 2) %>% print()

# visualize AUC versus n and y_c
grid.newpage()
visualize_AUC_across_n('simu_results_banded_c1_5') %>% print()

# convergence of intermediate estimators

metrics <- c('rho_i_dist', 'rho_ij_dist', 'g_ij_dist', 'C_HS', 'P_HS', 'auc')
grid.newpage()
visualize_metrics('simu_results_banded_c1_5', metrics, 1, 2) %>% print()






# values are -10, -7.5, -5.2, -1.95, 4.2
load('simu_results_banded_v2_3/n_100.RData')
grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_1$g_22)

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`4.2`$g_112)

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`4.2`$g_113)

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`-1.95`$g_112)

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`-1.95`$g_113)



grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`-5.2`$g_112)

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`-5.2`$g_113)


grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`-7.5`$g_112)

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`-7.5`$g_113)



grid.newpage()
pdf('temp2.pdf', width = 20, height = 14)
grid.draw(graph_results_i$estimated_graphs_part_2$`4`$g_84)
dev.off()