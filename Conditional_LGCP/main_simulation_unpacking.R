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

# choices: rho_i_dist, rho_ij_dist, g_ij_dist, C_HS, V_HS, P_HS, sens, spec, auc, accuracy

metrics <- c('rho_i_dist', 'rho_ij_dist', 'g_ij_dist', 'C_HS', 'P_HS', 'auc')
grid.newpage()
visualize_metrics('simu_results_banded_c0_5', metrics, 1, 2) %>% print()
# convergence of intermediate estimators





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