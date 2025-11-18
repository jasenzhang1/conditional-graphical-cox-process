# unpack
library(grid)
library(gridExtra)
library(dplyr)
source('functions/28_Simulation_Visualization.R')
source('functions/28b_Simulation_Visualization_2.R')
source('functions/28c_Grand_Visualizations.R')
source('functions/00d_debugging.R')
source('functions/28e_Finite_Basis_Visualization.R')

args <- commandArgs(trailingOnly = TRUE)

base_folder <- args[1]    # base_folder <- 'simu_results'
adj_type <- args[2]       # adj_type <- 'block_banded_v2'
method <- args[3]         # method <- 'CPGM'

data_folder <- 'simu_data'
n_large <- 1024
n <- 512

# heatmaps of certain metrics 

graph_ids <- c('22', '25', '31', '91')

truth_file_name <- paste0(data_folder, '/', adj_type, '_n_', n_large, '_truths.RData')
estimates_file_name <- paste0(base_folder, '/', adj_type, '/', method, '/', method, '_n_', n, '.RData')
estimates_file_name <- paste0(base_folder, '/', adj_type, '/', method, '/', method, '_n_', n, '_finer.RData')

g_heatmaps <- visualize_finite_basis(truth_file_name, estimates_file_name, graph_ids)

# convergence of intermediate estimators

metrics <- c('rho_i_dist', 'rho_ij_dist', 'g_ij_dist', 'C_HS', 'C_HS_v2', 'P_HS', 'auc')
results_folder <- paste0(base_folder, '/', adj_type, '/', method)
metrics_summary <- visualize_metrics(results_folder, metrics, i, j)


# print table and figure
write.csv(metrics_summary$metric_table, paste0(results_folder, "/metrics_summary.csv"), row.names = FALSE)


pdf(paste0(results_folder, "/metrics_summary.pdf"), width = 8, height = 6)  # open PDF file
grid.arrange(metrics_summary$metric_graph)                                 # draw the grob/layout
dev.off()                                                                  # close the file


# 2) ground truth vs estimated p x p matrix

max_n <- max(get_ns(results_folder))

n_y_c_query <- length(get_y_c_query(results_folder, max_n))

prec_mat_comparison <- visualize_prec_mat_over_time(results_folder, max_n)
pdf(paste0(results_folder, "/prec_mat_over_time.pdf"), width = 3 * (n_y_c_query-2), height = 3)  # open PDF file
grid.arrange(prec_mat_comparison$graph)                                                            # draw the grob/layout
dev.off()  

  
  
  
  