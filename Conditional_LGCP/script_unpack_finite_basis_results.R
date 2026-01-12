# main script to unpack finite basis results

# 1) packages and functions 

library(grid)
library(gridExtra)
library(dplyr)
library(patchwork)
library(ggplot2)

source('functions/00_function_wrapper.R')
source('functions/28_Simulation_Visualization.R')
source('functions/28b_Simulation_Visualization_2.R')
source('functions/28c_Grand_Visualizations.R')
source('functions/00d_debugging.R')
source('functions/28e_Finite_Basis_Visualization.R')

# 2) arguments for simulation settings

args <- commandArgs(trailingOnly = TRUE)

data_folder <- args[1]    # data_folder <- 'simu_data'
base_folder <- args[2]    # base_folder <- 'simu_results'
adj_type    <- args[3]    # adj_type <- 'block_banded_v2'
method      <- args[4]    # method <- 'CPGM'
X_truth     <- args[5]
beta_truth  <- args[6]

# data_folder <- 'simu_data'
# base_folder <- 'simu_results'
# adj_type <- 'flexible_block_banded_c0'
# method <- 'CPGM'
# X_truth <- T
# beta_truth <- T


n_large <- 1500
n <- 1000

# 3) heatmaps of certain metrics 

graph_ids <- c('11', '12',
               '22',    # rho_i
               '25',    # rho_ij [1,2]
               '32',    # g_ij [1,2]
               '112',   # w_mat  (P_HS)
               '113')   # roc

truth_file_name <- paste0(data_folder, '/', adj_type, '_n_', n_large, '_truths.RData')
estimates_file_name <- paste0(base_folder, '/', adj_type, '/', method, '/', adj_type, '_n_', n, '.RData')

# 28e - visualize results over y_c
g_heatmaps <- visualize_finite_basis(truth_file_name, estimates_file_name, graph_ids, beta_truth, X_truth)

# ------------------------------------------------------------------------------
# 4) convergence of intermediate estimators (28b)

i <- 1
j <- 2
results_folder <- paste0(base_folder, '/', adj_type, '/', method)
metrics_summary <- visualize_metrics_finite_basis(truth_file_name, results_folder, i, j)

results_folder_2 <- paste0(results_folder, '/export')

if (!dir.exists(results_folder_2)) {
  dir.create(results_folder_2)
}

# print table and figure
write.csv(metrics_summary$metric_table, paste0(results_folder_2, "/metrics_summary.csv"), row.names = FALSE)


pdf(paste0(results_folder_2, "/point_metrics_summary.pdf"), width = 8, height = 6)        # open PDF file
grid.arrange(metrics_summary$point_metrics_graph)                                 # draw the grob/layout
dev.off()                                                                         # close the file

pdf(paste0(results_folder_2, "/eval_metrics_summary.pdf"), width = 8, height = 6)        # open PDF file
grid.arrange(metrics_summary$eval_metrics_graph)                                  # draw the grob/layout
dev.off()                                                                         # close the file

# 2) ground truth vs estimated p x p matrix

# max_n <- max(get_ns(results_folder))
# 
# n_y_c_query <- length(get_y_c_query(results_folder, max_n))
# 
# prec_mat_comparison <- visualize_prec_mat_over_time(results_folder, max_n)
# pdf(paste0(results_folder, "/prec_mat_over_time.pdf"), width = 3 * (n_y_c_query-2), height = 3)  # open PDF file
# grid.arrange(prec_mat_comparison$graph)                                                            # draw the grob/layout
# dev.off()  

  
  
  
  