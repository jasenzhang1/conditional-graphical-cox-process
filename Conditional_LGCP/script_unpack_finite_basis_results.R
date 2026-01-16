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

data_folder <- args[1]    
base_folder <- args[2]  
method      <- args[3]    
X_truth     <- args[4]
beta_truth  <- args[5]
adj_type    <- args[6]
adj_params <- as.numeric(args[7:length(args)])


# data_folder <- 'simu_data'
# base_folder <- 'simu_results'
# adj_type <- 'hub_block_c0'
# method <- 'CPGM'
# X_truth <- T
# beta_truth <- T
# eigen_troubleshoot <- T

n_large <- 100
n <- 100


# 3) heatmaps of certain metrics 

exploratory_ids <- c('01', '02', '22', '23', '29', '41', '44', 
                     '50', '51' # check on beta_truth
                     ) 

bivariate_ids <- c('25', '32')  # rho_ij and g_ij

final_ids <- c('58',   # KL_cor assembled, 58b = KL_cov
               '95',   # C_HS 
               '112',  # w_mat
               '113',  # roc
               '114')  # adj_mat


truth_file_name <- paste0(data_folder, '/', adj_type, '_n_', n_large, '_truths.RData')
results_folder <- paste0(base_folder, '/', adj_type, '/', method)
estimates_file_name <- paste0(results_folder, '/', adj_type, '_n_', n, '.RData')

# 28e - visualize results over y_c
g_heatmaps_exploratory <- visualize_finite_basis(truth_file_name, estimates_file_name, exploratory_ids, beta_truth, X_truth, eigen_troubleshoot)
g_heatmaps_bivariate   <- visualize_finite_basis(truth_file_name, estimates_file_name, bivariate_ids,   beta_truth, X_truth, eigen_troubleshoot)
g_heatmaps_final       <- visualize_finite_basis(truth_file_name, estimates_file_name, final_ids,       beta_truth, X_truth, eigen_troubleshoot)

# 3b) gif

g_gif <- visualize_precision_gif(data_folder, p, d, adj_type, adj_params, nframes = 50, fps = 5)

# ------------------------------------------------------------------------------
# 4) convergence of intermediate estimators (28b)

i <- 1
j <- 2
metrics_summary <- visualize_metrics_finite_basis(truth_file_name, results_folder, i, j)

# ------------------------------------------------------------------------------
# 5) print table and figure
results_folder_2 <- paste0(results_folder, '/export')

if (!dir.exists(results_folder_2)) {
  dir.create(results_folder_2)
}


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

  
  
  
  