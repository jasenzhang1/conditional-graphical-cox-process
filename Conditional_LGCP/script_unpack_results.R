# unpack
library(grid)
library(gridExtra)
library(tidyverse)
source('functions/28_Simulation_Visualization.R')
source('functions/28b_Simulation_Visualization_2.R')
source('functions/00d_debugging.R')


args <- commandArgs(trailingOnly = TRUE)

base_folder <- args[1]
adj_type <- args[2]
method <- args[3]
i <- args[4]
j <- args[5]

# convergence of intermediate estimators

metrics <- c('rho_i_dist', 'rho_ij_dist', 'g_ij_dist', 'C_HS', 'P_HS', 'auc')
results_folder <- paste0(base_folder, '/', adj_type, '/', method)
metrics_summary <- visualize_metrics(results_folder, metrics, i, j)


# print table and figure
write.csv(metrics_summary$metric_table, paste0(results_folder, "/metrics_summary.csv"), row.names = FALSE)


pdf(paste0(results_folder, "metrics_summary.pdf"), width = 8, height = 6)  # open PDF file
grid.arrange(metrics_summary$metric_graph)                                 # draw the grob/layout
dev.off()                                                                  # close the file