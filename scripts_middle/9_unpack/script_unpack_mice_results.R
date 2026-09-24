# ------------------------------------------------------------------------------
#
# Per-mouse diagnostic figures for the data analysis (Section 6).
#
# Called by scripts_exec/run_data_analysis.sh after each mouse x stratum fit:
#
#   Rscript scripts_middle/9_unpack/script_unpack_mice_results.R \
#     <ID> <time_scale> <method> <movement> <VR> <eigen_setting> <results folder> [<results folder> ...]
#
# One results folder per lwGIC connectivity floor, e.g.
#   mice_results/paper/week_only_bw_001_min_05/CPGM_BOTH_100_NORMALIZED
# Figures go to <results folder>/export/<ID>_<stratum>_t<time_scale>_*.pdf
#
# ------------------------------------------------------------------------------

source('functions/00_function_wrapper.R')
source('functions/20_simulation_function_wrapper.R')

args <- commandArgs(trailingOnly = TRUE)

ID              <- args[1]                   # ID <- 'Tau3'
time_scale      <- as.numeric(args[2])       # time_scale <- 10   (each replicate is 10 seconds)
method          <- args[3]                   # method <- 'CPGM'
movement        <- as.numeric(args[4])       # movement <- 0
VR              <- as.numeric(args[5])       # VR <- 1
eigen_setting   <- args[6]                   # only_joint, trig_and_joint, trig_simple
results_folders <- args[7:length(args)]      # results_folders <- 'mice_results/paper/week_only_bw_001_min_05/CPGM_BOTH_100_NORMALIZED'

discrete_level <- paste0('m', movement, 'vr', VR)

# 3) heatmaps of certain metrics 

exploratory_ids <- c('01', '02', '22', '23', '29', '41', '44')

bivariate_ids <- c('25', '32')  # rho_ij and g_ij


KL_ids <- c('58',   # KL_cor assembled, 58b = KL_cov
            '59')   # KL_prec


HS_ids <- c('95',   # C_HS, 95b = specific entries only
            '111')  # w_mat

final_ids <- c('114')  # adj_mat

tau_ids <- c('121')  # tau_c and tau_p




for (results_folder in results_folders) {

results_folder_2 <- paste0(results_folder, '/export')
if (!dir.exists(results_folder_2)) dir.create(results_folder_2)

# ------------------------------------------------------------------------------
# 28e - visualize results over y_c

print('Plotting 28e Figures')

beta_truth <- F
X_truth <- F
ground_truth <- F

estimates_file_name <- paste0(results_folder, '/', ID, '_', discrete_level, '_t', time_scale, '.RData')

graph_results_i <- load_file(estimates_file_name)

g_exploratory <- visualize_over_time(graph_results_i, exploratory_ids, ground_truth, beta_truth, X_truth, eigen_setting) #28c
g_bivariate   <- visualize_over_time(graph_results_i, bivariate_ids,   ground_truth, beta_truth, X_truth, eigen_setting) 
g_KL          <- visualize_over_time(graph_results_i, KL_ids,          ground_truth, beta_truth, X_truth, eigen_setting) 
g_HS          <- visualize_over_time(graph_results_i, HS_ids,          ground_truth, beta_truth, X_truth, eigen_setting) 
g_final       <- visualize_over_time(graph_results_i, final_ids,       ground_truth, beta_truth, X_truth, eigen_setting) 
g_tau         <- visualize_over_time(graph_results_i, tau_ids,         ground_truth, beta_truth, X_truth, eigen_setting) 


g_name_exploratory <- paste0(results_folder_2, '/', ID, '_', discrete_level, '_t', time_scale, '_exploratory.pdf')
g_name_bivariate   <- paste0(results_folder_2, '/', ID, '_', discrete_level, '_t', time_scale, '_bivariate.pdf')
g_name_KL          <- paste0(results_folder_2, '/', ID, '_', discrete_level, '_t', time_scale, '_KL.pdf')
g_name_HS          <- paste0(results_folder_2, '/', ID, '_', discrete_level, '_t', time_scale, '_HS.pdf')
g_name_final       <- paste0(results_folder_2, '/', ID, '_', discrete_level, '_t', time_scale, '_final.pdf')
g_name_tau         <- paste0(results_folder_2, '/', ID, '_', discrete_level, '_t', time_scale, '_tau.pdf')

# print regular graphs

reg_graphs <- list(g_exploratory, g_bivariate, g_KL, g_HS, g_final, g_tau)
reg_graph_names <- c(g_name_exploratory, g_name_bivariate, g_name_KL, g_name_HS, g_name_final, g_name_tau)

# 1. Use 1:length() to iterate through all lists
for(k in 1:length(reg_graphs)){
  
  g <- reg_graphs[[k]]
  
  # Ensure the filename has a .pdf extension
  current_filename <- reg_graph_names[k]
  if(!grepl("\\.pdf$", current_filename)) {
    current_filename <- paste0(current_filename, ".pdf")
  }
  
  pdf(current_filename, width = 15, height = 15)
  
  # 2. Iterate through the NAMES of the list 'g'
  for (name in names(g)) {
    
    # 3. Pull the ACTUAL object from the list using the name
    obj <- g[[name]]
    
    # Check if it's a gtable (like g_29)
    if (inherits(obj, "gtable")) {
      grid.newpage()
      grid.draw(obj)
    } 
    # Check if it's a sub-list of plots (like g_50/g_51)
    else if (is.list(obj) && !inherits(obj, "ggplot") && !inherits(obj, "patchwork")) {
      for (sub_obj in obj) {
        if (inherits(sub_obj, "ggplot") || inherits(sub_obj, "patchwork")) {
          print(sub_obj)
        }
      }
    }
    # Handle standard ggplot and patchwork objects
    else {
      # This handles cases where obj might be NULL or non-plot types
      if (inherits(obj, "ggplot") || inherits(obj, "patchwork")) {
        print(obj)
      }
    }
  }
  
  dev.off()
}

}  # results folder


# Figures that compare mice and strata are drawn once every fit is done:
# scripts_figures/data_analysis_figures.R
