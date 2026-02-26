# main script to unpack finite basis results



source('functions/00_function_wrapper.R')
source('functions/20_simulation_function_wrapper.R')

# 2) arguments for simulation settings

args <- commandArgs(trailingOnly = TRUE)

n_large  <- as.numeric(args[1])
ns       <- as.numeric(unlist(strsplit(args[2], " ")))
n_reps      <- as.numeric(args[3])
method      <- args[4]    
X_truth     <- as.logical(args[5])
beta_truth  <- as.logical(args[6])
eigen_setting  <- args[7]
adj_type    <- args[8]
adj_params <- as.numeric(args[9:length(args)])


# n_large <- 10000
# n <- 10000
# n_reps <- 50
# method <- 'CPGM'
# X_truth <- T
# beta_truth <- T
# eigen_troubleshoot <- T
# adj_type <- 'hub_block_v2'
# adj_params <- c()

 
data_folder <- paste0('simu_data/', adj_type, '_n_', n_large, '_rep_1') # simu_data/hub_block_v2_n_100_rep_1
base_folder <- 'simu_results'

eigen_troubleshoot <- (eigen_setting == 'trig_and_joint')

ground_truth <- T
make_gif <- F

# 3) heatmaps of certain metrics 

exploratory_ids <- c('01', '02', '22', '23', '29', '41', '44')

bivariate_ids <- c('25', '32')  # rho_ij and g_ij

beta_ids <- c('50', '51')

KL_ids <- c('58',   # KL_cor assembled, 58b = KL_cov
            '59')   # KL_prec


HS_ids <- c('95',   # C_HS, 95b = specific entries only
            '111',  # w_mat
            '112')  # w_mat

final_ids <- c('113',  # roc
               '114',  # adj_mat
               '120')  # performance metrics

tau_ids <- c('121')  # tau_c and tau_p

truth_file_name <- paste0(data_folder, '/', adj_type, '_n_', n_large, '_rep_1_truths.RData')
results_folder <- paste0(base_folder, '/', adj_type, '/', method, '/rep1')  # simu_results/hub_block_v2/CPGM/rep1


results_folder_2 <- paste0(results_folder, '/export')

if (!dir.exists(results_folder_2)) {
  dir.create(results_folder_2)
}

# ------------------------------------------------------------------------------
# 28e - visualize results over y_c

print('Plotting 28e Figures')

for(n in ns){
  

  estimates_file_name <- paste0(results_folder, '/', adj_type, '_n_', n, '_rep_1.RData')
  
  g_exploratory <- visualize_finite_basis(truth_file_name, estimates_file_name, exploratory_ids, ground_truth, beta_truth, X_truth, eigen_troubleshoot)
  g_betas       <- visualize_finite_basis(truth_file_name, estimates_file_name, beta_ids,        ground_truth, beta_truth, X_truth, eigen_troubleshoot)
  g_bivariate   <- visualize_finite_basis(truth_file_name, estimates_file_name, bivariate_ids,   ground_truth, beta_truth, X_truth, eigen_troubleshoot)
  g_KL          <- visualize_finite_basis(truth_file_name, estimates_file_name, KL_ids,          ground_truth, beta_truth, X_truth, eigen_troubleshoot)
  g_HS          <- visualize_finite_basis(truth_file_name, estimates_file_name, HS_ids,          ground_truth, beta_truth, X_truth, eigen_troubleshoot)
  g_final       <- visualize_finite_basis(truth_file_name, estimates_file_name, final_ids,       ground_truth, beta_truth, X_truth, eigen_troubleshoot)
  g_tau         <- visualize_finite_basis(truth_file_name, estimates_file_name, tau_ids,         ground_truth, beta_truth, X_truth, eigen_troubleshoot)
  
  g_name_exploratory <- paste0(results_folder_2, '/', adj_type, '_n_', n, '_exploratory.pdf')
  g_name_betas       <- paste0(results_folder_2, '/', adj_type, '_n_', n, '_betas.pdf')
  g_name_bivariate   <- paste0(results_folder_2, '/', adj_type, '_n_', n, '_bivariate.pdf')
  g_name_KL          <- paste0(results_folder_2, '/', adj_type, '_n_', n, '_KL.pdf')
  g_name_HS          <- paste0(results_folder_2, '/', adj_type, '_n_', n, '_HS.pdf')
  g_name_final       <- paste0(results_folder_2, '/', adj_type, '_n_', n, '_final.pdf')
  g_name_tau         <- paste0(results_folder_2, '/', adj_type, '_n_', n, '_tau.pdf')
  
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
  
  # print beta graphs
  pdf(g_name_betas, width = 8, height = 5)
  
  # 2. Iterate through g_50 and g_51 within g_betas
  for (beta_name in c("g_50", "g_51")) {
    
    # Extract the sub-list (e.g., g_betas$g_50)
    sub_list <- g_betas[[beta_name]]
    
    # 3. Iterate through each covariate plot (0.125, 0.375, 0.625, 0.875)
    for (covariate_name in names(sub_list)) {
      
      plot_obj <- sub_list[[covariate_name]]
      
      # Check if it is a valid plot object before printing
      if (inherits(plot_obj, "ggplot") || inherits(plot_obj, "patchwork")) {
        
        # Optional: Add a title so you know which beta and covariate you are looking at
        plot_to_print <- plot_obj + 
          ggtitle(paste("Source:", beta_name, "| Covariate:", covariate_name))
        
        print(plot_to_print)
      }
    }
  }
  
  # 4. Close the device
  dev.off()
  
  
}

# -------------------------------------
# 3b) gif

if(make_gif){
  
  
  # adj_suffixes <- c('c2', 'v2', 'j2')
  # adj_1 <- as.list(paste0('flexible_block_banded_', adj_suffixes))
  # adj_2 <- as.list(paste0('hub_block_', adj_suffixes))
  # adj_3 <- as.list(paste0('complete_block_', adj_suffixes))
  # adj_type_list <- list(adj_1, adj_2, adj_3)
  # 
  # adj_params_1 <- list(c(0, 1, 2, 2, 0.7, 0.7),
  #                      c(0, 1, 2, 2, 0.5, 0.9, 0.5, 0.9),
  #                      c(0, 1, 0.5, 2, 2, 0.7, 0.7))
  # adj_params_2 <- list(c(0, 1, 4, 2, 2, 0.7, 0.7),
  #                      c(0, 1, 4, 2, 2, 0.5, 0.9, 0.5, 0.9),
  #                      c(0, 1, 4, 0.5, 2, 2, 0.7, 0.7))
  # adj_params_3 <- list(c(0, 1, 4, 2, 2, 0.7, 0.7),
  #                      c(0, 1, 4, 2, 2, 0.5, 0.9, 0.5, 0.9),
  #                      c(0, 1, 4, 0.5, 2, 2, 0.7, 0.7))
  # 
  # adj_params_list <- list(adj_params_1, adj_params_2, adj_params_3)
  
  
  matrix_type <- 'prec_norm'
  p <- 12
  d <- 2
  
  # x_lab_vec <- c('Constant', 'Linear', 'Jump')
  # y_lab_vec <- c('Banded', 'Hub', 'Complete')
  
  adj_type_list <- list(adj_type)
  adj_params_list <- list(adj_params)
  
  g_gif <- visualize_precision_gif(data_folder, p, d, adj_type, adj_params, nframes = 50, fps = 5)
  
  # saved in data_folder = 'simu_data'
  g_grid_gif <-       visualize_precision_grid_gif(data_folder, p, d, adj_type_list, adj_params_list, matrix_type, nframes = 50, fps = 5)
  g_grid_gif_v2 <- visualize_precision_grid_gif_v2(data_folder, p, d, adj_type_list, adj_params_list, matrix_type, nframes = 10, fps = 2)
}




# ------------------------------------------------------------------------------
# 4) convergence of intermediate estimators (28b)

print('Plotting Convergences')

i <- 1
j <- 2
metrics_summary <- visualize_metrics_finite_basis(truth_file_name, results_folder, i, j)


# 4b) print table and figure

print('Saving Results')

# raw dataframe
write.csv(metrics_summary$metric_df, paste0(results_folder_2, "/metrics_raw.csv"), row.names = FALSE)  

# stratified table - rds
saveRDS(metrics_summary$metric_tables, paste0(results_folder_2,  "/metrics_tabular.rds"))

# stratified table - csv
flat_table <- metrics_summary$metric_tables %>% 
  unnest(cols = c('row_metric', 'subrow_n')) 

write.csv(flat_table, paste0(results_folder_2, "/metrics_tabular.csv"), row.names = FALSE)

# pdfs 
pdf(paste0(results_folder_2, "/point_metrics_summary.pdf"), width = 12, height = 9)        # open PDF file
grid.arrange(metrics_summary$point_metrics_graph)                                 # draw the grob/layout
dev.off()                                                                         # close the file

pdf(paste0(results_folder_2, "/eval_metrics_summary.pdf"), width = 8, height = 6)        # open PDF file
grid.arrange(metrics_summary$eval_metrics_graph)                                  # draw the grob/layout
grid.arrange(metrics_summary$eval_metrics_graph2)                                  
dev.off()                                                                         # close the file

# ------------------------------------------------------------------------------
# 5) 28e) CI's for accuracy/f1/sens/spec etc for each n across reps
  
print('Getting CI Results')

results_folder <- paste0(base_folder, '/', adj_type, '/', method)

visualize_metrics_CI(results_folder, n_reps, adj_type)

