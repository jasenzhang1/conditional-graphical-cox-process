source('functions/28z_Visualization_helpers.R')
source('functions/28y_Visualization_Blocks.R')
source('functions/28_Simulation_Visualization.R')
source('functions/28c_Grand_Visualizations.R')

visualize_finite_basis <- function(truth_file_name, estimates_file_name, graph_ids, ground_truth, beta_truth, X_truth, eigen_setting){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: visualize intermediate plots of the finite basis simulation setting 
  #
  # - just merging `truths` and `estimates` from the finite basis setting
  # - then using `visualize over time` from 28_c
  #
  # inputs:
  # 
  # - truth_file_name       (string)              'simu_data/block_banded_v2_n_100_truths.RData'
  # - estimates_file_name   (string)              'simu_results/block_banded_v2/CPGM/CPGM_n_100.RData'
  # - graph_ids             (vector of strings)   which graphs do we want to see
  # - i
  # - j
  # - ground_truth          (boolean)             do our results have the overall truth?
  # - beta_truth            (boolean)             do our results have beta_truth values? 
  # - X_truth               (boolean)             do our results have X_truth values? 
  # - eigen_setting         (string)              only_joint, trig_and_joint, trig_simple
  #
  # outputs:
  #
  # - g_comparisons    (list)  list of graphs that display heatmaps (estimate, truth) over queries
  #
  # ----------------------------------------------------------------------------
  
  # 1) merge truths and estimates
  
  truths <- load_file(truth_file_name)
  estimates <- load_file(estimates_file_name)
  
  #merged <- convergence_metrics_part1(truth_file_name, estimates_file_name)
  merged <- merge_lists_recursive(truths, estimates)

  # 2) now, utilize the code in 28c

  g_comparisons <- visualize_over_time(merged, graph_ids, ground_truth, beta_truth, X_truth, eigen_setting)
  
  return(g_comparisons)
  

  
  
}


visualize_metrics_CI <- function(results_folder, n_reps, adj_type) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Visualize final prediction metrics across multiple reps 
  #
  # 
  # inputs:
  #
  # - results_folder    (string)    simu_results/hub_block_v2/CPGM
  # - n_reps            (integer)   how many reps?
  # - adj_type          (string)    hub_block_v2
  #
  # 
  #
  #
  # ----------------------------------------------------------------------------
  
  results_list <- list()
  
  # 1) Loop through each replicate folder
  for (i in 1:n_reps) {
    current_rep_path <- file.path(results_folder, paste0("rep", i))
    
    # Use the logic from our previous discussion to get Ns in this specific folder
    files <- list.files(current_rep_path, pattern = "\\.RData$", full.names = TRUE)
    
    # Extract N from the filename: "adj_type_n_X_rep_i.RData"
    # We use a regex that looks specifically for the pattern in your file names
    ns_in_folder <- as.numeric(sub(".*_n_(\\d+)_rep_.*", "\\1", files))
    
    # 2) Loop through each file (each n) within the rep folder
    for (j in seq_along(files)) {
      
      obj_name <- load(files[j])
      res <- get(obj_name)
      
      
      # Access step_12
      # Assuming the loaded object is named 'res' - change to your actual object name
      s12 <- res$step_12 
      y_c_names <- names(s12)
      
      for (y_query in y_c_names) {
        current_y_item <- s12[[y_query]]
        
        # Find items starting with "roc_KL_"
        suffix_items <- grep("^roc_KL_", names(current_y_item), value = TRUE)
        
        for (item_name in suffix_items) {
          suffix <- sub("roc_KL_", "", item_name)
          metrics <- current_y_item[[item_name]]
          
          # Create a single row for this specific combination
          row_data <- data.frame(
            adj_type = adj_type,
            n = ns_in_folder[j],
            y_c_query = y_query,
            suffix = suffix,
            rep = i,
            accuracy = metrics$accuracy,
            f1_score = metrics$f1_score,
            sensitivity = metrics$sensitivity,
            specificity = metrics$specificity,
            ppv = metrics$ppv,
            npv = metrics$npv,
            auc = metrics$auc,
            stringsAsFactors = FALSE
          )
          results_list[[length(results_list) + 1]] <- row_data
        }
      }
    }
  }
  
  # Combine everything into one master dataframe
  full_df <- do.call(rbind, results_list)
  full_df$n <- as.factor(full_df$n) # Factorize for better x-axis grouping
  
  # SAVE AS CSV
  csv_path <- paste0(results_folder, '/', adj_type, "_metrics_CI_results.csv")
  write.csv(full_df, csv_path, row.names = FALSE)
  
  # 3) Pivot Long for ggplot (one row per metric)
  plot_df <- full_df %>%
    pivot_longer(cols = c(accuracy, f1_score, sensitivity, specificity, ppv, npv, auc),
                 names_to = "metric_name",
                 values_to = "value")
  
  # 4) Plotting
  # Each metric gets a figure
  
  # 3) Generate PDF with one page per metric


  metrics_vec <- unique(plot_df$metric_name)
  
  # Increased height (14) is crucial for vertical stacking to look good
  pdf(paste0(results_folder, '/', adj_type, "_metrics_CI_plots.pdf"), width = 10, height = 14) 
  
  for (m in metrics_vec) {
    p <- plot_df %>%
      filter(metric_name == m) %>%
      ggplot(aes(x = n, y = value, fill = suffix)) +
      geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
      # This exact syntax (Variable ~ .) is what forces the vertical stack
      facet_grid(y_c_query ~ ., labeller = label_both) + 
      coord_cartesian(ylim = c(0, 1)) +
      labs(
        title = paste("Metric Analysis:", m),
        subtitle = paste("Adjustment Type:", adj_type),
        x = "Sample Size (n)",
        y = m,
        fill = "Suffix"
      ) +
      theme_bw() + # Better borders for stacked facets
      theme(
        legend.position = "bottom",
        # Rotates the vertical labels so they are readable (horizontal)
        strip.text.y = element_text(angle = 0, face = "bold"),
        # Adds space between the stacked rows
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = "gray95")
      )
    
    print(p)
  }
  
  dev.off()
  
}