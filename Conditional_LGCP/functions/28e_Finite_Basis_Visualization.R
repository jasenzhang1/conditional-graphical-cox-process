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



visualize_accuracy_CI_across_yc <- function(results_folder, mode, n_reps, adj_type, vert_dashed_line){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Visualize accuracy across reps at each sample size. Allow adding a vertically dashed line for the jump cases
  #
  # 
  # inputs:
  #
  # - results_folder    (string)    simu_results/hub_block_v2/CPGM
  # - mode              (string)    'method', 'local', 'hybrid', 'global'
  # - n_reps            (integer)   how many reps?
  # - adj_type          (string)    hub_block_v2
  # - vert_dashed_line  (boolean)   add a vertical dashed line at x = 0.5?
  #
  # 
  #
  #
  # ----------------------------------------------------------------------------
  
  results_list <- list()
  
  # 1) Data Extraction
  for (i in 1:n_reps) {
    current_rep_path <- file.path(results_folder, paste0("rep", i))
    if(!dir.exists(current_rep_path)) next
    
    files <- list.files(current_rep_path, pattern = "\\.RData$", full.names = TRUE)
    
    for (f in seq_along(files)) {
      obj_name <- load(files[f])
      res <- get(obj_name)
      s12 <- res$step_12 
      
      # Extract n from filename
      n_val <- as.numeric(sub(".*_n_(\\d+)_rep_.*", "\\1", files[f]))
      
      for (y_query in names(s12)) {
        current_y_item <- s12[[y_query]]
        suffix_items <- grep("^roc_KL_", names(current_y_item), value = TRUE)
        
        for (item_name in suffix_items) {
          suffix <- sub("roc_KL_", "", item_name)
          metrics <- current_y_item[[item_name]]
          
          results_list[[length(results_list) + 1]] <- data.frame(
            n = n_val,
            y_yc = as.numeric(y_query), 
            suffix = suffix,
            rep = i,
            accuracy = metrics$accuracy,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  full_df <- bind_rows(results_list) %>%
    filter(suffix != 'est_eig1') %>%
    mutate(suffix = recode(suffix, 
                           "GIC_local_est_eig1" = "Local", 
                           "GIC_hybrid_est_eig1" = "Hybrid",
                           "GIC_global_est_eig1" = "Global"))
  
  # ----------------------------------------------------------------------------
  # 2) Plotting Logic based on MODE
  # ----------------------------------------------------------------------------
  
  if (mode == "method") {
    # ORIGINAL LOGIC: One plot per 'n', comparing suffixes (Local/Hybrid/Global)
    n_values <- sort(unique(full_df$n))
    
    for (current_n in n_values) {
      plot_data <- full_df %>% filter(n == current_n)
      
      p <- ggplot(plot_data, aes(x = y_yc, y = accuracy, color = suffix, fill = suffix)) +
        geom_smooth(method = "loess", alpha = 0.3, linewidth = 1.2, level = 0.95, se = TRUE) +
        labs(title = paste("Accuracy at n =", current_n), x = "Time", y = "Accuracy")
      
      p <- apply_beamer_theme(p, vert_dashed_line)
      
      file_name <- paste0(adj_type, "_accuracy_n", current_n, ".png")
      ggsave(file.path(results_folder, file_name), plot = p, width = 8, height = 6, dpi = 300)
    }
    
  } else {
    # NEW LOGIC: One plot for a specific method, comparing different 'n' values
    target_suffix <- tools::toTitleCase(mode) # e.g., 'local' -> 'Local'
    
    plot_data <- full_df %>% 
      filter(suffix == target_suffix) %>%
      mutate(n_factor = as.factor(n)) # n must be discrete for the color scale
    
    p <- ggplot(plot_data, aes(x = y_yc, y = accuracy, color = n_factor, fill = n_factor)) +
      geom_smooth(method = "loess", alpha = 0.2, size = 1.2, level = 0.95, se = TRUE) +
      labs(x = "Time", y = "Accuracy", color = "Sample Size (n)", fill = "Sample Size (n)")
    
    p <- apply_beamer_theme(p, vert_dashed_line)
    
    file_name <- paste0(adj_type, "_", mode, "_comparison_across_n.png")
    ggsave(file.path(results_folder, file_name), plot = p, width = 8, height = 6, dpi = 300)
  }
  
  message("Visualization complete for mode: ", mode)
}

visualize_accuracy_CI_across_yc_faceted <- function(results_folder, mode, n_reps, row_names, col_names, vert_dashed_line, output_folder, fig_title = NULL){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Visualize loess accuracy of all 6 settings in one figure!!!
  #
  # 
  # inputs:
  #
  # - results_folder    (list)      named list of chr vectors: list(Linear=c(...), Jump=c(...))
  # - mode              (string)    'method', 'local', 'hybrid', 'global'
  # - n_reps            (integer)   how many reps?
  # - row_names         e.g. c("Linear", "Jump")
  # - col_names         e.g. c("Banded", "Hub", "Complete")
  # - vert_dashed_line  (list)      nmaed list of logical vectors: add a vertical dashed line at x = 0.5?
  # - output_folder      where to save the combined figure
  # - fig_title = NULL   optional overall title
  # 
  # 
  # ----------------------------------------------------------------------------

  
  results_list <- list()
  
  # --------------------------------------------------------------------------
  # 1) Data Extraction — iterate over all 6 (row x col) combinations
  # --------------------------------------------------------------------------
  for (row_idx in seq_along(row_names)) {
    row_label <- row_names[row_idx]
    
    for (col_idx in seq_along(col_names)) {
      col_label   <- col_names[col_idx]
      folder      <- results_folder[[row_idx]][col_idx]
      vdl         <- vert_dashed_line[[row_idx]][col_idx]
      
      for (i in 1:n_reps) {
        current_rep_path <- file.path(folder, paste0("rep", i))
        if (!dir.exists(current_rep_path)) next
        
        files <- list.files(current_rep_path, pattern = "\\.RData$", full.names = TRUE)
        
        for (f in seq_along(files)) {
          obj_name <- load(files[f])
          res      <- get(obj_name)
          s12      <- res$step_12
          
          n_val <- as.numeric(sub(".*_n_(\\d+)_rep_.*", "\\1", files[f]))
          
          for (y_query in names(s12)) {
            current_y_item <- s12[[y_query]]
            suffix_items   <- grep("^roc_KL_", names(current_y_item), value = TRUE)
            
            for (item_name in suffix_items) {
              suffix  <- sub("roc_KL_", "", item_name)
              metrics <- current_y_item[[item_name]]
              
              results_list[[length(results_list) + 1]] <- data.frame(
                n         = n_val,
                y_yc      = as.numeric(y_query),
                suffix    = suffix,
                rep       = i,
                accuracy  = metrics$accuracy,
                row_label = row_label,
                col_label = col_label,
                vdl       = vdl,          # carry per-panel flag through
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    }
  }
  
  full_df <- bind_rows(results_list) %>%
    filter(suffix != "est_eig1") %>%
    mutate(
      suffix    = recode(suffix,
                         "GIC_local_est_eig1"  = "Local",
                         "GIC_hybrid_est_eig1" = "Hybrid",
                         "GIC_global_est_eig1" = "Global"),
      # Enforce display order for facet axes
      row_label = factor(row_label, levels = row_names),
      col_label = factor(col_label, levels = col_names)
    )
  
  full_df <- full_df %>%
    mutate(xintercept = ifelse(vdl, 0.5, NA_real_))
  

  
  # --------------------------------------------------------------------------
  # 3) Plotting logic based on MODE
  # --------------------------------------------------------------------------
  if (mode == "method") {
    p <- ggplot(full_df, aes(x = y_yc, y = accuracy, color = suffix, fill = suffix)) +
      geom_smooth(method = "loess", alpha = 0.3, size = 1.2, level = 0.95, se = TRUE)
    
  } else {
    target_suffix <- tools::toTitleCase(mode)
    
    plot_data <- full_df %>%
      filter(suffix == target_suffix) %>%
      mutate(n_factor = as.factor(n))
    
    p <- ggplot(plot_data, aes(x = y_yc, y = accuracy, color = n_factor, fill = n_factor)) +
      geom_smooth(method = "loess", alpha = 0.2, size = 1.2, level = 0.95, se = TRUE) +
      labs(color = "Sample Size (n)", fill = "Sample Size (n)")
  }
  
  # --------------------------------------------------------------------------
  # 4) Add shared layers, facet, and theme
  # --------------------------------------------------------------------------
  p <- p +
    geom_vline(
      aes(xintercept = xintercept),
      linetype  = "dashed",
      color     = "gray60",
      size      = 0.8,
      na.rm     = TRUE
    ) +
    facet_grid(
      rows = vars(row_label),
      cols = vars(col_label)
    ) +
    labs(
      title = fig_title,
      x     = "Time",
      y     = "Accuracy"
    )
  
  p <- apply_beamer_theme_faceted(p) 
  

  
  # --------------------------------------------------------------------------
  # 5) Save
  # --------------------------------------------------------------------------
  file_name <- paste0(mode, "_accuracy_faceted_2x3.png")
  ggsave(
    file.path(output_folder, file_name),
    plot   = p,
    width  = 10,   # wider to accommodate 3 columns
    height = 7,    # taller to accommodate 2 rows
    dpi    = 300
  )
  
  message("Faceted visualization complete for mode: ", mode)
}


# Helper function to maintain consistent Beamer aesthetics
apply_beamer_theme <- function(p, vert_dashed_line) {
  p <- p +
    scale_y_continuous(limits = c(0.4, 1.05), breaks = seq(0.4, 1, 0.1)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    theme_classic() +  
    theme(
      legend.position = c(0.05, 0.05),
      legend.justification = c(0, 0),
      legend.title = element_text(size = 14, face = "bold"),
      legend.background = element_rect(fill = alpha("white", 0.5), color = NA),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.8),
      text = element_text(size = 18) 
    )
  
  if (vert_dashed_line) {
    p <- p + geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", alpha = 0.5)
  }
  return(p)
}

apply_beamer_theme_faceted <- function(p, vert_dashed_line = FALSE) {
  p <- p +
    scale_y_continuous(breaks = seq(0.6, 1.0, 0.2)) +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) + 
    theme_classic() +
    theme(
      panel.border         = element_rect(color = "black", fill = NA, size = 0.8),
      axis.line            = element_blank(),
      legend.position      = "bottom",
      legend.justification = "center",
      legend.title         = element_text(size = 14, face = "bold"),
      legend.background    = element_rect(fill = alpha("white", 0.5), color = NA),
      legend.margin        = margin(t = 0),
      panel.grid.major     = element_blank(),
      panel.grid.minor     = element_blank(),
      text                 = element_text(size = 18)
    ) +
    coord_cartesian(xlim = c(-0.02, 1.02), ylim = c(0.58, 1.05), expand = FALSE)
  
  if (vert_dashed_line) {
    p <- p + geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray60", alpha = 0.5)
  }
  
  return(p)
}

visualize_accuracy_heatmap_across_yc <- function(results_folder, truth_file, adj_type, mode = 'local', eigen_setting = "trig_simple") {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Heatmaps with dots to visualize accuracy for a single rep at different sample sizes. Select one mode only.
  #
  # 
  # inputs:
  #
  # - results_folder    (string)    result folder of this specific replicate
  # - truth_file        (string)    name of truth file for this particular replicate
  # - adj_type          (string)    hub_block_v2
  # - mode              (string)    'method', 'local', 'hybrid', 'global'
  # - eigen_setting     (string)    useful to get eig_est_x
  #
  # 
  # outputs:
  #
  # - none, just save a png
  #
  # ----------------------------------------------------------------------------
  
  
  
  # ----------------------------------------------------------------------------
  # 1. File Discovery and Sorting
  # ----------------------------------------------------------------------------
  estimate_files <- list.files(results_folder, 
                               pattern = "\\.RData$", 
                               full.names = TRUE)
  
  truths <- load_file(truth_file)
  
  
  # Helper to get Ns (Assumes get_ns_with_rep_unsorted is defined in your environment)
  original_ns <- get_ns_with_rep_unsorted(results_folder)
  
  # Index + Sort to ensure n=250 comes before n=500, etc.
  idxs <- order(original_ns)
  ns_sorted <- original_ns[idxs]
  estimate_files <- estimate_files[idxs]
  
  # Filter for the specific mode (local, hybrid, global)
  # This pattern matches the internal suffix in your naming convention
  mode_pattern <- switch(tolower(mode),
                         "local"  = "GIC_local_",
                         "hybrid" = "GIC_hybrid_",
                         "global" = "GIC_global_",
                         stop("Mode must be 'local', 'hybrid', or 'global'"))
  
  # ----------------------------------------------------------------------------
  # 2. Main Processing Loop
  # ----------------------------------------------------------------------------
  step_12b_v4 <- list()
  
  for (i in seq_along(estimate_files)) {
    current_n <- ns_sorted[i]
    obj_name <- load(estimate_files[i])
    res <- get(obj_name)
    
    merged <- merge_lists_recursive(truths, res)
    # Access the specific list item requested
    step_12b <- merged$step_12b 
    if (is.null(step_12b)) next
    
    # Mapping eigen_setting to regex pattern
    suffix_pattern <- switch(eigen_setting,
                             "trig_simple"    = "est_eig1$",
                             "trig_and_joint" = "est_eig[123]$",
                             "only_joint"     = "est_eig3$",
                             "est_eig[123]$")
    
    # Identify indices
    step_12b_names <- names(step_12b[[1]])
    suffixes       <- sub("^adj_mat_", "", step_12b_names)
    
    # Filter by the requested MODE and the EIGEN pattern
    target_indices <- grep(paste0(mode_pattern, suffix_pattern), suffixes)
    truth_indices  <- which(grepl("truth$", suffixes))
    
    kept_indices   <- sort(unique(c(truth_indices, target_indices)))
    
    # Map names and transform
    kept_names <- step_12b_names[kept_indices]
    
    # Rename the mode to indicate the sample size: "n = 250"
    n_label <- paste0("n = ", current_n)
    name_map <- setNames(rep(n_label, length(kept_names)), kept_names)
    
    # Rename truth specifically if it needs to stay distinct, 
    # but based on prompt we rename the mode for comparison
    if (length(truth_indices) > 0) {
      truth_orig_name <- step_12b_names[truth_indices]
      name_map[truth_orig_name] <- "truth"
    }
    
    # Rename items of step_12b (the y queries) as "Time: 0.2", and also "Time: 0.0" as opposed to "t = 0"
    names(step_12b) <- paste0("Time: ", sprintf("%.1f", as.numeric(names(step_12b))))
    
    # Filter and Rename the model types inside each y-item
    step_12b_v2 <- lapply(step_12b, function(x) {
      x_filtered <- x[kept_indices]
      names(x_filtered) <- name_map[names(x_filtered)]
      return(x_filtered)
    })
    
    # Final filtering for non-GIC estimates (those renamed to "n = X")
    v2_names <- names(step_12b_v2[[1]])
    est_v2_names <- v2_names[v2_names != "none"] # Keep n-labels and truth
    
    step_12b_v3 <- lapply(step_12b_v2, function(x) x[names(x) %in% est_v2_names])
    
    # Store these estimates in the global aggregation object
    # We combine lists across different n values
    if (length(step_12b_v4) == 0) {
      step_12b_v4 <- step_12b_v3
    } else {
      # Merge based on y-query names
      for (y_name in names(step_12b_v3)) {
        step_12b_v4[[y_name]] <- c(step_12b_v4[[y_name]], step_12b_v3[[y_name]])
        # Ensure 'truth' is only present once in the combined list per y-query
        # We keep the first truth encountered or re-filter
        step_12b_v4[[y_name]] <- step_12b_v4[[y_name]][!duplicated(names(step_12b_v4[[y_name]]))]
      }
    }
  }
  
  # ----------------------------------------------------------------------------
  # 3. Generate Heatmap
  # ----------------------------------------------------------------------------
  # Plug the combined information into the heatmap function
  # This will now display the comparison of different n's against the truth
  
  ordering_vec <- paste0('n = ', ns_sorted)
  p <- result_heatmap_mismatch_x(step_12b_v4, 'truth', 1.2, 14, rm_diag = TRUE, zmid = 0, ordering_vec = ordering_vec)
  
  file_name <- paste0(adj_type, "_", mode, "_heatmap_combined_n.png")
  export_folder <- paste0(results_folder, '/export')
  if (!dir.exists(export_folder)) {
    dir.create(export_folder)
  }
  
  full_path <- file.path(export_folder, file_name)
  
  # Open PNG device
  png(filename = full_path, width = 12, height = 8, units = "in", res = 300)
  print(p)
  dev.off()
  
  message("Saved aggregated heatmap to: ", full_path)
  
}
