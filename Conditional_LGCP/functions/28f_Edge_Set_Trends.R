plot_edge_proportion_all_mice <- function(results_folder, time_scale, discrete_levels) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: For each discrete stratum, plot edge proportion over weeks for all
  #       available mice discovered automatically from results_folder.
  #       Returns a named list of ggplots, one per discrete level.
  #
  # inputs:
  #
  # - results_folder   (string)
  # - time_scale       (integer)  e.g. 10
  # - discrete_levels  (vector of strings)  e.g. c('m0vr0', 'm1vr1', ...)
  #
  # returns: named list of ggplot objects, one per discrete level
  #
  # ----------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Discover all IDs for a given discrete level
  # --------------------------------------------------------------------------
  discover_IDs <- function(d_level) {
    pattern   <- paste0("^(.+)_", d_level, "_t", time_scale, "\\.RData$")
    all_files <- list.files(results_folder, pattern = pattern, full.names = FALSE)
    ids       <- sub(paste0("_", d_level, "_t", time_scale, "\\.RData$"), "", all_files)
    tau_ids   <- sort(ids[grepl("^Tau", ids)])
    wt_ids    <- sort(ids[grepl("^WT",  ids)])
    other     <- sort(ids[!grepl("^(Tau|WT)", ids)])
    c(tau_ids, wt_ids, other)
  }
  
  # --------------------------------------------------------------------------
  # Load edge proportions for one mouse across all discrete levels
  # --------------------------------------------------------------------------
  load_proportions <- function(ID) {
    
    result <- list()
    
    for (d_level in discrete_levels) {
      
      file_path <- paste0(results_folder, '/', ID, '_', d_level, '_t', time_scale, '.RData')
      
      if (!file.exists(file_path)) next
      
      tmp_env <- new.env()
      load(file_path, envir = tmp_env)
      res_i <- tmp_env$graph_results_i
      
      step_data <- res_i$step_12b
      y_c_weeks <- as.numeric(res_i$y_c_query)
      proportions <- numeric(length(y_c_weeks))
      
      for (idx in seq_along(y_c_weeks)) {
        adj_mat <- step_data[[idx]][["adj_mat_KL_GIC_local_est_eig3"]]
        if (is.null(adj_mat)) {
          proportions[idx] <- NA
          next
        }
        m                <- nrow(adj_mat)
        max_edges        <- m * (m - 1) / 2
        n_edges          <- sum(adj_mat) / 2
        proportions[idx] <- n_edges / max_edges
      }
      
      result[[d_level]] <- data.frame(
        week       = y_c_weeks,
        proportion = proportions,
        mouse      = ID
      )
      
      rm(tmp_env, res_i, step_data)
    }
    
    return(result)
  }
  
  # --------------------------------------------------------------------------
  # Collect data across all mice for each discrete level
  # --------------------------------------------------------------------------
  
  # Get union of all IDs across all discrete levels
  all_IDs <- unique(unlist(lapply(discrete_levels, discover_IDs)))
  
  # Load proportions for every mouse
  all_data <- lapply(all_IDs, load_proportions)
  names(all_data) <- all_IDs
  
  # --------------------------------------------------------------------------
  # Build one plot per discrete level
  # --------------------------------------------------------------------------
  plot_list <- list()
  
  for (d_level in discrete_levels) {
    
    # Collect rows from all mice that have data for this stratum
    df_list <- lapply(all_IDs, function(ID) all_data[[ID]][[d_level]])
    df_list <- Filter(Negate(is.null), df_list)
    
    if (length(df_list) == 0) next
    
    plot_df        <- do.call(rbind, df_list)
    plot_df$mouse  <- factor(plot_df$mouse, levels = all_IDs)  # Tau first, then WT
    
    # Color palette: Tau in reds/oranges, WT in blues
    tau_mice <- all_IDs[grepl("^Tau", all_IDs)]
    wt_mice  <- all_IDs[grepl("^WT",  all_IDs)]
    tau_cols <- setNames(scales::hue_pal(h = c(0, 60))(length(tau_mice)),  tau_mice)
    wt_cols  <- setNames(scales::hue_pal(h = c(200, 260))(length(wt_mice)), wt_mice)
    color_map <- c(tau_cols, wt_cols)
    
    g <- ggplot(plot_df, aes(x = week, y = proportion, color = mouse, group = mouse)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      scale_color_manual(values = color_map) +
      scale_x_continuous(breaks = 17:38) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
      labs(
        title  = d_level,
        x      = "Week",
        y      = "Edge Proportion",
        color  = "Mouse"
      ) +
      theme_bw() +
      theme(
        axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      )
    
    plot_list[[d_level]] <- g
  }
  
  return(plot_list)
}

plot_edge_instability_all_mice <- function(results_folder, time_scale, discrete_levels) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: For each discrete stratum, plot edge instability over weeks for all
  #       available mice discovered automatically from results_folder.
  #       Instability at week t is defined as:
  #
  #         |E(t) △ E(t+1)| / max(|E(t)|, |E(t+1)|)
  #
  #       i.e. the fraction of edges that changed relative to the larger
  #       of the two edge sets. Value of 0 = identical graphs, 1 = fully
  #       disjoint. Only computed for consecutive present weeks (gap = 1).
  #       All computations use upper triangle only to exclude diagonal.
  #
  # inputs:
  #
  # - results_folder   (string)
  # - time_scale       (integer)  e.g. 10
  # - discrete_levels  (vector of strings)  e.g. c('m0vr0', 'm1vr1', ...)
  #
  # returns: named list of ggplot objects, one per discrete level
  #
  # ----------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Discover all IDs for a given discrete level
  # --------------------------------------------------------------------------
  discover_IDs <- function(d_level) {
    pattern   <- paste0("^(.+)_", d_level, "_t", time_scale, "\\.RData$")
    all_files <- list.files(results_folder, pattern = pattern, full.names = FALSE)
    ids       <- sub(paste0("_", d_level, "_t", time_scale, "\\.RData$"), "", all_files)
    tau_ids   <- sort(ids[grepl("^Tau", ids)])
    wt_ids    <- sort(ids[grepl("^WT",  ids)])
    other     <- sort(ids[!grepl("^(Tau|WT)", ids)])
    c(tau_ids, wt_ids, other)
  }
  
  # --------------------------------------------------------------------------
  # Load edge instability for one mouse across all discrete levels
  # --------------------------------------------------------------------------
  load_instability <- function(ID) {
    
    result <- list()
    
    for (d_level in discrete_levels) {
      
      file_path <- paste0(results_folder, '/', ID, '_', d_level, '_t', time_scale, '.RData')
      if (!file.exists(file_path)) next
      
      tmp_env <- new.env()
      load(file_path, envir = tmp_env)
      res_i <- tmp_env$graph_results_i
      
      step_data <- res_i$step_12b
      y_c_weeks <- as.numeric(res_i$y_c_query)
      
      # Build named list of adjacency matrices keyed by week
      adj_by_week <- list()
      for (idx in seq_along(y_c_weeks)) {
        adj_mat <- step_data[[idx]][["adj_mat_KL_GIC_local_est_eig3"]]
        if (!is.null(adj_mat)) {
          adj_by_week[[as.character(y_c_weeks[idx])]] <- adj_mat
        }
      }
      
      present_weeks     <- sort(as.numeric(names(adj_by_week)))
      instability_vals  <- c()
      instability_weeks <- c()
      
      for (w_idx in seq_len(length(present_weeks) - 1)) {
        w_curr <- present_weeks[w_idx]
        w_next <- present_weeks[w_idx + 1]
        
        # Only compute for consecutive weeks (gap of exactly 1)
        if (w_next - w_curr != 1) next
        
        A_curr <- adj_by_week[[as.character(w_curr)]]
        A_next <- adj_by_week[[as.character(w_next)]]
        
        # Align dimensions if neuron count differs across weeks
        n_min <- min(nrow(A_curr), nrow(A_next))
        if (n_min == 0) next
        A_curr <- A_curr[1:n_min, 1:n_min]
        A_next <- A_next[1:n_min, 1:n_min]
        
        # Use upper triangle only — excludes diagonal (self-loops) and
        # avoids double counting from symmetry
        upper       <- upper.tri(A_curr)
        edges_curr  <- sum(A_curr[upper])
        edges_next  <- sum(A_next[upper])
        union_edges <- sum(A_curr[upper] | A_next[upper])
        denom       <- union_edges
        
        if (denom == 0) {
          instability_vals <- c(instability_vals, 1)
          message(sprintf("    week %d -> %d: edges_curr=%g, edges_next=%g — both empty, instability=0",
                          w_curr, w_next, edges_curr, edges_next))
        } else {
          sym_diff        <- sum(abs(A_curr[upper] - A_next[upper]))
          instability_val <- 1 - sym_diff / denom
          instability_vals <- c(instability_vals, instability_val)
          message(sprintf("    week %d -> %d: edges_curr=%g, edges_next=%g, gained=%g, lost=%g, sym_diff=%g, denom=%g, instability=%.4f",
                          w_curr, w_next,
                          edges_curr, edges_next,
                          sum(A_next[upper] > A_curr[upper]),   # edges gained
                          sum(A_curr[upper] > A_next[upper]),   # edges lost
                          sym_diff, denom,
                          instability_val))
        }
        
        instability_weeks <- c(instability_weeks, w_curr)
      }
      
      if (length(instability_weeks) > 0) {
        result[[d_level]] <- data.frame(
          week        = instability_weeks,
          instability = instability_vals,
          mouse       = ID
        )
      }
      
      rm(tmp_env, res_i, step_data)
    }
    
    return(result)
  }
  
  # --------------------------------------------------------------------------
  # Collect data across all mice
  # --------------------------------------------------------------------------
  all_IDs  <- unique(unlist(lapply(discrete_levels, discover_IDs)))
  all_data <- lapply(all_IDs, load_instability)
  names(all_data) <- all_IDs
  
  # --------------------------------------------------------------------------
  # Color palette: Tau in warm tones, WT in cool tones
  # --------------------------------------------------------------------------
  tau_mice  <- all_IDs[grepl("^Tau", all_IDs)]
  wt_mice   <- all_IDs[grepl("^WT",  all_IDs)]
  tau_cols  <- setNames(scales::hue_pal(h = c(0, 60))(length(tau_mice)),   tau_mice)
  wt_cols   <- setNames(scales::hue_pal(h = c(200, 260))(length(wt_mice)), wt_mice)
  color_map <- c(tau_cols, wt_cols)
  
  # --------------------------------------------------------------------------
  # Build one plot per discrete level
  # --------------------------------------------------------------------------
  plot_list      <- list()
  plot_list_sqrt <- list()
  
  for (d_level in discrete_levels) {
    
    df_list <- lapply(all_IDs, function(ID) all_data[[ID]][[d_level]])
    df_list <- Filter(Negate(is.null), df_list)
    
    if (length(df_list) == 0) next
    
    plot_df       <- do.call(rbind, df_list)
    plot_df$mouse <- factor(plot_df$mouse, levels = all_IDs)
    
    base_plot <- ggplot(plot_df, aes(x = week, y = instability, color = mouse, group = mouse)) +
      geom_line(linewidth = 0.8) +
      scale_color_manual(values = color_map) +
      scale_x_continuous(breaks = c(20, 25, 30, 35)) +
      labs(
        x     = "Age (Weeks)",
        y     = "Jaccard Similarity",
        color = "Mouse"
      ) +
      guides(color = guide_legend(nrow = 1)) +
      theme_bw(base_size = 16) +
      theme(
        panel.grid      = element_blank(),
        legend.position = "bottom"
      )
    
    g <- base_plot +
      scale_y_continuous(
        breaks = c(0, 0.25, 0.5, 0.75, 1),
        limits = c(0, 1)
      )
    
    g_sqrt <- base_plot +
      scale_y_continuous(
        trans   = "sqrt",
        breaks  = c(0, 0.25, 0.5, 0.75, 1),
        limits  = c(0, 1)
      )
    
    plot_list[[d_level]]      <- g
    plot_list_sqrt[[d_level]] <- g_sqrt
  }
  
  return(list(linear = plot_list, sqrt = plot_list_sqrt))
}

# jaccard distance between edge sets of the same mouse, but different discrete strata
plot_strata_instability_all_mice <- function(results_folder, time_scale, discrete_levels) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Supply two discrete strata and plot Jaccard similarity over weeks for
  #       all available mice discovered automatically from results_folder.
  #       Jaccard similarity at week t is defined as:
  #
  #         |E_i(t) ∩ E_j(t)| / |E_i(t) ∪ E_j(t)|
  #
  #       where i and j are the two different strata (discrete_levels[1] and
  #       discrete_levels[2]). Value of 1 = identical graphs, 0 = fully disjoint.
  #       When both edge sets are empty at week t, similarity is defined as 1
  #       (both graphs agree on having no edges).
  #       All computations use upper triangle only to exclude diagonal.
  #
  # inputs:
  #
  # - results_folder   (string)
  # - time_scale       (integer)  e.g. 10
  # - discrete_levels  (vector of exactly two strings)  e.g. c('m0vr1', 'm1vr1')
  #
  # returns: named list with two elements:
  #   - linear : named list of ggplot objects (linear y-axis), one per mouse
  #   - sqrt   : named list of ggplot objects (sqrt y-axis), one per mouse
  #   Both sub-lists also contain a 'combined' entry with all mice on one plot.
  #
  # ----------------------------------------------------------------------------
  
  if (length(discrete_levels) != 2) {
    stop("plot_strata_instability_all_mice: discrete_levels must be exactly length 2")
  }
  
  d_level_A <- discrete_levels[1]
  d_level_B <- discrete_levels[2]
  
  # --------------------------------------------------------------------------
  # Discover all IDs that have files for BOTH discrete levels
  # --------------------------------------------------------------------------
  discover_IDs <- function() {
    pattern_A <- paste0("^(.+)_", d_level_A, "_t", time_scale, "\\.RData$")
    pattern_B <- paste0("^(.+)_", d_level_B, "_t", time_scale, "\\.RData$")
    
    files_A <- list.files(results_folder, pattern = pattern_A, full.names = FALSE)
    files_B <- list.files(results_folder, pattern = pattern_B, full.names = FALSE)
    
    ids_A <- sub(paste0("_", d_level_A, "_t", time_scale, "\\.RData$"), "", files_A)
    ids_B <- sub(paste0("_", d_level_B, "_t", time_scale, "\\.RData$"), "", files_B)
    
    # Only keep IDs present in BOTH strata
    ids <- intersect(ids_A, ids_B)
    
    tau_ids <- sort(ids[grepl("^Tau", ids)])
    wt_ids  <- sort(ids[grepl("^WT",  ids)])
    other   <- sort(ids[!grepl("^(Tau|WT)", ids)])
    c(tau_ids, wt_ids, other)
  }
  
  # --------------------------------------------------------------------------
  # Load adjacency matrices for one mouse and one discrete level,
  # keyed by week number
  # --------------------------------------------------------------------------
  load_adj_by_week <- function(ID, d_level) {
    file_path <- paste0(results_folder, '/', ID, '_', d_level, '_t', time_scale, '.RData')
    if (!file.exists(file_path)) return(NULL)
    
    tmp_env <- new.env()
    load(file_path, envir = tmp_env)
    res_i <- tmp_env$graph_results_i
    
    step_data <- res_i$step_12b
    y_c_weeks <- as.numeric(res_i$y_c_query)
    
    adj_by_week <- list()
    for (idx in seq_along(y_c_weeks)) {
      adj_mat <- step_data[[idx]][["adj_mat_KL_GIC_local_est_eig3"]]
      if (!is.null(adj_mat)) {
        adj_by_week[[as.character(y_c_weeks[idx])]] <- adj_mat
      }
    }
    
    rm(tmp_env, res_i, step_data)
    return(adj_by_week)
  }
  
  # --------------------------------------------------------------------------
  # Compute cross-stratum Jaccard similarity for one mouse across all weeks
  # --------------------------------------------------------------------------
  load_strata_similarity <- function(ID) {
    
    message(sprintf("  Loading mouse %s ...", ID))
    
    adj_A <- load_adj_by_week(ID, d_level_A)
    adj_B <- load_adj_by_week(ID, d_level_B)
    
    if (is.null(adj_A) || is.null(adj_B)) return(NULL)
    
    # Weeks present in both strata
    weeks_A       <- as.numeric(names(adj_A))
    weeks_B       <- as.numeric(names(adj_B))
    common_weeks  <- sort(intersect(weeks_A, weeks_B))
    
    if (length(common_weeks) == 0) {
      message(sprintf("    No common weeks found for mouse %s — skipping.", ID))
      return(NULL)
    }
    
    similarity_vals  <- numeric(length(common_weeks))
    
    for (w_idx in seq_along(common_weeks)) {
      w      <- common_weeks[w_idx]
      A_curr <- adj_A[[as.character(w)]]
      B_curr <- adj_B[[as.character(w)]]
      
      # Align dimensions if neuron count differs between strata
      n_min <- min(nrow(A_curr), nrow(B_curr))
      if (n_min == 0) {
        similarity_vals[w_idx] <- NA
        next
      }
      A_curr <- A_curr[1:n_min, 1:n_min]
      B_curr <- B_curr[1:n_min, 1:n_min]
      
      # Use upper triangle only — excludes diagonal (self-loops) and
      # avoids double counting from symmetry
      upper        <- upper.tri(A_curr)
      edges_A      <- sum(A_curr[upper])
      edges_B      <- sum(B_curr[upper])
      intersection <- sum(A_curr[upper] & B_curr[upper])
      union_edges  <- sum(A_curr[upper] | B_curr[upper])
      
      if (union_edges == 0) {
        # Both graphs have no edges — define similarity as 1 (they agree)
        jaccard_val <- 1
        message(sprintf("    week %d: edges_A=%g, edges_B=%g — both empty, similarity=1",
                        w, edges_A, edges_B))
      } else {
        jaccard_val <- intersection / union_edges
        message(sprintf("    week %d: edges_A=%g, edges_B=%g, intersection=%g, union=%g, similarity=%.4f",
                        w, edges_A, edges_B, intersection, union_edges, jaccard_val))
      }
      
      similarity_vals[w_idx] <- jaccard_val
    }
    
    data.frame(
      week       = common_weeks,
      similarity = similarity_vals,
      mouse      = ID
    )
  }
  
  # --------------------------------------------------------------------------
  # Collect data across all mice
  # --------------------------------------------------------------------------
  all_IDs  <- discover_IDs()
  
  if (length(all_IDs) == 0) {
    message(sprintf("No mice found with files for both %s and %s.", d_level_A, d_level_B))
    return(list(linear = list(), sqrt = list()))
  }
  
  message(sprintf("Computing cross-stratum Jaccard similarity: %s vs %s", d_level_A, d_level_B))
  message(sprintf("Mice found: %s", paste(all_IDs, collapse = ", ")))
  
  all_data <- lapply(all_IDs, load_strata_similarity)
  names(all_data) <- all_IDs
  
  # Remove mice with no data
  all_data <- Filter(Negate(is.null), all_data)
  present_IDs <- names(all_data)
  
  if (length(present_IDs) == 0) {
    message("No valid data found for any mouse.")
    return(list(linear = list(), sqrt = list()))
  }
  
  # --------------------------------------------------------------------------
  # Color palette: Tau in warm tones, WT in cool tones
  # --------------------------------------------------------------------------
  tau_mice  <- present_IDs[grepl("^Tau", present_IDs)]
  wt_mice   <- present_IDs[grepl("^WT",  present_IDs)]
  tau_cols  <- setNames(scales::hue_pal(h = c(0, 60))(max(length(tau_mice), 1)),   tau_mice)
  wt_cols   <- setNames(scales::hue_pal(h = c(200, 260))(max(length(wt_mice), 1)), wt_mice)
  color_map <- c(tau_cols, wt_cols)
  
  # --------------------------------------------------------------------------
  # Build combined data frame for the all-mice overlay plot
  # --------------------------------------------------------------------------
  combined_df        <- do.call(rbind, all_data)
  combined_df$mouse  <- factor(combined_df$mouse, levels = present_IDs)
  
  plot_title <- sprintf("%s vs %s", d_level_A, d_level_B)
  
  # --------------------------------------------------------------------------
  # Helper: build base ggplot for a given data frame
  # --------------------------------------------------------------------------
  make_base_plot <- function(df) {
    ggplot(df, aes(x = week, y = similarity, color = mouse, group = mouse)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      scale_color_manual(values = color_map) +
      scale_x_continuous(breaks = c(20, 25, 30, 35)) +
      labs(
        title = plot_title,
        x     = "Age (Weeks)",
        y     = "Jaccard Similarity",
        color = "Mouse"
      ) +
      guides(color = guide_legend(nrow = 1)) +
      theme_bw(base_size = 16) +
      theme(
        panel.grid      = element_blank(),
        legend.position = "bottom"
      )
  }
  
  
  
  # --------------------------------------------------------------------------
  # Build combined (all-mice overlay) plots
  # --------------------------------------------------------------------------
  base_combined <- make_base_plot(combined_df)
  
  plot_list_linear <- list()
  plot_list_sqrt   <- list()
  
  plot_list_linear[["combined"]] <- base_combined +
    scale_y_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      limits = c(0, 1)
    )
  
  plot_list_sqrt[["combined"]] <- base_combined +
    scale_y_continuous(
      trans  = "sqrt",
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      limits = c(0, 1)
    )
  
  return(list(linear = plot_list_linear, sqrt = plot_list_sqrt))
}

# plot % of neurons that connect HIP-HIP, HIP-EHC, or EHC-EHC
plot_edge_regional_proportion_all_mice <- function(results_folder, time_scale, discrete_levels) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: For each discrete stratum, plot the proportion of existing edges that
  #       are HIP-HIP, HIP-EHC, or EHC-EHC across all discovered mice.
  #
  #       n_hip inferred per mouse from res_i$recovery_params$kept_neuron_regions
  #       (same source used by visualize_strata_all_mice_v2).
  #       All computations use upper.tri() only.
  #
  # inputs:
  #   results_folder   (string)
  #   time_scale       (integer)  e.g. 10
  #   discrete_levels  (vector of strings)  e.g. c('m0vr0', 'm1vr1', ...)
  #
  # returns: named list of ggplot objects, keyed as e.g. "m0vr1_HIPHIP"
  #
  # ----------------------------------------------------------------------------
  
  discover_IDs <- function(d_level) {
    pattern   <- paste0("^(.+)_", d_level, "_t", time_scale, "\\.RData$")
    all_files <- list.files(results_folder, pattern = pattern, full.names = FALSE)
    ids       <- sub(paste0("_", d_level, "_t", time_scale, "\\.RData$"), "", all_files)
    tau_ids   <- sort(ids[grepl("^Tau", ids)])
    wt_ids    <- sort(ids[grepl("^WT",  ids)])
    other     <- sort(ids[!grepl("^(Tau|WT)", ids)])
    c(tau_ids, wt_ids, other)
  }
  
  region_pairs <- c("HIP-HIP", "HIP-EHC", "EHC-EHC")
  
  all_IDs <- unique(unlist(lapply(discrete_levels, discover_IDs)))
  
  # Color palette: Tau in warm tones, WT in cool tones
  tau_mice  <- all_IDs[grepl("^Tau", all_IDs)]
  wt_mice   <- all_IDs[grepl("^WT",  all_IDs)]
  tau_cols  <- setNames(scales::hue_pal(h = c(0, 60))(max(length(tau_mice), 1)),   tau_mice)
  wt_cols   <- setNames(scales::hue_pal(h = c(200, 260))(max(length(wt_mice), 1)), wt_mice)
  color_map <- c(tau_cols, wt_cols)
  
  plot_list <- list()
  
  for (lvl in discrete_levels) {
    
    ids <- discover_IDs(lvl)
    message(sprintf("[plot_edge_regional_proportion] level=%s, found %d IDs: %s",
                    lvl, length(ids), paste(ids, collapse = ", ")))
    
    all_data <- data.frame()
    
    for (ID in ids) {
      
      file_path <- paste0(results_folder, '/', ID, '_', lvl, '_t', time_scale, '.RData')
      if (!file.exists(file_path)) {
        message(sprintf("  [SKIP] not found: %s", file_path))
        next
      }
      
      tmp_env <- new.env()
      load(file_path, envir = tmp_env)
      res_i <- tmp_env$graph_results_i
      
      step_data <- res_i$step_12b
      y_c_weeks <- as.numeric(res_i$y_c_query)
      
      # Infer n_hip from the same source as visualize_strata_all_mice_v2
      regions <- res_i$recovery_params$kept_neuron_regions
      n_hip   <- sum(regions == regions[1])  # count of first region label
      p       <- length(regions)
      n_ehc   <- p - n_hip
      message(sprintf("  [%s] p=%d, n_hip=%d, n_ehc=%d, weeks=%d",
                      ID, p, n_hip, n_ehc, length(y_c_weeks)))
      
      for (idx in seq_along(y_c_weeks)) {
        
        week    <- y_c_weeks[idx]
        adj_mat <- step_data[[idx]][["adj_mat_KL_GIC_local_est_eig3"]]
        
        if (is.null(adj_mat)) {
          all_data <- rbind(all_data, data.frame(
            week = week, HIP_HIP = NA_real_, HIP_EHC = NA_real_, EHC_EHC = NA_real_, mouse = ID
          ))
          next
        }
        
        ut <- which(upper.tri(adj_mat) & adj_mat != 0, arr.ind = TRUE)
        
        if (nrow(ut) == 0) {
          all_data <- rbind(all_data, data.frame(
            week = week, HIP_HIP = NA_real_, HIP_EHC = NA_real_, EHC_EHC = NA_real_, mouse = ID
          ))
        } else {
          is_hip  <- function(i) i <= n_hip
          hip_hip <- sum( is_hip(ut[,1]) &  is_hip(ut[,2]))
          hip_ehc <- sum( is_hip(ut[,1]) & !is_hip(ut[,2])) +
            sum(!is_hip(ut[,1]) &  is_hip(ut[,2]))
          ehc_ehc <- sum(!is_hip(ut[,1]) & !is_hip(ut[,2]))
          total   <- hip_hip + hip_ehc + ehc_ehc
          
          all_data <- rbind(all_data, data.frame(
            week    = week,
            HIP_HIP = hip_hip / total,
            HIP_EHC = hip_ehc / total,
            EHC_EHC = ehc_ehc / total,
            mouse   = ID
          ))
        }
      }
      
      rm(tmp_env, res_i, step_data)
    }
    
    if (nrow(all_data) == 0) {
      message(sprintf("  [WARN] no data for level=%s", lvl))
      next
    }
    
    all_data$mouse <- factor(all_data$mouse, levels = all_IDs)
    
    for (rp in region_pairs) {
      
      col_name <- switch(rp,
                         "HIP-HIP" = "HIP_HIP",
                         "HIP-EHC" = "HIP_EHC",
                         "EHC-EHC" = "EHC_EHC")
      
      p_plot <- ggplot(all_data, aes(x = week, y = .data[[col_name]],
                                     color = mouse, group = mouse)) +
        geom_line(linewidth = 0.8) +
        geom_point(size = 2) +
        scale_color_manual(values = color_map) +
        scale_x_continuous(breaks = c(20, 25, 30, 35)) +
        scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
        labs(
          title = paste(lvl, "|", rp),
          x     = "Age (Weeks)",
          y     = paste("Proportion of", rp, "edges"),
          color = "Mouse"
        ) +
        guides(color = guide_legend(nrow = 1)) +
        theme_bw(base_size = 16) +
        theme(
          panel.grid      = element_blank(),
          legend.position = "bottom"
        )
      
      plot_key <- paste0(lvl, "_", gsub("-", "", rp))
      plot_list[[plot_key]] <- p_plot
      message(sprintf("  [OK] plot created: %s", plot_key))
    }
  }
  
  return(plot_list)
}

plot_median_nonzero_degree_all_mice <- function(results_folder, time_scale, discrete_levels) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: For each discrete stratum, and for each mouse, compute the 25th
  #       percentile, median, and 75th percentile of degree across neurons
  #       that have at least one edge (positive degree) at each week.
  #       Plots one loess curve per mouse (fit to the median), with a shaded
  #       band spanning the 25th-75th percentile loess fits. The 50% CI band
  #       reflects the IQR of the non-zero degree distribution at each week.
  #
  #       Returns one ggplot per discrete level.
  #
  # inputs:
  #   results_folder   (string)
  #   time_scale       (integer)  e.g. 10
  #   discrete_levels  (vector of strings)  e.g. c('m0vr0', 'm1vr1', ...)
  #
  # returns: named list of ggplot objects, one per discrete level
  #
  # ----------------------------------------------------------------------------
  
  discover_IDs <- function(d_level) {
    pattern   <- paste0("^(.+)_", d_level, "_t", time_scale, "\\.RData$")
    all_files <- list.files(results_folder, pattern = pattern, full.names = FALSE)
    ids       <- sub(paste0("_", d_level, "_t", time_scale, "\\.RData$"), "", all_files)
    tau_ids   <- sort(ids[grepl("^Tau", ids)])
    wt_ids    <- sort(ids[grepl("^WT",  ids)])
    other     <- sort(ids[!grepl("^(Tau|WT)", ids)])
    c(tau_ids, wt_ids, other)
  }
  
  all_IDs  <- unique(unlist(lapply(discrete_levels, discover_IDs)))
  tau_mice <- all_IDs[grepl("^Tau", all_IDs)]
  wt_mice  <- all_IDs[grepl("^WT",  all_IDs)]
  tau_cols <- setNames(scales::hue_pal(h = c(0, 60))(max(length(tau_mice), 1)),   tau_mice)
  wt_cols  <- setNames(scales::hue_pal(h = c(200, 260))(max(length(wt_mice), 1)), wt_mice)
  color_map <- c(tau_cols, wt_cols)
  
  plot_list <- list()
  
  for (lvl in discrete_levels) {
    
    ids <- discover_IDs(lvl)
    message(sprintf("[plot_median_nonzero_degree] level=%s, found %d IDs: %s",
                    lvl, length(ids), paste(ids, collapse = ", ")))
    
    all_data <- data.frame()
    
    for (ID in ids) {
      
      file_path <- paste0(results_folder, '/', ID, '_', lvl, '_t', time_scale, '.RData')
      if (!file.exists(file_path)) {
        message(sprintf("  [SKIP] not found: %s", file_path))
        next
      }
      
      tmp_env <- new.env()
      load(file_path, envir = tmp_env)
      res_i <- tmp_env$graph_results_i
      
      step_data <- res_i$step_12b
      y_c_weeks <- as.numeric(res_i$y_c_query)
      
      message(sprintf("  [%s] weeks=%d", ID, length(y_c_weeks)))
      
      for (idx in seq_along(y_c_weeks)) {
        
        week    <- y_c_weeks[idx]
        adj_mat <- step_data[[idx]][["adj_mat_KL_GIC_local_est_eig3"]]
        
        if (is.null(adj_mat)) next
        
        # Degree = row sums of upper + lower triangle (symmetric), excluding diagonal
        diag(adj_mat) <- 0
        degree <- rowSums(adj_mat) / (nrow(adj_mat) - 1)
        
        # Keep only neurons with positive degree
        nonzero_degree <- degree[degree > 0]
        
        if (length(nonzero_degree) == 0) next
        
        all_data <- rbind(all_data, data.frame(
          week   = week,
          q25    = quantile(nonzero_degree, 0.25),
          median = quantile(nonzero_degree, 0.50),
          q75    = quantile(nonzero_degree, 0.75),
          mouse  = ID
        ))
      }
      
      rm(tmp_env, res_i, step_data)
    }
    
    if (nrow(all_data) == 0) {
      message(sprintf("  [WARN] no data for level=%s", lvl))
      next
    }
    
    all_data$mouse <- factor(all_data$mouse, levels = all_IDs)
    
    p_plot <- ggplot(all_data, aes(x = week, color = mouse, fill = mouse, group = mouse)) +
      # Shaded IQR band via loess on q25 and q75
      geom_ribbon(
        aes(ymin = q25, ymax = q75),
        stat = "smooth", method = "loess", formula = y ~ x,
        alpha = 0.15, color = NA
      ) +
      # Median loess line
      geom_smooth(
        aes(y = median),
        method = "loess", formula = y ~ x,
        se = FALSE, linewidth = 0.9
      ) +
      scale_color_manual(values = color_map) +
      scale_fill_manual(values = color_map) +
      scale_x_continuous(breaks = c(20, 25, 30, 35)) +
      scale_y_continuous(limits = c(0, NA)) +
      labs(
        title = lvl,
        x     = "Age (Weeks)",
        y     = "Normalized Degree (non-zero neurons)",
        color = "Mouse",
        fill  = "Mouse"
      ) +
      guides(color = guide_legend(nrow = 1),
             fill  = guide_legend(nrow = 1)) +
      theme_bw(base_size = 16) +
      theme(
        panel.grid      = element_blank(),
        legend.position = "bottom"
      )
    
    plot_list[[lvl]] <- p_plot
    message(sprintf("  [OK] plot created: %s", lvl))
  }
  
  return(plot_list)
}