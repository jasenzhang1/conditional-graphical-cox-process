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
      geom_line(size = 0.8) +
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

# temporal jaccard distance - dots and lines
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
  # returns: named list with two ggplot objects (linear, sqrt), each a
  #          single faceted plot with all discrete strata stacked vertically,
  #          one raw geom_line per mouse
  #
  # ----------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Stratum display labels
  # --------------------------------------------------------------------------
  stratum_labels <- c(m0vr1 = "Resting", m1vr1 = "Running")
  
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
  # Combine all strata into a single long data frame, applying display labels
  # where available (e.g. m0vr1 -> "Resting", m1vr1 -> "Running").
  # Strata without a label entry keep their original name.
  # The stratum column is an ordered factor to control facet ordering
  # (Resting on top, Running below).
  # --------------------------------------------------------------------------
  all_rows <- list()
  for (d_level in discrete_levels) {
    df_list <- lapply(all_IDs, function(ID) all_data[[ID]][[d_level]])
    df_list <- Filter(Negate(is.null), df_list)
    if (length(df_list) == 0) next
    df           <- do.call(rbind, df_list)
    df$stratum   <- ifelse(d_level %in% names(stratum_labels),
                           stratum_labels[[d_level]],
                           d_level)
    all_rows[[d_level]] <- df
  }
  
  plot_df       <- do.call(rbind, all_rows)
  plot_df$mouse <- factor(plot_df$mouse, levels = all_IDs)
  
  # Derive ordered factor levels from discrete_levels order, mapped through labels
  stratum_level_order <- ifelse(discrete_levels %in% names(stratum_labels),
                                stratum_labels[discrete_levels],
                                discrete_levels)
  plot_df$stratum <- factor(plot_df$stratum, levels = stratum_level_order)
  
  # --------------------------------------------------------------------------
  # Build the shared base plot: one raw line per mouse (geom_line + geom_point),
  # faceted by stratum and stacked vertically (ncol = 1).
  # Legend placed below all panels.
  # --------------------------------------------------------------------------
  base_plot <- ggplot(plot_df, aes(x = week, y = instability, color = mouse, group = mouse)) +
    geom_line(size = 0.8) +
    geom_point(size = 1.5) +
    facet_wrap(~ stratum, ncol = 1) +
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
  
  # Linear y-axis: range 0–0.5 with ticks only at 0.00, 0.25, 0.50
  g <- base_plot +
    scale_y_continuous(
      breaks = c(0, 0.25, 0.50),
      limits = c(0, 0.50)
    )
  
  # Sqrt-transformed y-axis: same tick marks, same range
  g_sqrt <- base_plot +
    scale_y_continuous(
      trans   = "sqrt",
      breaks  = c(0, 0.25, 0.50),
      limits  = c(0, 0.50)
    )
  
  return(list(linear = g, sqrt = g_sqrt))
}

# temporal jaccard distance - loess
plot_edge_instability_all_mice_v2 <- function(results_folder, time_scale, discrete_levels) {
  
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
  #       Loess is pre-computed separately per mouse x stratum to guarantee
  #       within-panel fitting (geom_smooth pools across facets).
  #
  # inputs:
  #
  # - results_folder   (string)
  # - time_scale       (integer)  e.g. 10
  # - discrete_levels  (vector of strings)  e.g. c('m0vr0', 'm1vr1', ...)
  #
  # returns: named list with two ggplot objects (linear, sqrt), each a
  #          single faceted plot with all discrete strata stacked vertically
  #
  # ----------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Stratum display labels
  # --------------------------------------------------------------------------
  stratum_labels <- c(m0vr1 = "Resting", m1vr1 = "Running")
  
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
      
      stratum_label <- ifelse(d_level %in% names(stratum_labels),
                              stratum_labels[[d_level]],
                              d_level)
      
      file_path <- paste0(results_folder, '/', ID, '_', d_level, '_t', time_scale, '.RData')
      if (!file.exists(file_path)) {
        message(sprintf("  [SKIP] mouse=%s, stratum=%s: file not found", ID, stratum_label))
        next
      }
      
      message(sprintf("  [LOAD] mouse=%s, stratum=%s", ID, stratum_label))
      
      tmp_env <- new.env()
      load(file_path, envir = tmp_env)
      res_i <- tmp_env$graph_results_i
      
      step_data <- res_i$step_12b
      y_c_weeks <- as.numeric(res_i$y_c_query)
      
      # Build named list of adjacency matrices keyed by week;
      # keys stored as character but sorted numerically below
      adj_by_week <- list()
      for (idx in seq_along(y_c_weeks)) {
        adj_mat <- step_data[[idx]][["adj_mat_KL_GIC_local_est_eig3"]]
        if (!is.null(adj_mat)) {
          adj_by_week[[as.character(y_c_weeks[idx])]] <- adj_mat
        }
      }
      
      # Numeric sort is essential — string sort of week keys can misordering
      # (e.g. "9" > "38" lexicographically), creating spurious consecutive pairs
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
        
        # Guard against NULL lookup (e.g. week stored with unexpected key format)
        if (is.null(A_curr) || is.null(A_next)) {
          message(sprintf("    [%s | %s] week %d or %d: NULL adj_mat — skipping pair",
                          ID, stratum_label, w_curr, w_next))
          next
        }
        
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
          # Both graphs have no edges — instability is 0 (they agree)
          instability_vals <- c(instability_vals, 0)
          message(sprintf("    [%s | %s] week %d -> %d: edges_curr=%g, edges_next=%g — both empty, instability=0",
                          ID, stratum_label, w_curr, w_next, edges_curr, edges_next))
        } else {
          sym_diff        <- sum(abs(A_curr[upper] - A_next[upper]))
          instability_val <- sym_diff / denom
          instability_vals <- c(instability_vals, instability_val)
          message(sprintf("    [%s | %s] week %d -> %d: edges_curr=%g, edges_next=%g, gained=%g, lost=%g, sym_diff=%g, denom=%g, instability=%.4f",
                          ID, stratum_label, w_curr, w_next,
                          edges_curr, edges_next,
                          sum(A_next[upper] > A_curr[upper]),
                          sum(A_curr[upper] > A_next[upper]),
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
      } else {
        message(sprintf("  [WARN] mouse=%s, stratum=%s: no consecutive week pairs found",
                        ID, stratum_label))
      }
      
      rm(tmp_env, res_i, step_data)
    }
    
    return(result)
  }
  
  # --------------------------------------------------------------------------
  # Collect data across all mice
  # --------------------------------------------------------------------------
  all_IDs  <- unique(unlist(lapply(discrete_levels, discover_IDs)))
  
  message(sprintf("[plot_edge_instability] mice found: %s", paste(all_IDs, collapse = ", ")))
  
  all_data        <- lapply(all_IDs, load_instability)
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
  # Combine all strata into a single long data frame, applying display labels
  # where available (e.g. m0vr1 -> "Resting", m1vr1 -> "Running").
  # Strata without a label entry keep their original name.
  # Panel labels assigned as plain character before rbind, then converted to
  # ordered factor after, to prevent level contamination across blocks.
  # --------------------------------------------------------------------------
  all_rows <- list()
  for (d_level in discrete_levels) {
    df_list <- lapply(all_IDs, function(ID) all_data[[ID]][[d_level]])
    df_list <- Filter(Negate(is.null), df_list)
    if (length(df_list) == 0) next
    df           <- do.call(rbind, df_list)
    df$stratum   <- ifelse(d_level %in% names(stratum_labels),
                           stratum_labels[[d_level]],
                           d_level)
    all_rows[[d_level]] <- df
  }
  
  plot_df        <- do.call(rbind, all_rows)
  plot_df$mouse  <- factor(plot_df$mouse, levels = all_IDs)
  
  # Derive ordered factor levels from discrete_levels order, mapped through labels
  stratum_level_order <- ifelse(discrete_levels %in% names(stratum_labels),
                                stratum_labels[discrete_levels],
                                discrete_levels)
  plot_df$stratum <- factor(plot_df$stratum, levels = stratum_level_order)
  
  # Sanity check: row counts per mouse x stratum
  message("[plot_edge_instability] row counts per mouse x stratum:")
  print(table(plot_df$mouse, plot_df$stratum, useNA = "always"))
  
  # --------------------------------------------------------------------------
  # Pre-compute loess fits separately per mouse x stratum.
  # Necessary because geom_smooth pools all data before faceting,
  # fitting loess across all strata rather than within each.
  # --------------------------------------------------------------------------
  loess_lines <- do.call(rbind, lapply(
    split(plot_df, list(plot_df$mouse, plot_df$stratum), drop = TRUE),
    function(df) {
      df <- df[!is.na(df$instability), ]
      if (nrow(df) < 4) return(NULL)
      week_seq <- seq(min(df$week), max(df$week), length.out = 200)
      fit <- tryCatch(
        predict(loess(instability ~ week, data = df), newdata = data.frame(week = week_seq)),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      data.frame(
        week        = week_seq,
        instability = fit,
        mouse       = df$mouse[1],
        stratum     = df$stratum[1]
      )
    }
  ))
  loess_lines$mouse   <- factor(loess_lines$mouse,   levels = all_IDs)
  loess_lines$stratum <- factor(loess_lines$stratum, levels = stratum_level_order)
  
  message("[plot_edge_instability] loess row counts per mouse x stratum:")
  print(table(loess_lines$mouse, loess_lines$stratum, useNA = "always"))
  
  # --------------------------------------------------------------------------
  # Build the shared base plot: loess curves drawn from pre-computed
  # predictions via geom_line, faceted by stratum and stacked vertically
  # (ncol = 1). Legend placed below all panels.
  # --------------------------------------------------------------------------
  base_plot <- ggplot(loess_lines,
                      aes(x = week, y = instability, color = mouse, group = mouse)) +
    geom_line(size = 0.9) +
    facet_wrap(~ stratum, ncol = 1) +
    scale_color_manual(values = color_map) +
    scale_x_continuous(breaks = c(20, 25, 30, 35)) +
    labs(
      x     = "Age (Weeks)",
      y     = "Edge Instability",
      color = "Mouse"
    ) +
    guides(color = guide_legend(nrow = 1)) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid      = element_blank(),
      legend.position = "bottom"
    )
  
  # Linear y-axis: range 0–0.5 with ticks only at 0.00, 0.25, 0.50
  g <- base_plot +
    scale_y_continuous(
      breaks = c(0, 0.2, 0.4, 0.6),
      limits = c(0, 0.7)
    )
  
  # Sqrt-transformed y-axis: same tick marks, same range
  g_sqrt <- base_plot +
    scale_y_continuous(
      trans  = "sqrt",
      breaks = c(0, 0.2, 0.4, 0.6),
      limits = c(0, 0.7)
    )
  
  return(list(linear = g, sqrt = g_sqrt))
}

# strata jaccard distance - dots and lines
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

# strata jaccard distance - loess
plot_strata_instability_all_mice_v2 <- function(results_folder, time_scale, discrete_levels) {
  
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
  #       Loess is pre-computed per mouse to avoid geom_smooth pooling across
  #       groups.
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
  # keyed by week number (numeric, stored as character key)
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
    
    # Weeks present in both strata — sorted numerically
    weeks_A      <- as.numeric(names(adj_A))
    weeks_B      <- as.numeric(names(adj_B))
    common_weeks <- sort(intersect(weeks_A, weeks_B))
    
    if (length(common_weeks) == 0) {
      message(sprintf("    No common weeks found for mouse %s — skipping.", ID))
      return(NULL)
    }
    
    similarity_vals <- numeric(length(common_weeks))
    
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
  all_IDs <- discover_IDs()
  
  if (length(all_IDs) == 0) {
    message(sprintf("No mice found with files for both %s and %s.", d_level_A, d_level_B))
    return(list(linear = list(), sqrt = list()))
  }
  
  message(sprintf("Computing cross-stratum Jaccard similarity: %s vs %s", d_level_A, d_level_B))
  message(sprintf("Mice found: %s", paste(all_IDs, collapse = ", ")))
  
  all_data <- lapply(all_IDs, load_strata_similarity)
  names(all_data) <- all_IDs
  
  # Remove mice with no data
  all_data    <- Filter(Negate(is.null), all_data)
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
  combined_df       <- do.call(rbind, all_data)
  combined_df$mouse <- factor(combined_df$mouse, levels = present_IDs)
  
  plot_title <- sprintf("%s vs %s", d_level_A, d_level_B)
  
  # --------------------------------------------------------------------------
  # Pre-compute loess fits separately per mouse.
  # Avoids geom_smooth pooling across groups; also ensures NA rows (weeks
  # where both strata had no edges) are dropped cleanly before fitting.
  # --------------------------------------------------------------------------
  loess_lines <- do.call(rbind, lapply(
    split(combined_df, combined_df$mouse, drop = TRUE),
    function(df) {
      df <- df[!is.na(df$similarity), ]
      if (nrow(df) < 4) return(NULL)
      week_seq <- seq(min(df$week), max(df$week), length.out = 200)
      fit <- tryCatch(
        predict(loess(similarity ~ week, data = df), newdata = data.frame(week = week_seq)),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      data.frame(
        week       = week_seq,
        similarity = fit,
        mouse      = df$mouse[1]
      )
    }
  ))
  loess_lines$mouse <- factor(loess_lines$mouse, levels = present_IDs)
  
  message(sprintf("[plot_strata_instability] loess computed for mice: %s",
                  paste(unique(as.character(loess_lines$mouse)), collapse = ", ")))
  
  # --------------------------------------------------------------------------
  # Helper: build base ggplot for a given loess data frame
  # --------------------------------------------------------------------------
  make_base_plot <- function(df) {
    ggplot(df, aes(x = week, y = similarity, color = mouse, group = mouse)) +
      geom_line(size = 0.9) +
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
  base_combined <- make_base_plot(loess_lines)
  
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

# jaccard distance both temporally and cross-sectionally!!
# Jaccard similarity across consecutive weeks (within stratum) and
# between strata at the same week, all mice, all in one faceted plot
plot_edge_stability_combined <- function(results_folder, time_scale, discrete_levels) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Produce a single faceted plot with three panels stacked vertically:
  #
  #   Panel 1 ("Temporal, Resting"):
  #     Jaccard similarity between consecutive-week edge sets within m0vr1.
  #     J(t) = |E(t) ∩ E(t+1)| / |E(t) ∪ E(t+1)|
  #     Only computed for weeks with gap = 1.
  #
  #   Panel 2 ("Temporal, Running"):
  #     Same as Panel 1 but within m1vr1.
  #
  #   Panel 3 ("Cross-Stratum, Same Week"):
  #     Jaccard similarity between Resting and Running edge sets at
  #     the same week t.
  #     J(t) = |E_rest(t) ∩ E_run(t)| / |E_rest(t) ∪ E_run(t)|
  #     Only computed for weeks present in both strata.
  #
  #   All computations use upper triangle only to exclude diagonal
  #   (self-loops) and avoid double counting from symmetry.
  #   When both edge sets are empty, similarity is defined as 1
  #   (both graphs agree on having no edges).
  #   Loess is pre-computed separately per mouse x panel to guarantee
  #   within-panel fitting (geom_smooth pools across facets).
  #   Panel labels are assigned after rbind to prevent factor-level
  #   contamination across row blocks.
  #
  # inputs:
  #
  # - results_folder   (string)
  # - time_scale       (integer)  e.g. 10
  # - discrete_levels  (vector of exactly two strings) c('m0vr1', 'm1vr1')
  #
  # returns: named list with two ggplot objects:
  #   - linear : faceted plot with linear y-axis
  #   - sqrt   : faceted plot with sqrt-transformed y-axis
  #
  # ----------------------------------------------------------------------------
  
  if (length(discrete_levels) != 2) {
    stop("plot_edge_stability_combined: discrete_levels must be exactly length 2")
  }
  
  d_level_A <- discrete_levels[1]   # Resting  (m0vr1)
  d_level_B <- discrete_levels[2]   # Running  (m1vr1)
  
  # --------------------------------------------------------------------------
  # Panel display labels
  # --------------------------------------------------------------------------
  panel_labels <- c(
    A_temporal = "Temporal, Resting",
    B_temporal = "Temporal, Running",
    cross      = "Cross-Stratum"
  )
  panel_order <- unname(panel_labels)   # controls facet stacking order
  
  # --------------------------------------------------------------------------
  # Discover all IDs present in results_folder for a given discrete level
  # --------------------------------------------------------------------------
  discover_IDs_for_level <- function(d_level) {
    pattern   <- paste0("^(.+)_", d_level, "_t", time_scale, "\\.RData$")
    all_files <- list.files(results_folder, pattern = pattern, full.names = FALSE)
    ids       <- sub(paste0("_", d_level, "_t", time_scale, "\\.RData$"), "", all_files)
    tau_ids   <- sort(ids[grepl("^Tau", ids)])
    wt_ids    <- sort(ids[grepl("^WT",  ids)])
    other     <- sort(ids[!grepl("^(Tau|WT)", ids)])
    c(tau_ids, wt_ids, other)
  }
  
  all_IDs_A <- discover_IDs_for_level(d_level_A)
  all_IDs_B <- discover_IDs_for_level(d_level_B)
  all_IDs   <- unique(c(all_IDs_A, all_IDs_B))
  
  # --------------------------------------------------------------------------
  # Load adjacency matrices for one mouse and one discrete level,
  # keyed by week number (numeric, stored as character key)
  # --------------------------------------------------------------------------
  load_adj_by_week <- function(ID, d_level) {
    file_path <- paste0(results_folder, '/', ID, '_', d_level, '_t', time_scale, '.RData')
    if (!file.exists(file_path)) return(NULL)
    
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
    
    rm(tmp_env, res_i, step_data)
    return(adj_by_week)
  }
  
  # --------------------------------------------------------------------------
  # Compute temporal Jaccard similarity (consecutive weeks, gap == 1) for
  # one mouse within a single discrete level. Returns a data.frame or NULL.
  # Week index assigned at w_curr (the earlier of the two weeks).
  # present_weeks is sorted numerically to prevent lexicographic misordering.
  # --------------------------------------------------------------------------
  compute_temporal_jaccard <- function(ID, d_level) {
    
    adj_by_week <- load_adj_by_week(ID, d_level)
    if (is.null(adj_by_week) || length(adj_by_week) < 2) return(NULL)
    
    # Numeric sort is essential — string sort of week keys can misordering
    # (e.g. "9" > "38" lexicographically), creating spurious consecutive pairs
    present_weeks    <- sort(as.numeric(names(adj_by_week)))
    similarity_vals  <- c()
    similarity_weeks <- c()
    
    for (w_idx in seq_len(length(present_weeks) - 1)) {
      w_curr <- present_weeks[w_idx]
      w_next <- present_weeks[w_idx + 1]
      
      # Only compute for strictly consecutive weeks (gap of exactly 1)
      if (w_next - w_curr != 1) next
      
      A_curr <- adj_by_week[[as.character(w_curr)]]
      A_next <- adj_by_week[[as.character(w_next)]]
      
      # Guard against NULL lookup (e.g. week stored with unexpected key format)
      if (is.null(A_curr) || is.null(A_next)) {
        message(sprintf("    [%s %s] week %d or %d: NULL adj_mat — skipping pair",
                        ID, d_level, w_curr, w_next))
        next
      }
      
      # Align dimensions if neuron count differs across weeks
      n_min <- min(nrow(A_curr), nrow(A_next))
      if (n_min == 0) next
      A_curr <- A_curr[1:n_min, 1:n_min]
      A_next <- A_next[1:n_min, 1:n_min]
      
      # Use upper triangle only — excludes diagonal (self-loops) and
      # avoids double counting from symmetry
      upper        <- upper.tri(A_curr)
      edges_curr   <- sum(A_curr[upper])
      edges_next   <- sum(A_next[upper])
      intersection <- sum(A_curr[upper] & A_next[upper])
      union_edges  <- sum(A_curr[upper] | A_next[upper])
      
      if (union_edges == 0) {
        # Both graphs have no edges — define similarity as 1 (they agree)
        jaccard_val <- 1
        message(sprintf("    [%s %s] week %d -> %d: both empty, similarity=1",
                        ID, d_level, w_curr, w_next))
      } else {
        jaccard_val <- intersection / union_edges
        message(sprintf("    [%s %s] week %d -> %d: edges_curr=%g, edges_next=%g, intersection=%g, union=%g, similarity=%.4f",
                        ID, d_level, w_curr, w_next,
                        edges_curr, edges_next, intersection, union_edges, jaccard_val))
      }
      
      similarity_vals  <- c(similarity_vals,  jaccard_val)
      similarity_weeks <- c(similarity_weeks, w_curr)
    }
    
    if (length(similarity_weeks) == 0) return(NULL)
    
    data.frame(
      week       = similarity_weeks,
      similarity = similarity_vals,
      mouse      = ID
    )
  }
  
  # --------------------------------------------------------------------------
  # Compute cross-stratum Jaccard similarity (same week, two strata) for
  # one mouse. Returns a data.frame or NULL.
  # --------------------------------------------------------------------------
  compute_cross_stratum_jaccard <- function(ID) {
    
    adj_A <- load_adj_by_week(ID, d_level_A)
    adj_B <- load_adj_by_week(ID, d_level_B)
    
    if (is.null(adj_A) || is.null(adj_B)) return(NULL)
    
    # Weeks present in both strata — sorted numerically
    common_weeks <- sort(intersect(as.numeric(names(adj_A)),
                                   as.numeric(names(adj_B))))
    
    if (length(common_weeks) == 0) {
      message(sprintf("    [%s] no common weeks between %s and %s — skipping",
                      ID, d_level_A, d_level_B))
      return(NULL)
    }
    
    similarity_vals <- numeric(length(common_weeks))
    
    for (w_idx in seq_along(common_weeks)) {
      w      <- common_weeks[w_idx]
      A_curr <- adj_A[[as.character(w)]]
      B_curr <- adj_B[[as.character(w)]]
      
      # Guard against NULL lookup
      if (is.null(A_curr) || is.null(B_curr)) {
        message(sprintf("    [%s] week %d: NULL adj_mat in one stratum — skipping",
                        ID, w))
        similarity_vals[w_idx] <- NA
        next
      }
      
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
        message(sprintf("    [%s] week %d: both strata empty, similarity=1", ID, w))
      } else {
        jaccard_val <- intersection / union_edges
        message(sprintf("    [%s] week %d: edges_A=%g, edges_B=%g, intersection=%g, union=%g, similarity=%.4f",
                        ID, w, edges_A, edges_B, intersection, union_edges, jaccard_val))
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
  # Collect all three panels' data across all mice.
  # Panel labels are assigned AFTER rbind to prevent factor-level
  # contamination: if labels were set before rbind as factors with different
  # level sets, rows could silently receive wrong panel assignments.
  # --------------------------------------------------------------------------
  message(sprintf("Computing temporal Jaccard: %s", d_level_A))
  rows_A <- do.call(rbind, Filter(Negate(is.null),
                                  lapply(all_IDs_A, compute_temporal_jaccard, d_level = d_level_A)))
  
  message(sprintf("Computing temporal Jaccard: %s", d_level_B))
  rows_B <- do.call(rbind, Filter(Negate(is.null),
                                  lapply(all_IDs_B, compute_temporal_jaccard, d_level = d_level_B)))
  
  message(sprintf("Computing cross-stratum Jaccard: %s vs %s", d_level_A, d_level_B))
  cross_IDs  <- intersect(all_IDs_A, all_IDs_B)
  rows_cross <- do.call(rbind, Filter(Negate(is.null),
                                      lapply(cross_IDs, compute_cross_stratum_jaccard)))
  
  # Assign panel labels as plain character AFTER each block is assembled,
  # then combine — this guarantees no row inherits a label from a different block
  if (!is.null(rows_A))     rows_A$panel     <- panel_labels[["A_temporal"]]
  if (!is.null(rows_B))     rows_B$panel     <- panel_labels[["B_temporal"]]
  if (!is.null(rows_cross)) rows_cross$panel <- panel_labels[["cross"]]
  
  plot_df       <- do.call(rbind, list(rows_A, rows_B, rows_cross))
  # Convert to factors only after rbind, with verified level sets
  plot_df$mouse <- factor(plot_df$mouse, levels = all_IDs)
  plot_df$panel <- factor(plot_df$panel, levels = panel_order)
  
  # Sanity check: print row counts per mouse x panel to confirm no contamination
  message("[plot_edge_stability_combined] row counts per mouse x panel:")
  print(table(plot_df$mouse, plot_df$panel, useNA = "always"))
  
  # --------------------------------------------------------------------------
  # Color palette: Tau in warm tones, WT in cool tones
  # --------------------------------------------------------------------------
  tau_mice  <- all_IDs[grepl("^Tau", all_IDs)]
  wt_mice   <- all_IDs[grepl("^WT",  all_IDs)]
  tau_cols  <- setNames(scales::hue_pal(h = c(0, 60))(max(length(tau_mice), 1)),   tau_mice)
  wt_cols   <- setNames(scales::hue_pal(h = c(200, 260))(max(length(wt_mice), 1)), wt_mice)
  color_map <- c(tau_cols, wt_cols)
  
  # --------------------------------------------------------------------------
  # Pre-compute loess fits separately per mouse x panel.
  # This is necessary because geom_smooth pools all data before faceting,
  # fitting loess across all three panels rather than within each.
  # NA similarity rows are dropped before fitting.
  # --------------------------------------------------------------------------
  loess_lines <- do.call(rbind, lapply(
    split(plot_df, list(plot_df$mouse, plot_df$panel), drop = TRUE),
    function(df) {
      df <- df[!is.na(df$similarity), ]
      if (nrow(df) < 4) return(NULL)
      week_seq <- seq(min(df$week), max(df$week), length.out = 200)
      fit <- tryCatch(
        predict(loess(similarity ~ week, data = df), newdata = data.frame(week = week_seq)),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      data.frame(
        week       = week_seq,
        similarity = fit,
        mouse      = df$mouse[1],
        panel      = df$panel[1]
      )
    }
  ))
  loess_lines$mouse <- factor(loess_lines$mouse, levels = all_IDs)
  loess_lines$panel <- factor(loess_lines$panel, levels = panel_order)
  
  # Confirm loess split respected panel boundaries
  message("[plot_edge_stability_combined] loess row counts per mouse x panel:")
  print(table(loess_lines$mouse, loess_lines$panel, useNA = "always"))
  
  # --------------------------------------------------------------------------
  # Build the shared base plot: loess curves drawn from pre-computed
  # predictions via geom_line, faceted by panel and stacked vertically
  # (ncol = 1). Legend placed below all panels.
  # --------------------------------------------------------------------------
  base_plot <- ggplot(loess_lines,
                      aes(x = week, y = similarity, color = mouse, group = mouse)) +
    geom_line(size = 0.9) +
    facet_wrap(~ panel, ncol = 1) +
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
  
  # Linear y-axis
  g <- base_plot +
    scale_y_continuous(
      breaks = c(0, 0.2, 0.4, 0.6),
      limits = c(-0.1, 0.7)
    )
  
  # Sqrt-transformed y-axis
  g_sqrt <- base_plot +
    scale_y_continuous(
      trans  = "sqrt",
      breaks = c(0, 0.2, 0.4, 0.6),
      limits = c(-0.1, 0.7)
    )
  
  return(list(linear = g, sqrt = g_sqrt))
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
        geom_line(size = 0.8) +
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

# facet wrapped both strata + 3 regional connections into 6 plots in one.
plot_edge_regional_proportion_all_mice_v2 <- function(results_folder, time_scale, discrete_levels) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: For each discrete stratum, plot the proportion of existing edges that
  #       are HIP-HIP, HIP-EHC, or EHC-EHC across all discovered mice.
  #
  #       n_hip inferred per mouse from res_i$recovery_params$kept_neuron_regions
  #       (same source used by visualize_strata_all_mice_v2).
  #       All computations use upper.tri() only.
  #
  #       Facet layout: rows = region pair (HIP-HIP, HIP-EHC, EHC-EHC),
  #                     cols = discrete stratum (e.g. Resting, Running).
  #       Loess is pre-computed separately per mouse x stratum x region_pair
  #       to guarantee within-panel fitting (geom_smooth pools across facets).
  #
  # inputs:
  #   results_folder   (string)
  #   time_scale       (integer)  e.g. 10
  #   discrete_levels  (vector of strings)  e.g. c('m0vr1', 'm1vr1')
  #
  # returns: a single ggplot object (faceted grid)
  #
  # ----------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Stratum display labels
  # --------------------------------------------------------------------------
  stratum_labels <- c(m0vr1 = "Resting", m1vr1 = "Running")
  
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
  
  # --------------------------------------------------------------------------
  # Collect data across all discrete levels into one long data frame
  # --------------------------------------------------------------------------
  all_data <- data.frame()
  
  for (lvl in discrete_levels) {
    
    ids <- discover_IDs(lvl)
    message(sprintf("[plot_edge_regional_proportion] level=%s, found %d IDs: %s",
                    lvl, length(ids), paste(ids, collapse = ", ")))
    
    # Apply display label if available, otherwise keep original level name
    stratum_label <- ifelse(lvl %in% names(stratum_labels),
                            stratum_labels[[lvl]],
                            lvl)
    
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
            week = week, HIP_HIP = NA_real_, HIP_EHC = NA_real_, EHC_EHC = NA_real_,
            mouse = ID, stratum = stratum_label
          ))
          next
        }
        
        ut <- which(upper.tri(adj_mat) & adj_mat != 0, arr.ind = TRUE)
        
        if (nrow(ut) == 0) {
          all_data <- rbind(all_data, data.frame(
            week = week, HIP_HIP = NA_real_, HIP_EHC = NA_real_, EHC_EHC = NA_real_,
            mouse = ID, stratum = stratum_label
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
            mouse   = ID,
            stratum = stratum_label
          ))
        }
      }
      
      rm(tmp_env, res_i, step_data)
    }
  }
  
  if (nrow(all_data) == 0) {
    message("[WARN] no data found across all levels")
    return(NULL)
  }
  
  # --------------------------------------------------------------------------
  # Pivot to long format so region_pair becomes a facet dimension
  # --------------------------------------------------------------------------
  long_data <- reshape(
    all_data,
    varying      = c("HIP_HIP", "HIP_EHC", "EHC_EHC"),
    v.names      = "proportion",
    timevar      = "region_pair",
    times        = c("HIP-HIP", "HIP-EHC", "EHC-EHC"),
    direction    = "long"
  )
  
  long_data$mouse       <- factor(long_data$mouse,       levels = all_IDs)
  long_data$region_pair <- factor(long_data$region_pair, levels = region_pairs)
  
  # Stratum column order follows discrete_levels order, mapped through labels
  stratum_order <- ifelse(discrete_levels %in% names(stratum_labels),
                          stratum_labels[discrete_levels],
                          discrete_levels)
  long_data$stratum <- factor(long_data$stratum, levels = stratum_order)
  
  # --------------------------------------------------------------------------
  # Pre-compute loess fits separately per mouse x stratum x region_pair.
  # This is necessary because geom_smooth pools all data before faceting,
  # so it would fit loess across all strata and region pairs rather than
  # within each panel. Manually fitting guarantees within-panel curves.
  # NA rows (weeks with no edges) are dropped before fitting.
  # --------------------------------------------------------------------------
  loess_lines <- do.call(rbind, lapply(
    split(long_data, list(long_data$mouse, long_data$stratum, long_data$region_pair), drop = TRUE),
    function(df) {
      df <- df[!is.na(df$proportion), ]   # drop weeks with no edges
      if (nrow(df) < 4) return(NULL)
      week_seq <- seq(min(df$week), max(df$week), length.out = 200)
      fit <- tryCatch(
        predict(loess(proportion ~ week, data = df), newdata = data.frame(week = week_seq)),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      data.frame(
        week        = week_seq,
        proportion  = fit,
        mouse       = df$mouse[1],
        stratum     = df$stratum[1],
        region_pair = df$region_pair[1]
      )
    }
  ))
  loess_lines$mouse       <- factor(loess_lines$mouse,       levels = all_IDs)
  loess_lines$region_pair <- factor(loess_lines$region_pair, levels = region_pairs)
  loess_lines$stratum     <- factor(loess_lines$stratum,     levels = stratum_order)
  
  # --------------------------------------------------------------------------
  # Build faceted plot: rows = region pair, cols = stratum.
  # Curves drawn from manually pre-computed loess predictions (geom_line)
  # to ensure within-panel fitting. Legend below all panels.
  # --------------------------------------------------------------------------
  p_plot <- ggplot(loess_lines, aes(x = week, y = proportion,
                                    color = mouse, group = mouse)) +
    # Loess curves pre-computed per mouse x stratum x region_pair
    geom_line(size = 0.9) +
    facet_grid(region_pair ~ stratum) +
    scale_color_manual(values = color_map) +
    scale_x_continuous(breaks = c(20, 25, 30, 35)) +
    scale_y_continuous(
      limits = c(-0.05, 1),
      breaks = c(0, 0.5, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      x     = "Age (Weeks)",
      y     = "Proportion of Edges",
      color = "Mouse"
    ) +
    guides(color = guide_legend(nrow = 1)) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid      = element_blank(),
      legend.position = "bottom"
    )
  
  message("[OK] combined regional proportion plot created")
  return(p_plot)
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
    
    # Reset for each stratum
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
    
    message(sprintf("  [%s] all_data has %d rows, mice: %s",
                    lvl, nrow(all_data), paste(unique(all_data$mouse), collapse = ", ")))
    print(head(all_data))
    
    all_data$mouse <- factor(all_data$mouse, levels = all_IDs)
    
    # Pre-compute loess-smoothed q25/q75 per mouse for the ribbon.
    # geom_ribbon does not support stat="smooth" with ymin/ymax aesthetics
    # (it requires a y aesthetic), so we fit loess manually on a fine grid
    # and pass the smoothed band as a separate data frame.
    smooth_band <- do.call(rbind, lapply(split(all_data, all_data$mouse), function(df) {
      if (nrow(df) < 4) return(NULL)
      week_seq <- seq(min(df$week), max(df$week), length.out = 100)
      lo25 <- tryCatch(predict(loess(q25 ~ week, data = df), newdata = data.frame(week = week_seq)), error = function(e) NULL)
      lo75 <- tryCatch(predict(loess(q75 ~ week, data = df), newdata = data.frame(week = week_seq)), error = function(e) NULL)
      if (is.null(lo25) || is.null(lo75)) return(NULL)
      data.frame(week = week_seq, q25_smooth = lo25, q75_smooth = lo75, mouse = df$mouse[1])
    }))
    smooth_band$mouse <- factor(smooth_band$mouse, levels = all_IDs)
    
    p_plot <- ggplot(all_data, aes(x = week, color = mouse, fill = mouse, group = mouse)) +
      # Raw points to confirm data is present
      geom_point(aes(y = median), size = 1.5, alpha = 0.5) +
      # Shaded IQR band from manually pre-computed loess on q25 and q75
      geom_ribbon(
        data = smooth_band,
        aes(ymin = q25_smooth, ymax = q75_smooth),
        alpha = 0.15, color = NA
      ) +
      # Median loess line
      geom_smooth(
        aes(y = median),
        method = "loess", formula = y ~ x,
        se = FALSE, size = 0.9
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

plot_median_nonzero_degree_all_mice_v2 <- function(results_folder, time_scale, discrete_levels) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: For each discrete stratum, and for each mouse, compute the median
  #       degree across neurons that have at least one edge (positive degree)
  #       at each week. Plots one loess curve per mouse (fit to the median),
  #       with no confidence interval or IQR ribbon.
  #
  #       Loess is fit separately within each mouse x stratum combination
  #       to match the per-stratum fits of the v1 function.
  #       All discrete strata are combined into a single faceted plot,
  #       stacked vertically (ncol = 1).
  #
  # inputs:
  #   results_folder   (string)
  #   time_scale       (integer)  e.g. 10
  #   discrete_levels  (vector of strings)  e.g. c('m0vr1', 'm1vr1')
  #
  # returns: a single ggplot object (faceted vertically by stratum)
  #
  # ----------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Stratum display labels
  # --------------------------------------------------------------------------
  stratum_labels <- c(m0vr1 = "Resting", m1vr1 = "Running")
  
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
  
  # --------------------------------------------------------------------------
  # Collect data across all strata into one long data frame.
  # unname() on quantile() strips the "50%" name that would otherwise
  # attach as a row name and cause silent misalignment across rbind calls.
  # --------------------------------------------------------------------------
  all_data <- data.frame()
  
  for (lvl in discrete_levels) {
    
    ids <- discover_IDs(lvl)
    message(sprintf("[plot_median_nonzero_degree] level=%s, found %d IDs: %s",
                    lvl, length(ids), paste(ids, collapse = ", ")))
    
    # Apply display label if available, otherwise keep original level name
    stratum_label <- ifelse(lvl %in% names(stratum_labels),
                            stratum_labels[[lvl]],
                            lvl)
    
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
          week    = week,
          median  = unname(quantile(nonzero_degree, 0.50)),
          mouse   = ID,
          stratum = stratum_label
        ))
      }
      
      rm(tmp_env, res_i, step_data)
    }
  }
  
  if (nrow(all_data) == 0) {
    message("[WARN] no data found across all levels")
    return(NULL)
  }
  
  message(sprintf("[plot_median_nonzero_degree] all_data has %d rows, mice: %s",
                  nrow(all_data), paste(unique(all_data$mouse), collapse = ", ")))
  print(head(all_data))
  
  all_data$mouse <- factor(all_data$mouse, levels = all_IDs)
  
  # Derive ordered factor levels from discrete_levels order, mapped through labels
  stratum_order <- ifelse(discrete_levels %in% names(stratum_labels),
                          stratum_labels[discrete_levels],
                          discrete_levels)
  all_data$stratum <- factor(all_data$stratum, levels = stratum_order)
  
  # --------------------------------------------------------------------------
  # Pre-compute loess fits separately per mouse x stratum combination.
  # This is necessary because geom_smooth pools all data before faceting,
  # so it would fit loess across both strata rather than within each panel.
  # Manually fitting and passing predicted values guarantees each curve
  # matches the per-stratum v1 fits exactly.
  # --------------------------------------------------------------------------
  loess_lines <- do.call(rbind, lapply(
    split(all_data, list(all_data$mouse, all_data$stratum), drop = TRUE),
    function(df) {
      if (nrow(df) < 4) return(NULL)
      week_seq  <- seq(min(df$week), max(df$week), length.out = 200)
      fit       <- tryCatch(
        predict(loess(median ~ week, data = df), newdata = data.frame(week = week_seq)),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      data.frame(
        week    = week_seq,
        median  = fit,
        mouse   = df$mouse[1],
        stratum = df$stratum[1]
      )
    }
  ))
  loess_lines$mouse   <- factor(loess_lines$mouse,   levels = all_IDs)
  loess_lines$stratum <- factor(loess_lines$stratum, levels = stratum_order)
  
  # --------------------------------------------------------------------------
  # Build single faceted plot: stacked vertically by stratum (ncol = 1).
  # Curves drawn from manually pre-computed loess predictions (geom_line)
  # to ensure within-stratum fitting. No CI, no ribbon, no raw points.
  # Legend placed below all panels.
  # --------------------------------------------------------------------------
  p_plot <- ggplot(loess_lines, aes(x = week, y = median, color = mouse, group = mouse)) +
    # Loess curves pre-computed per mouse x stratum — not fit across strata
    geom_line(size = 0.9) +
    facet_wrap(~ stratum, ncol = 1) +
    scale_color_manual(values = color_map) +
    scale_x_continuous(breaks = c(20, 25, 30, 35)) +
    scale_y_continuous(
      limits = c(-0.05, 0.65),
      breaks = c(0.0, 0.3, 0.6)
    ) +
    labs(
      x     = "Age (Weeks)",
      y     = "Normalized Degree",
      color = "Mouse"
    ) +
    guides(color = guide_legend(nrow = 1)) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid      = element_blank(),
      legend.position = "bottom"
    )
  
  message("[OK] combined normalized degree plot created")
  return(p_plot)
}

# median degree including zeros. Just your basic calculation.
plot_median_degree_all_mice_v3 <- function(results_folder, time_scale, discrete_levels) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: For each discrete stratum, and for each mouse, compute the median
  #       normalized degree across ALL neurons (including zero-degree nodes)
  #       at each week. Degree is normalized by (n_neurons - 1). Plots one
  #       loess curve per mouse, no CI or ribbon. All discrete strata combined
  #       into a single faceted plot, stacked vertically (ncol = 1).
  #
  # inputs:
  #   results_folder   (string)
  #   time_scale       (integer)  e.g. 10
  #   discrete_levels  (vector of strings)  e.g. c('m0vr1', 'm1vr1')
  #
  # returns: a single ggplot object (faceted vertically by stratum)
  #
  # ----------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Stratum display labels
  # --------------------------------------------------------------------------
  stratum_labels <- c(m0vr1 = "Resting", m1vr1 = "Running")
  
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
  
  # --------------------------------------------------------------------------
  # Collect data across all strata into one long data frame.
  # All neurons included — zero-degree nodes are NOT filtered out.
  # unname() on quantile() strips the "50%" name.
  # --------------------------------------------------------------------------
  all_data <- data.frame()
  
  for (lvl in discrete_levels) {
    
    ids <- discover_IDs(lvl)
    message(sprintf("[plot_median_degree] level=%s, found %d IDs: %s",
                    lvl, length(ids), paste(ids, collapse = ", ")))
    
    stratum_label <- ifelse(lvl %in% names(stratum_labels),
                            stratum_labels[[lvl]],
                            lvl)
    
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
        
        # Degree = row sums, excluding diagonal; normalized by (n - 1)
        # All neurons included — zero-degree nodes kept
        diag(adj_mat) <- 0
        degree <- rowSums(adj_mat) / (nrow(adj_mat) - 1)
        
        all_data <- rbind(all_data, data.frame(
          week    = week,
          median  = unname(quantile(degree, 0.50)),
          mouse   = ID,
          stratum = stratum_label
        ))
      }
      
      rm(tmp_env, res_i, step_data)
    }
  }
  
  if (nrow(all_data) == 0) {
    message("[WARN] no data found across all levels")
    return(NULL)
  }
  
  message(sprintf("[plot_median_degree] all_data has %d rows, mice: %s",
                  nrow(all_data), paste(unique(all_data$mouse), collapse = ", ")))
  print(head(all_data))
  
  all_data$mouse <- factor(all_data$mouse, levels = all_IDs)
  
  stratum_order <- ifelse(discrete_levels %in% names(stratum_labels),
                          stratum_labels[discrete_levels],
                          discrete_levels)
  all_data$stratum <- factor(all_data$stratum, levels = stratum_order)
  
  # --------------------------------------------------------------------------
  # Pre-compute loess fits separately per mouse x stratum combination.
  # --------------------------------------------------------------------------
  loess_lines <- do.call(rbind, lapply(
    split(all_data, list(all_data$mouse, all_data$stratum), drop = TRUE),
    function(df) {
      if (nrow(df) < 4) return(NULL)
      week_seq <- seq(min(df$week), max(df$week), length.out = 200)
      fit <- tryCatch(
        predict(loess(median ~ week, data = df), newdata = data.frame(week = week_seq)),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      data.frame(
        week    = week_seq,
        median  = fit,
        mouse   = df$mouse[1],
        stratum = df$stratum[1]
      )
    }
  ))
  loess_lines$mouse   <- factor(loess_lines$mouse,   levels = all_IDs)
  loess_lines$stratum <- factor(loess_lines$stratum, levels = stratum_order)
  
  # --------------------------------------------------------------------------
  # Build single faceted plot: stacked vertically by stratum (ncol = 1).
  # --------------------------------------------------------------------------
  p_plot <- ggplot(loess_lines, aes(x = week, y = median, color = mouse, group = mouse)) +
    geom_line(size = 0.9) +
    facet_wrap(~ stratum, ncol = 1) +
    scale_color_manual(values = color_map) +
    scale_x_continuous(breaks = c(20, 25, 30, 35)) +
    scale_y_continuous(limits = c(-0.05, 0.65), breaks = c(0.0, 0.2, 0.4, 0.6))
    labs(
      x     = "Age (Weeks)",
      y     = "Median Normalized Degree",
      color = "Mouse"
    ) +
    guides(color = guide_legend(nrow = 1)) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid      = element_blank(),
      legend.position = "bottom"
    )
  
  message("[OK] combined median degree plot created")
  return(p_plot)
}


# ridgeline graph for degree distribution (cut off the part of the histogram with 0's)

plot_degree_ridgeline_all_mice <- function(results_folder, time_scale, discrete_levels) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: For each mouse x stratum combination, plot a ridgeline (joy plot)
  #       of the non-zero normalized degree distribution over time. Each
  #       ridge corresponds to one week; ridges are stacked by week (oldest
  #       at bottom, most recent at top). Area under each ridge is proportional
  #       to the number of non-zero degree neurons at that week (i.e. density
  #       is NOT normalized to integrate to 1 — it is scaled by count).
  #       A gray vertical segment marks the median of each ridge.
  #       Faint gray dotted vertical reference lines at 0, 0.25, 0.5, 0.75, 1.
  #       A light shaded vertical band spans the IQR of per-week medians for
  #       each mouse x stratum panel.
  #
  #       Layout: rows = mice (Tau first, then WT), cols = strata.
  #       Uses ggridges::geom_density_ridges with stat="density" and
  #       after_stat(count) scaling.
  #
  # inputs:
  #   results_folder   (string)
  #   time_scale       (integer)  e.g. 10
  #   discrete_levels  (vector of strings)  e.g. c('m0vr1', 'm1vr1')
  #
  # returns: a single ggplot object (faceted by mouse x stratum)
  #
  # ----------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Stratum display labels
  # --------------------------------------------------------------------------
  stratum_labels <- c(m0vr1 = "Resting", m1vr1 = "Running")
  
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
  
  # --------------------------------------------------------------------------
  # Collect raw non-zero degree values (one row per neuron x week x mouse x
  # stratum). We keep the individual degree values — NOT per-week summaries —
  # so that ggridges can fit a kernel density estimate per week ridge.
  # --------------------------------------------------------------------------
  all_data <- data.frame()
  
  for (lvl in discrete_levels) {
    
    ids <- discover_IDs(lvl)
    message(sprintf("[plot_degree_ridgeline] level=%s, found %d IDs: %s",
                    lvl, length(ids), paste(ids, collapse = ", ")))
    
    stratum_label <- ifelse(lvl %in% names(stratum_labels),
                            stratum_labels[[lvl]],
                            lvl)
    
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
      
      message(sprintf("  [%s | %s] weeks=%d", ID, lvl, length(y_c_weeks)))
      
      for (idx in seq_along(y_c_weeks)) {
        
        week    <- y_c_weeks[idx]
        adj_mat <- step_data[[idx]][["adj_mat_KL_GIC_local_est_eig3"]]
        
        if (is.null(adj_mat)) next
        
        # Normalized degree for all neurons; keep only positive-degree neurons
        diag(adj_mat) <- 0
        degree         <- rowSums(adj_mat) / (nrow(adj_mat) - 1)
        nonzero_degree <- degree[degree > 0]
        
        if (length(nonzero_degree) == 0) next
        
        # One row per neuron — raw values needed for per-ridge KDE
        all_data <- rbind(all_data, data.frame(
          degree  = nonzero_degree,
          week    = week,
          mouse   = ID,
          stratum = stratum_label
        ))
      }
      
      rm(tmp_env, res_i, step_data)
    }
  }
  
  if (nrow(all_data) == 0) {
    message("[WARN] no data found across all levels")
    return(NULL)
  }
  
  message(sprintf("[plot_degree_ridgeline] all_data has %d rows, mice: %s",
                  nrow(all_data), paste(unique(all_data$mouse), collapse = ", ")))
  print(head(all_data))
  
  # Ordered factors: mice rows top-to-bottom (Tau1..Tau3, WT1..WT3)
  # week kept numeric so scale_y_continuous breaks work correctly
  all_data$mouse   <- factor(all_data$mouse, levels = all_IDs)
  all_data$week    <- as.numeric(as.character(all_data$week))
  
  stratum_order    <- ifelse(discrete_levels %in% names(stratum_labels),
                             stratum_labels[discrete_levels],
                             discrete_levels)
  all_data$stratum <- factor(all_data$stratum, levels = stratum_order)
  
  # --------------------------------------------------------------------------
  # Pre-compute per-week median per mouse x stratum for overlay segments.
  # geom_density_ridges stacks ridges at each y = week value; we draw a
  # short vertical segment at x = median, centered on the ridge baseline.
  # --------------------------------------------------------------------------
  median_data <- do.call(rbind, lapply(
    split(all_data, list(all_data$mouse, all_data$stratum, all_data$week), drop = TRUE),
    function(df) {
      data.frame(
        week    = df$week[1],
        median  = median(df$degree),
        mouse   = df$mouse[1],
        stratum = df$stratum[1]
      )
    }
  ))
  median_data$mouse   <- factor(median_data$mouse,   levels = all_IDs)
  median_data$stratum <- factor(median_data$stratum, levels = stratum_order)
  
  # --------------------------------------------------------------------------
  # Pre-compute IQR of per-week medians per mouse x stratum, for the shaded
  # vertical band. One row per panel (mouse x stratum); xmin/xmax are the
  # 25th and 75th percentiles of the weekly median values in that panel.
  # ymin/ymax span the full plot range so geom_rect fills the panel height.
  # --------------------------------------------------------------------------
  iqr_data <- do.call(rbind, lapply(
    split(median_data, list(median_data$mouse, median_data$stratum), drop = TRUE),
    function(df) {
      data.frame(
        xmin    = unname(quantile(df$median, 0.25)),
        xmax    = unname(quantile(df$median, 0.75)),
        mouse   = df$mouse[1],
        stratum = df$stratum[1]
      )
    }
  ))
  iqr_data$mouse   <- factor(iqr_data$mouse,   levels = all_IDs)
  iqr_data$stratum <- factor(iqr_data$stratum, levels = stratum_order)
  
  # --------------------------------------------------------------------------
  # Ridgeline plot with median overlay and IQR shaded band.
  # Layer order (back to front):
  #   1. geom_vline       — dotted reference lines at 0, 0.25, 0.5, 0.75, 1
  #   2. geom_rect        — shaded IQR band of per-week medians
  #   3. geom_density_ridges — ridge fills and outlines
  #   4. geom_segment     — per-ridge median tick
  # --------------------------------------------------------------------------
  p_plot <- ggplot(all_data,
                   aes(x      = degree,
                       y      = week,
                       fill   = mouse,
                       color  = mouse,
                       height = after_stat(count),
                       group  = week)) +
    # 1. Faint dotted reference lines at canonical normalized-degree values
    geom_vline(
      xintercept = c(0, 0.25, 0.5, 0.75, 1),
      color      = "gray50",
      linewidth  = 0.35,
      linetype   = "dashed"
    ) +
    # 2. Shaded vertical band spanning IQR of per-week medians per panel;
    #    ymin/ymax use -Inf/Inf so the band spans the full panel height
    geom_rect(
      data        = iqr_data,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill        = "gray70",
      alpha       = 0.25,
      inherit.aes = FALSE
    ) +
    # 3. Ridge densities scaled by count
    ggridges::geom_density_ridges(
      stat           = "density",
      alpha          = 0.4,
      linewidth      = 0.4,
      scale          = 0.9,
      rel_min_height = 0.01
    ) +
    # 4. Median vertical tick per ridge: segment rises from ridge baseline
    #    (y = week) to slightly above it (yend = week + 0.8)
    geom_segment(
      data      = median_data,
      aes(x = median, xend = median, y = week, yend = week + 0.8),
      color     = "gray40",
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    facet_grid(mouse ~ stratum) +
    scale_fill_manual(values  = color_map) +
    scale_color_manual(values = color_map) +
    scale_x_continuous(limits = c(0, NA)) +
    scale_y_continuous(breaks = c(20, 25, 30, 35)) +
    labs(
      x     = "Normalized Degree",
      y     = "Age (Weeks)",
      fill  = "Mouse",
      color = "Mouse"
    ) +
    guides(fill  = guide_legend(nrow = 1),
           color = guide_legend(nrow = 1)) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid      = element_blank(),
      legend.position = "bottom",
      strip.text      = element_text(size = 13)
    )
  
  message("[OK] ridgeline degree plot created")
  return(p_plot)
}


