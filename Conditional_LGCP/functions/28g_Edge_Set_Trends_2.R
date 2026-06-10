
# for one mouse, plot its edge sets across both strata

visualize_edge_set_one_mouse_all_strata <- function(results_folder, all_weeks, time_scale, discrete_levels, output, region_border, mouse_ID, display_weeks = all_weeks) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: For a single mouse, plot all available discrete strata as rows in a
  #       single facet_grid. Rows = stratum label, Columns = week.
  #       One set of axis labels. Stratum labels are human-readable
  #       (m0vr1 -> "Resting", m1vr1 -> "Running"; others passed through).
  #
  # inputs:
  #
  # - results_folder   (string)
  # - all_weeks        (vector)  full set of possible weeks, e.g. 17:38;
  #                              used only to compute absent_weeks
  # - time_scale       (integer)  e.g. 10
  # - discrete_levels  (vector of strings)  e.g. c("m0vr1", "m1vr1")
  # - output           (string)   'adj', 'P_HS', or 'C_HS'
  # - region_border    (boolean)
  # - mouse_ID         (string)   e.g. "Tau1", "WT2"
  # - display_weeks    (vector)  subset of weeks to actually show, e.g. 17:22;
  #                              defaults to all_weeks
  #
  # returns: named list with two elements:
  #   - plot         : single ggplot object (all strata as rows)
  #   - found_levels : character vector of discrete levels actually loaded
  #
  # Color scheme for output == 'adj':
  #   HIP-HIP edges  -> emerald teal  (#00C896)
  #   EHC-EHC edges  -> vivid purple  (#CC3FFF)
  #   HIP-EHC edges  -> periwinkle    (#6694CC)  [sqrt-mean-squares RGB fusion]
  #
  # Region assignment: nodes 1..boundary are HIP, nodes (boundary+1)..max_node
  # are EHC, where boundary = floor(first element of boundaries vector).
  # If no boundary is available, all nodes are treated as HIP-HIP.
  #
  # Col_ID factor levels are set to display_weeks so that exactly those columns
  # appear. Weeks outside display_weeks are never rendered.
  #
  # ----------------------------------------------------------------------------
  
  # Color constants for adj output
  COL_HIP_HIP <- "#00C896"   # emerald teal  (HIP)
  COL_EHC_EHC <- "#CC3FFF"   # vivid purple  (EHC)
  COL_HIP_EHC <- "#6694CC"   # periwinkle    (HIP-EHC sqrt-mean-squares fusion)
  
  # Human-readable stratum label lookup
  stratum_labels <- c(
    m0vr1 = "Resting",
    m1vr1 = "Running"
  )
  pretty_label <- function(d_level) {
    if (d_level %in% names(stratum_labels)) stratum_labels[[d_level]] else d_level
  }
  
  # --------------------------------------------------------------------------
  # Load sparse data for one mouse x one discrete level.
  # Only weeks in display_weeks are retained.
  # --------------------------------------------------------------------------
  load_sparse_data <- function(d_level) {
    file_path <- paste0(results_folder, '/', mouse_ID, '_', d_level, '_t', time_scale, '.RData')
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      return(NULL)
    }
    tmp_env <- new.env()
    load(file_path, envir = tmp_env)
    res_i <- tmp_env$graph_results_i
    if (region_border) {
      boundaries <- get_factor_boundaries(res_i$recovery_params$kept_neuron_regions)
    } else {
      boundaries <- numeric(0)
    }
    y_c_weeks    <- as.numeric(res_i$y_c_query)
    # absent_weeks: weeks in display_weeks that were not estimated
    absent_weeks <- setdiff(display_weeks, y_c_weeks)
    sparse_data  <- list()
    if (output == 'adj') {
      step_data <- res_i$step_12b
      for (idx in seq_along(y_c_weeks)) {
        current_week <- y_c_weeks[idx]
        # skip weeks outside the display window
        if (!(current_week %in% display_weeks)) next
        adj_mat <- step_data[[idx]][["adj_mat_KL_GIC_local_est_eig3"]]
        if (is.null(adj_mat)) next
        coords <- which(adj_mat == 1, arr.ind = TRUE)
        if (nrow(coords) > 0) coords <- coords[coords[, 2] > coords[, 1], , drop = FALSE]
        if (nrow(coords) > 0) {
          df_coords <- as.data.frame(coords)
          colnames(df_coords) <- c("Node_Row", "Node_Col")
          sparse_data[[as.character(current_week)]] <- df_coords
        } else {
          sparse_data[[as.character(current_week)]] <- data.frame(Node_Row = integer(0),
                                                                  Node_Col = integer(0))
        }
      }
    } else if (output == 'P_HS') {
      step_data <- res_i$step_11
      for (idx in seq_along(y_c_weeks)) {
        current_week <- y_c_weeks[idx]
        if (!(current_week %in% display_weeks)) next
        P_HS_log      <- log(step_data[[idx]][["w_mat_KL_est_eig3"]])
        if (is.null(P_HS_log)) next
        upper_tri_idx <- which(col(P_HS_log) >= row(P_HS_log), arr.ind = TRUE)
        df_coords     <- as.data.frame(upper_tri_idx)
        colnames(df_coords) <- c("Node_Row", "Node_Col")
        df_coords$Value <- P_HS_log[upper_tri_idx]
        sparse_data[[as.character(current_week)]] <- df_coords
      }
    } else if (output == 'C_HS') {
      step_data <- res_i$step_11b
      for (idx in seq_along(y_c_weeks)) {
        current_week <- y_c_weeks[idx]
        if (!(current_week %in% display_weeks)) next
        C_HS_log      <- log(step_data[[idx]][["C_HS_KL_est_eig3"]])
        if (is.null(C_HS_log)) next
        upper_tri_idx <- which(col(C_HS_log) >= row(C_HS_log), arr.ind = TRUE)
        df_coords     <- as.data.frame(upper_tri_idx)
        colnames(df_coords) <- c("Node_Row", "Node_Col")
        df_coords$Value <- C_HS_log[upper_tri_idx]
        sparse_data[[as.character(current_week)]] <- df_coords
      }
    } else {
      stop('error in visualize_edge_set_one_mouse_all_strata: unknown output type')
    }
    rm(tmp_env, res_i, step_data)
    return(list(sparse_data = sparse_data, absent_weeks = absent_weeks, boundaries = boundaries))
  }
  
  # --------------------------------------------------------------------------
  # Helper: get max node index for a mouse from its sparse data
  # --------------------------------------------------------------------------
  get_max_node <- function(sparse_data) {
    max_node <- 0
    for (wk in sparse_data) {
      if (nrow(wk) > 0) max_node <- max(max_node, max(wk$Node_Row, wk$Node_Col))
    }
    max(max_node, 1)
  }
  
  # --------------------------------------------------------------------------
  # Helper: classify an edge (row_node, col_node) given the HIP boundary.
  # All nodes <= hip_boundary are HIP; nodes > hip_boundary are EHC.
  # Returns one of: "HIP_HIP", "EHC_EHC", "HIP_EHC"
  # If hip_boundary is NA (no boundary available), returns "HIP_HIP" for all.
  # --------------------------------------------------------------------------
  classify_edge <- function(row_nodes, col_nodes, hip_boundary) {
    if (is.na(hip_boundary)) return(rep("HIP_HIP", length(row_nodes)))
    row_is_hip <- row_nodes <= hip_boundary
    col_is_hip <- col_nodes <= hip_boundary
    dplyr::case_when(
      row_is_hip &  col_is_hip ~ "HIP_HIP",
      !row_is_hip & !col_is_hip ~ "EHC_EHC",
      TRUE                      ~ "HIP_EHC"
    )
  }
  
  # --------------------------------------------------------------------------
  # Load all strata for this mouse
  # --------------------------------------------------------------------------
  loaded_data <- list()
  for (d_level in discrete_levels) {
    res <- load_sparse_data(d_level)
    if (!is.null(res)) {
      loaded_data[[d_level]] <- res
      message(sprintf("Loaded %s / %s", mouse_ID, d_level))
    } else {
      message(sprintf("Skipping %s / %s — file not found", mouse_ID, d_level))
    }
  }
  
  if (length(loaded_data) == 0) {
    warning(sprintf("No data found for mouse %s — returning NULL", mouse_ID))
    return(NULL)
  }
  
  found_levels <- names(loaded_data)
  
  # Pretty row labels: used as factor levels for the Row_ID facet dimension.
  # Map raw d_level -> human-readable label, preserving found_levels order.
  pretty_levels <- sapply(found_levels, pretty_label, USE.NAMES = FALSE)
  
  # ── per-stratum max_node and HIP boundary ──────────────────────────────────
  max_node_lookup <- setNames(
    sapply(found_levels, function(d) get_max_node(loaded_data[[d]]$sparse_data)),
    found_levels
  )
  hip_boundary_lookup <- setNames(
    sapply(found_levels, function(d) {
      b <- loaded_data[[d]]$boundaries
      if (length(b) == 0) NA_real_ else floor(b[1])
    }),
    found_levels
  )
  message("HIP boundary lookup: ",
          paste(found_levels, hip_boundary_lookup, sep = "=", collapse = ", "))
  
  # --------------------------------------------------------------------------
  # Build plot_data: tile centers at (2k-1)/(2*max_node), tile_size=1/max_node
  # so tile edges span exactly [0,1]. All strata share the same [0,1] coordinate
  # space so coord_fixed(ratio=1) enforces square panels identically for all.
  # Row_ID carries the pretty stratum label for facet strip display.
  # --------------------------------------------------------------------------
  plot_data_list <- list()
  
  for (d_level in found_levels) {
    max_node     <- max_node_lookup[[d_level]]
    hip_boundary <- hip_boundary_lookup[[d_level]]
    tile_size    <- 1 / max_node
    row_label    <- pretty_label(d_level)
    row_content  <- loaded_data[[d_level]]$sparse_data
    for (c_idx in seq_along(row_content)) {
      df_coords <- row_content[[c_idx]]
      if (is.null(df_coords) || nrow(df_coords) == 0) next
      original          <- df_coords
      mirrored          <- original
      mirrored$Node_Row <- original$Node_Col
      mirrored$Node_Col <- original$Node_Row
      # unique before rescaling so integer dedup is exact
      combined_df <- unique(rbind(original, mirrored))
      # assign region type before rescaling (node indices still integer here)
      if (output == "adj") {
        combined_df$Region_Type <- classify_edge(
          combined_df$Node_Row, combined_df$Node_Col, hip_boundary
        )
      }
      # tile center: (2k-1)/(2*max_node) maps node 1 to 1/(2n) and
      # node n to (2n-1)/(2n), so tile edges align exactly with [0,1]
      combined_df$Node_Row  <- (2 * combined_df$Node_Row - 1) / (2 * max_node)
      combined_df$Node_Col  <- (2 * combined_df$Node_Col - 1) / (2 * max_node)
      combined_df$tile_size <- tile_size
      combined_df$Row_ID    <- row_label
      combined_df$Col_ID    <- as.numeric(names(row_content)[c_idx])
      plot_data_list[[length(plot_data_list) + 1]] <- combined_df
    }
  }
  
  if (length(plot_data_list) > 0) {
    plot_data <- do.call(rbind, plot_data_list)
  } else {
    plot_data <- data.frame(Node_Row = NA, Node_Col = NA, tile_size = NA,
                            Row_ID = pretty_levels[1], Col_ID = display_weeks[1])
    if (output == "adj")                plot_data$Region_Type <- NA_character_
    if (output %in% c("P_HS", "C_HS")) plot_data$Value       <- NA
  }
  
  # Pad strata with no tiles (so all rows appear in facet_grid).
  # Anchor dummy Col_ID to display_weeks[1] — guaranteed to be a valid level.
  missing_labels <- setdiff(pretty_levels, unique(as.character(plot_data$Row_ID)))
  if (length(missing_labels) > 0) {
    dummy <- data.frame(Node_Row = NA, Node_Col = NA, tile_size = NA,
                        Row_ID = missing_labels, Col_ID = display_weeks[1])
    if (output == "adj")                dummy$Region_Type <- NA_character_
    if (output %in% c("P_HS", "C_HS")) dummy$Value       <- NA
    plot_data <- rbind(plot_data, dummy)
  }
  
  # Factor Col_ID over display_weeks only, with "Week N" labels.
  # This is the single source of truth for which columns appear.
  week_labels      <- setNames(paste("Week", display_weeks), display_weeks)
  plot_data$Row_ID <- factor(plot_data$Row_ID, levels = pretty_levels)
  plot_data$Col_ID <- factor(plot_data$Col_ID, levels = display_weeks)
  if (output == "adj") {
    plot_data$Region_Type <- factor(plot_data$Region_Type,
                                    levels = c("HIP_HIP", "EHC_EHC", "HIP_EHC"))
  }
  
  # --------------------------------------------------------------------------
  # Build bg_gray_data: one row per (stratum, absent week within display_weeks).
  # geom_rect with -Inf/Inf floods the entire facet panel with gray.
  # --------------------------------------------------------------------------
  bg_gray_list <- list()
  for (d_level in found_levels) {
    absent_weeks <- loaded_data[[d_level]]$absent_weeks   # already intersected with display_weeks
    row_label    <- pretty_label(d_level)
    if (length(absent_weeks) > 0)
      bg_gray_list[[d_level]] <- data.frame(Row_ID = row_label, Col_ID = absent_weeks)
  }
  
  if (length(bg_gray_list) > 0) {
    bg_gray_data        <- do.call(rbind, bg_gray_list)
    bg_gray_data$Row_ID <- factor(bg_gray_data$Row_ID, levels = pretty_levels)
    bg_gray_data$Col_ID <- factor(bg_gray_data$Col_ID, levels = display_weeks)
  } else {
    bg_gray_data <- NULL
  }
  
  # --------------------------------------------------------------------------
  # Build per-stratum boundary data: boundaries are midpoints between regions
  # (e.g. 4.5 = between node 4 and 5). Apply same rescaling as tile centers:
  # b -> (2b-1)/(2*max_node), which places the line at the exact midpoint
  # between the two flanking tile centers in [0,1] space.
  # Col_ID is factored over display_weeks so panels align with plot_data.
  # --------------------------------------------------------------------------
  boundary_data_list <- list()
  for (d_level in found_levels) {
    boundaries_i <- loaded_data[[d_level]]$boundaries
    max_node     <- max_node_lookup[[d_level]]
    row_label    <- pretty_label(d_level)
    if (length(boundaries_i) == 0) next
    if (max_node == 1) next
    absent_i        <- loaded_data[[d_level]]$absent_weeks
    present_weeks_i <- setdiff(display_weeks, absent_i)
    if (length(present_weeks_i) == 0) next
    bd <- expand.grid(
      Row_ID   = row_label,
      Col_ID   = present_weeks_i,
      boundary = (2 * boundaries_i - 1) / (2 * max_node),
      stringsAsFactors = FALSE
    )
    boundary_data_list[[d_level]] <- bd
  }
  
  if (length(boundary_data_list) > 0) {
    boundary_data        <- do.call(rbind, boundary_data_list)
    boundary_data$Row_ID <- factor(boundary_data$Row_ID, levels = pretty_levels)
    # factor over display_weeks — must match plot_data$Col_ID levels exactly
    boundary_data$Col_ID <- factor(boundary_data$Col_ID, levels = display_weeks)
  } else {
    boundary_data <- NULL
  }
  
  # --------------------------------------------------------------------------
  # Build the plot
  # --------------------------------------------------------------------------
  new_palette <- hcl.colors(3, palette = 'Blue-Red 2')
  c_low  <- new_palette[1]
  c_mid  <- new_palette[2]
  c_high <- new_palette[3]
  
  g <- ggplot() +
    
    # Layer 1: gray background for absent weeks — floods entire panel
    { if (!is.null(bg_gray_data) && nrow(bg_gray_data) > 0)
      geom_rect(data    = bg_gray_data,
                mapping = aes(group = interaction(Row_ID, Col_ID)),
                xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                fill = "gray80", alpha = 0.8)
    }
  
  # Layer 2: edge tiles — tile_size varies per stratum so all matrices fill
  # [0,1] corner to corner regardless of adjacency matrix dimension.
  # For adj output, Region_Type drives fill color via scale_fill_manual.
  # For P_HS / C_HS output, continuous Value drives fill via scale_fill_gradient2.
  if (output == "adj") {
    g <- g + geom_tile(data = plot_data,
                       aes(x = Node_Col, y = Node_Row,
                           width = tile_size, height = tile_size,
                           fill = Region_Type)) +
      scale_fill_manual(
        name   = "Edge type",
        values = c(HIP_HIP = COL_HIP_HIP,
                   EHC_EHC = COL_EHC_EHC,
                   HIP_EHC = COL_HIP_EHC),
        labels = c(HIP_HIP = "HIP\u2013HIP",
                   EHC_EHC = "EHC\u2013EHC",
                   HIP_EHC = "HIP\u2013EHC"),
        na.value = "transparent"
      ) +
      theme(legend.position = "right")
  } else {
    zmin <- min(plot_data$Value, na.rm = TRUE)
    zmax <- max(plot_data$Value, na.rm = TRUE)
    zmax <- zmax + 0.1 * (zmax - zmin)
    zmin <- zmin - 0.1 * (zmax - zmin)
    g <- g + geom_tile(data = plot_data,
                       aes(x = Node_Col, y = Node_Row,
                           width = tile_size, height = tile_size,
                           fill = Value)) +
      scale_fill_gradient2(low = c_low, mid = c_mid, high = c_high,
                           midpoint = 0, limits = c(zmin, zmax)) +
      theme(legend.position = "right")
  }
  
  # Layer 3: per-stratum region borders in rescaled [0,1] space
  if (!is.null(boundary_data) && nrow(boundary_data) > 0) {
    g <- g +
      geom_vline(data = boundary_data,
                 aes(xintercept = boundary),
                 color = "gray40", alpha = 1, linewidth = 0.3) +
      geom_hline(data = boundary_data,
                 aes(yintercept = boundary),
                 color = "gray40", alpha = 1, linewidth = 0.3)
  }
  
  # Layer 4: facet + formatting.
  # All strata share [0,1] coordinate space so scales = "fixed" works and
  # coord_fixed(ratio=1) enforces square panels.
  # scale_x/y with expand=c(0,0) removes ggplot's default padding so tiles
  # fill the panel exactly edge to edge.
  # scale_y_reverse: node 1 at top, conventional matrix layout.
  # plot title is blank (mouse ID removed).
  # labeller: column strips show "Week N"; row strips show pretty stratum label.
  g <- g +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_reverse(limits = c(1, 0), expand = c(0, 0)) +
    facet_grid(Row_ID ~ Col_ID,
               drop     = TRUE,
               scales   = "fixed",
               labeller = labeller(Col_ID = week_labels)) +
    coord_fixed(ratio = 1) +
    ggtitle("") +
    theme_minimal(base_size = 15) +
    theme(
      axis.text        = element_blank(),
      axis.title       = element_blank(),
      axis.ticks       = element_blank(),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background  = element_rect(fill = "transparent", color = NA),
      plot.title       = element_text(face = "bold", size = rel(1.4), hjust = 0.5),
      strip.background = element_rect(fill = "gray95"),
      strip.text.x     = element_text(face = "bold", size = rel(0.9)),
      strip.text.y     = element_text(face = "bold", size = rel(1.1)),
      legend.text      = element_text(size = rel(1.1)),
      legend.title     = element_text(size = rel(1.1)),
      legend.key.size  = unit(1.2, "lines")
    )
  
  return(list(
    plot         = g,
    found_levels = found_levels
  ))
}

visualize_edge_set_one_mouse_all_strata_v2 <- function(results_folder, time_scale, discrete_levels, output, region_border, mouse_ID) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: For a single mouse, plot all available discrete strata as rows in a
  #       single facet_grid. Rows = stratum label, Columns = week.
  #       One set of axis labels. Stratum labels are human-readable
  #       (m0vr1 -> "Resting", m1vr1 -> "Running"; others passed through).
  #
  #       Week wrapping: all_weeks are split into two halves. The first half
  #       occupies wrap_group = 1 and the second half occupies wrap_group = 2.
  #       The facet is facet_grid(wrap_group + Row_ID ~ Col_ID) so that both
  #       halves have their own column and row label strips, with a vertical
  #       gap between them. This converts a 2x22 layout into a 4x11 layout
  #       (for 2 strata and 22 weeks).
  #
  # inputs:
  #
  # - results_folder   (string)
  # - time_scale       (integer)  e.g. 10
  # - discrete_levels  (vector of strings)  e.g. c("m0vr1", "m1vr1")
  # - output           (string)   'adj', 'P_HS', or 'C_HS'
  # - region_border    (boolean)
  # - mouse_ID         (string)   e.g. "Tau1", "WT2"
  #
  # returns: named list with two elements:
  #   - plot         : single ggplot object (all strata as rows, weeks wrapped)
  #   - found_levels : character vector of discrete levels actually loaded
  #
  # Color scheme for output == 'adj':
  #   HIP-HIP edges  -> emerald teal  (#00C896)
  #   EHC-EHC edges  -> vivid purple  (#CC3FFF)
  #   HIP-EHC edges  -> periwinkle    (#6694CC)  [sqrt-mean-squares RGB fusion]
  #
  # Region assignment: nodes 1..boundary are HIP, nodes (boundary+1)..max_node
  # are EHC, where boundary = floor(first element of boundaries vector).
  # If no boundary is available, all nodes are treated as HIP-HIP.
  #
  # ----------------------------------------------------------------------------
  
  all_weeks <- 17:38
  n_weeks   <- length(all_weeks)
  half      <- ceiling(n_weeks / 2)
  
  # Week -> wrap_group mapping (1 = first half, 2 = second half)
  week_wrap_group <- setNames(
    ifelse(seq_along(all_weeks) <= half, 1L, 2L),
    as.character(all_weeks)
  )
  # Within each wrap_group, weeks are re-labelled by their position in that half
  # so Col_ID factor levels reset to 1..half for both groups
  week_col_id <- setNames(
    ifelse(seq_along(all_weeks) <= half,
           all_weeks[seq_len(half)],
           c(all_weeks[(half + 1):n_weeks], rep(NA, half - (n_weeks - half)))[seq_len(n_weeks - half)]),
    as.character(all_weeks)
  )
  # Simpler: Col_ID is just the actual week; wrap_group partitions rows.
  # facet_grid(wrap_group + Row_ID ~ Col_ID, drop=FALSE) will show only the
  # weeks present in each wrap_group because of the nesting.
  first_half_weeks  <- all_weeks[seq_len(half)]
  second_half_weeks <- all_weeks[(half + 1):n_weeks]
  
  # Color constants for adj output
  COL_HIP_HIP <- "#00C896"   # emerald teal   (HIP)
  COL_EHC_EHC <- "#CC3FFF"   # vivid purple   (EHC)
  COL_HIP_EHC <- "#6694CC"   # periwinkle     (HIP-EHC sqrt-mean-squares fusion)
  
  # Human-readable stratum label lookup
  stratum_labels <- c(
    m0vr1 = "Resting",
    m1vr1 = "Running"
  )
  pretty_label <- function(d_level) {
    if (d_level %in% names(stratum_labels)) stratum_labels[[d_level]] else d_level
  }
  
  # --------------------------------------------------------------------------
  # Load sparse data for one mouse x one discrete level
  # --------------------------------------------------------------------------
  load_sparse_data <- function(d_level) {
    file_path <- paste0(results_folder, '/', mouse_ID, '_', d_level, '_t', time_scale, '.RData')
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      return(NULL)
    }
    tmp_env <- new.env()
    load(file_path, envir = tmp_env)
    res_i <- tmp_env$graph_results_i
    if (region_border) {
      boundaries <- get_factor_boundaries(res_i$recovery_params$kept_neuron_regions)
    } else {
      boundaries <- numeric(0)
    }
    y_c_weeks    <- as.numeric(res_i$y_c_query)
    absent_weeks <- setdiff(all_weeks, y_c_weeks)
    sparse_data  <- list()
    if (output == 'adj') {
      step_data <- res_i$step_12b
      for (idx in seq_along(y_c_weeks)) {
        current_week <- y_c_weeks[idx]
        adj_mat <- step_data[[idx]][["adj_mat_KL_GIC_local_est_eig3"]]
        if (is.null(adj_mat)) next
        coords <- which(adj_mat == 1, arr.ind = TRUE)
        if (nrow(coords) > 0) coords <- coords[coords[, 2] > coords[, 1], , drop = FALSE]
        if (nrow(coords) > 0) {
          df_coords <- as.data.frame(coords)
          colnames(df_coords) <- c("Node_Row", "Node_Col")
          sparse_data[[as.character(current_week)]] <- df_coords
        } else {
          sparse_data[[as.character(current_week)]] <- data.frame(Node_Row = integer(0),
                                                                  Node_Col = integer(0))
        }
      }
    } else if (output == 'P_HS') {
      step_data <- res_i$step_11
      for (idx in seq_along(y_c_weeks)) {
        current_week  <- y_c_weeks[idx]
        P_HS_log      <- log(step_data[[idx]][["w_mat_KL_est_eig3"]])
        if (is.null(P_HS_log)) next
        upper_tri_idx <- which(col(P_HS_log) >= row(P_HS_log), arr.ind = TRUE)
        df_coords     <- as.data.frame(upper_tri_idx)
        colnames(df_coords) <- c("Node_Row", "Node_Col")
        df_coords$Value <- P_HS_log[upper_tri_idx]
        sparse_data[[as.character(current_week)]] <- df_coords
      }
    } else if (output == 'C_HS') {
      step_data <- res_i$step_11b
      for (idx in seq_along(y_c_weeks)) {
        current_week  <- y_c_weeks[idx]
        C_HS_log      <- log(step_data[[idx]][["C_HS_KL_est_eig3"]])
        if (is.null(C_HS_log)) next
        upper_tri_idx <- which(col(C_HS_log) >= row(C_HS_log), arr.ind = TRUE)
        df_coords     <- as.data.frame(upper_tri_idx)
        colnames(df_coords) <- c("Node_Row", "Node_Col")
        df_coords$Value <- C_HS_log[upper_tri_idx]
        sparse_data[[as.character(current_week)]] <- df_coords
      }
    } else {
      stop('error in visualize_edge_set_one_mouse_all_strata: unknown output type')
    }
    rm(tmp_env, res_i, step_data)
    return(list(sparse_data = sparse_data, absent_weeks = absent_weeks, boundaries = boundaries))
  }
  
  # --------------------------------------------------------------------------
  # Helper: get max node index for a mouse from its sparse data
  # --------------------------------------------------------------------------
  get_max_node <- function(sparse_data) {
    max_node <- 0
    for (wk in sparse_data) {
      if (nrow(wk) > 0) max_node <- max(max_node, max(wk$Node_Row, wk$Node_Col))
    }
    max(max_node, 1)
  }
  
  # --------------------------------------------------------------------------
  # Helper: classify an edge (row_node, col_node) given the HIP boundary.
  # All nodes <= hip_boundary are HIP; nodes > hip_boundary are EHC.
  # Returns one of: "HIP_HIP", "EHC_EHC", "HIP_EHC"
  # If hip_boundary is NA (no boundary available), returns "HIP_HIP" for all.
  # --------------------------------------------------------------------------
  classify_edge <- function(row_nodes, col_nodes, hip_boundary) {
    if (is.na(hip_boundary)) return(rep("HIP_HIP", length(row_nodes)))
    row_is_hip <- row_nodes <= hip_boundary
    col_is_hip <- col_nodes <= hip_boundary
    dplyr::case_when(
      row_is_hip &  col_is_hip ~ "HIP_HIP",
      !row_is_hip & !col_is_hip ~ "EHC_EHC",
      TRUE                      ~ "HIP_EHC"
    )
  }
  
  # --------------------------------------------------------------------------
  # Load all strata for this mouse
  # --------------------------------------------------------------------------
  loaded_data <- list()
  for (d_level in discrete_levels) {
    res <- load_sparse_data(d_level)
    if (!is.null(res)) {
      loaded_data[[d_level]] <- res
      message(sprintf("Loaded %s / %s", mouse_ID, d_level))
    } else {
      message(sprintf("Skipping %s / %s — file not found", mouse_ID, d_level))
    }
  }
  
  if (length(loaded_data) == 0) {
    warning(sprintf("No data found for mouse %s — returning NULL", mouse_ID))
    return(NULL)
  }
  
  found_levels <- names(loaded_data)
  
  # Pretty row labels: used as factor levels for the Row_ID facet dimension.
  # Map raw d_level -> human-readable label, preserving found_levels order.
  pretty_levels <- sapply(found_levels, pretty_label, USE.NAMES = FALSE)
  
  # ── per-stratum max_node and HIP boundary ──────────────────────────────────
  max_node_lookup <- setNames(
    sapply(found_levels, function(d) get_max_node(loaded_data[[d]]$sparse_data)),
    found_levels
  )
  hip_boundary_lookup <- setNames(
    sapply(found_levels, function(d) {
      b <- loaded_data[[d]]$boundaries
      if (length(b) == 0) NA_real_ else floor(b[1])
    }),
    found_levels
  )
  message("HIP boundary lookup: ",
          paste(found_levels, hip_boundary_lookup, sep = "=", collapse = ", "))
  
  # --------------------------------------------------------------------------
  # Build plot_data: tile centers at (2k-1)/(2*max_node), tile_size=1/max_node
  # so tile edges span exactly [0,1]. All strata share the same [0,1] coordinate
  # space so coord_fixed(ratio=1) enforces square panels identically for all.
  # Row_ID carries the pretty stratum label for facet strip display.
  # Wrap_Group assigns each week to either the first or second panel half.
  # --------------------------------------------------------------------------
  plot_data_list <- list()
  
  for (d_level in found_levels) {
    max_node     <- max_node_lookup[[d_level]]
    hip_boundary <- hip_boundary_lookup[[d_level]]
    tile_size    <- 1 / max_node
    row_label    <- pretty_label(d_level)
    row_content  <- loaded_data[[d_level]]$sparse_data
    for (c_idx in seq_along(row_content)) {
      df_coords <- row_content[[c_idx]]
      if (is.null(df_coords) || nrow(df_coords) == 0) next
      original          <- df_coords
      mirrored          <- original
      mirrored$Node_Row <- original$Node_Col
      mirrored$Node_Col <- original$Node_Row
      # unique before rescaling so integer dedup is exact
      combined_df <- unique(rbind(original, mirrored))
      # assign region type before rescaling (node indices still integer here)
      if (output == "adj") {
        combined_df$Region_Type <- classify_edge(
          combined_df$Node_Row, combined_df$Node_Col, hip_boundary
        )
      }
      # tile center: (2k-1)/(2*max_node) maps node 1 to 1/(2n) and
      # node n to (2n-1)/(2n), so tile edges align exactly with [0,1]
      combined_df$Node_Row  <- (2 * combined_df$Node_Row - 1) / (2 * max_node)
      combined_df$Node_Col  <- (2 * combined_df$Node_Col - 1) / (2 * max_node)
      combined_df$tile_size <- tile_size
      combined_df$Row_ID    <- row_label
      current_week          <- as.numeric(names(row_content)[c_idx])
      combined_df$Col_ID    <- current_week
      combined_df$Wrap_Group <- week_wrap_group[as.character(current_week)]
      plot_data_list[[length(plot_data_list) + 1]] <- combined_df
    }
  }
  
  # All (Wrap_Group, Row_ID, Col_ID) combinations needed for drop=FALSE padding
  # Col_ID factor: first half uses first_half_weeks, second half uses second_half_weeks
  # We keep Col_ID as the actual week number; facet_grid nesting with Wrap_Group
  # naturally restricts each half-block to its own weeks.
  all_col_levels <- all_weeks   # full set; panels for wrong half will be empty/dropped
  
  if (length(plot_data_list) > 0) {
    plot_data <- do.call(rbind, plot_data_list)
  } else {
    plot_data <- data.frame(Node_Row = NA, Node_Col = NA, tile_size = NA,
                            Row_ID = pretty_levels[1], Col_ID = all_weeks[1],
                            Wrap_Group = 1L)
    if (output == "adj")                plot_data$Region_Type <- NA_character_
    if (output %in% c("P_HS", "C_HS")) plot_data$Value       <- NA
  }
  
  # Pad strata with no tiles (so all rows appear in facet_grid)
  missing_labels <- setdiff(pretty_levels, unique(as.character(plot_data$Row_ID)))
  if (length(missing_labels) > 0) {
    dummy <- data.frame(Node_Row = NA, Node_Col = NA, tile_size = NA,
                        Row_ID = missing_labels, Col_ID = all_weeks[1],
                        Wrap_Group = 1L)
    if (output == "adj")                dummy$Region_Type <- NA_character_
    if (output %in% c("P_HS", "C_HS")) dummy$Value       <- NA
    plot_data <- rbind(plot_data, dummy)
  }
  
  plot_data$Row_ID     <- factor(plot_data$Row_ID,     levels = pretty_levels)
  plot_data$Col_ID     <- factor(plot_data$Col_ID,     levels = all_weeks)
  plot_data$Wrap_Group <- factor(plot_data$Wrap_Group, levels = c(1L, 2L))
  if (output == "adj") {
    plot_data$Region_Type <- factor(plot_data$Region_Type,
                                    levels = c("HIP_HIP", "EHC_EHC", "HIP_EHC"))
  }
  
  # --------------------------------------------------------------------------
  # Build bg_gray_data: one row per (stratum, absent week) — geom_rect with
  # -Inf/Inf floods the entire facet panel with gray for that cell
  # --------------------------------------------------------------------------
  bg_gray_list <- list()
  for (d_level in found_levels) {
    absent_weeks <- loaded_data[[d_level]]$absent_weeks
    row_label    <- pretty_label(d_level)
    if (length(absent_weeks) > 0) {
      bg_df <- data.frame(
        Row_ID     = row_label,
        Col_ID     = absent_weeks,
        Wrap_Group = factor(week_wrap_group[as.character(absent_weeks)], levels = c(1L, 2L))
      )
      bg_gray_list[[d_level]] <- bg_df
    }
  }
  
  if (length(bg_gray_list) > 0) {
    bg_gray_data            <- do.call(rbind, bg_gray_list)
    bg_gray_data$Row_ID     <- factor(bg_gray_data$Row_ID,     levels = pretty_levels)
    bg_gray_data$Col_ID     <- factor(bg_gray_data$Col_ID,     levels = all_weeks)
  } else {
    bg_gray_data <- NULL
  }
  
  # --------------------------------------------------------------------------
  # Build per-stratum boundary data: boundaries are midpoints between regions
  # (e.g. 4.5 = between node 4 and 5). Apply same rescaling as tile centers:
  # b -> (2b-1)/(2*max_node), which places the line at the exact midpoint
  # between the two flanking tile centers in [0,1] space.
  # --------------------------------------------------------------------------
  boundary_data_list <- list()
  for (d_level in found_levels) {
    boundaries_i <- loaded_data[[d_level]]$boundaries
    max_node     <- max_node_lookup[[d_level]]
    row_label    <- pretty_label(d_level)
    if (length(boundaries_i) == 0) next
    if (max_node == 1) next
    absent_i        <- loaded_data[[d_level]]$absent_weeks
    present_weeks_i <- setdiff(all_weeks, absent_i)
    bd <- expand.grid(
      Row_ID   = row_label,
      Col_ID   = present_weeks_i,
      boundary = (2 * boundaries_i - 1) / (2 * max_node),
      stringsAsFactors = FALSE
    )
    bd$Wrap_Group <- factor(week_wrap_group[as.character(bd$Col_ID)], levels = c(1L, 2L))
    boundary_data_list[[d_level]] <- bd
  }
  
  if (length(boundary_data_list) > 0) {
    boundary_data            <- do.call(rbind, boundary_data_list)
    boundary_data$Row_ID     <- factor(boundary_data$Row_ID, levels = pretty_levels)
    boundary_data$Col_ID     <- factor(boundary_data$Col_ID, levels = all_weeks)
  } else {
    boundary_data <- NULL
  }
  
  # --------------------------------------------------------------------------
  # Custom labeller: suppress the Wrap_Group numeric label in the row strip,
  # show only Row_ID (pretty stratum label).
  # --------------------------------------------------------------------------
  wrap_labeller <- labeller(
    Wrap_Group = function(x) rep("", length(x)),
    Row_ID     = label_value
  )
  
  # --------------------------------------------------------------------------
  # Build the plot
  # --------------------------------------------------------------------------
  new_palette <- hcl.colors(3, palette = 'Blue-Red 2')
  c_low  <- new_palette[1]
  c_mid  <- new_palette[2]
  c_high <- new_palette[3]
  
  # Plot title: strip mouse ID of any trailing digits so e.g. "Tau2" -> "Tau"
  # (more precisely: remove trailing digit group from the ID string)
  plot_title <- sub("\\d+$", "", mouse_ID)
  
  g <- ggplot() +
    
    # Layer 1: gray background for absent weeks — floods entire panel
    { if (!is.null(bg_gray_data) && nrow(bg_gray_data) > 0)
      geom_rect(data    = bg_gray_data,
                mapping = aes(group = interaction(Row_ID, Col_ID)),
                xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                fill = "gray80", alpha = 0.8)
    }
  
  # Layer 2: edge tiles — tile_size varies per stratum so all matrices fill
  # [0,1] corner to corner regardless of adjacency matrix dimension.
  # For adj output, Region_Type drives fill color via scale_fill_manual.
  # For P_HS / C_HS output, continuous Value drives fill via scale_fill_gradient2.
  if (output == "adj") {
    g <- g + geom_tile(data = plot_data,
                       aes(x = Node_Col, y = Node_Row,
                           width = tile_size, height = tile_size,
                           fill = Region_Type)) +
      scale_fill_manual(
        name   = "Edge type",
        values = c(HIP_HIP = COL_HIP_HIP,
                   EHC_EHC = COL_EHC_EHC,
                   HIP_EHC = COL_HIP_EHC),
        labels = c(HIP_HIP = "HIP\u2013HIP",
                   EHC_EHC = "EHC\u2013EHC",
                   HIP_EHC = "HIP\u2013EHC"),
        na.value = "transparent"
      ) +
      guides(fill = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  } else {
    zmin <- min(plot_data$Value, na.rm = TRUE)
    zmax <- max(plot_data$Value, na.rm = TRUE)
    zmax <- zmax + 0.1 * (zmax - zmin)
    zmin <- zmin - 0.1 * (zmax - zmin)
    g <- g + geom_tile(data = plot_data,
                       aes(x = Node_Col, y = Node_Row,
                           width = tile_size, height = tile_size,
                           fill = Value)) +
      scale_fill_gradient2(low = c_low, mid = c_mid, high = c_high,
                           midpoint = 0, limits = c(zmin, zmax)) +
      theme(legend.position = "bottom")
  }
  
  # Layer 3: per-stratum region borders in rescaled [0,1] space
  if (!is.null(boundary_data) && nrow(boundary_data) > 0) {
    g <- g +
      geom_vline(data = boundary_data,
                 aes(xintercept = boundary),
                 color = "gray40", alpha = 1, linewidth = 0.3) +
      geom_hline(data = boundary_data,
                 aes(yintercept = boundary),
                 color = "gray40", alpha = 1, linewidth = 0.3)
  }
  
  # Layer 4: facet + formatting.
  # facet_grid(Wrap_Group + Row_ID ~ Col_ID): nesting Wrap_Group with Row_ID
  # splits weeks into two stacked blocks, each with its own column and row
  # strips. panel.spacing.y controls the gap between the two week-halves.
  # wrap_labeller suppresses the numeric Wrap_Group strip, showing only Row_ID.
  g <- g +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_reverse(limits = c(1, 0), expand = c(0, 0)) +
    facet_grid(Wrap_Group + Row_ID ~ Col_ID,
               drop      = TRUE,
               scales    = "fixed",
               labeller  = wrap_labeller) +
    coord_fixed(ratio = 1) +
    ggtitle(plot_title) +
    theme_minimal(base_size = 15) +
    theme(
      axis.text        = element_blank(),
      axis.title       = element_blank(),
      axis.ticks       = element_blank(),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background  = element_rect(fill = "transparent", color = NA),
      plot.title       = element_text(face = "bold", size = rel(1.4), hjust = 0.5),
      strip.background = element_rect(fill = "gray95"),
      strip.text       = element_text(face = "bold", size = rel(2)),
      # gap between the two wrap-group blocks
      panel.spacing.y  = unit(0.8, "lines"),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold"),
      legend.key.width = unit(1.5, "lines")
    )
  
  return(list(
    plot         = g,
    found_levels = found_levels
  ))
}