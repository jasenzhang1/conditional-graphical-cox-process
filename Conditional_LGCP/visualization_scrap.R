library(ggplot2)

# ------------------------------------------------------------
# Toy dataset
# ------------------------------------------------------------
set.seed(42)

make_adj_edges <- function(n, density = 0.3) {
  mat <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (runif(1) < density) {
        mat[i, j] <- 1
        mat[j, i] <- 1
      }
    }
  }
  coords <- which(mat == 1 & col(mat) > row(mat), arr.ind = TRUE)
  if (nrow(coords) == 0) return(data.frame(Node_Row = integer(0), Node_Col = integer(0)))
  df <- as.data.frame(coords)
  colnames(df) <- c("Node_Row", "Node_Col")
  df
}

# Each mouse: list of weeks -> edge data frames
toy_data <- list(
  Tau1 = list(
    sparse_data  = list("1" = make_adj_edges(2),
                        "2" = make_adj_edges(2),
                        "3" = make_adj_edges(2),
                        "4" = make_adj_edges(2)),
    absent_weeks = integer(0),
    boundaries   = c(2.5)        # one boundary after node 2 (only 2 nodes so at the edge)
  ),
  Tau2 = list(
    sparse_data  = list("1" = make_adj_edges(5),
                        "2" = make_adj_edges(5),
                        "3" = make_adj_edges(5),
                        "4" = make_adj_edges(5)),
    absent_weeks = integer(0),
    boundaries   = c(3.5)        # one boundary splitting nodes 1-3 | 4-5
  ),
  WT1 = list(
    sparse_data  = list("1" = make_adj_edges(10),
                        "2" = make_adj_edges(10),
                        "3" = make_adj_edges(10),
                        "4" = make_adj_edges(10)),
    absent_weeks = integer(0),
    boundaries   = c(4.5, 7.5)    # two boundaries: 1-4 | 5-7 | 8-10
  ),
  WT2 = list(
    sparse_data  = list("1" = make_adj_edges(20),
                        "2" = make_adj_edges(20),
                        "3" = make_adj_edges(20),
                        "4" = make_adj_edges(20)),
    absent_weeks = integer(0),
    boundaries   = c(8.5, 14.5)   # two boundaries: 1-8 | 9-14 | 15-20
  ),
  WT3 = list(
    # weeks 1 and 3 absent; weeks 2 and 4 present but no edges (empty data frames)
    sparse_data  = list("2" = data.frame(Node_Row = integer(0), Node_Col = integer(0)),
                        "4" = data.frame(Node_Row = integer(0), Node_Col = integer(0))),
    absent_weeks = c(1, 3),
    boundaries   = numeric(0)
  )
)

# ------------------------------------------------------------
# Toy version of the plotting function
# ------------------------------------------------------------
plot_toy_adj <- function(loaded_data) {
  
  all_weeks   <- 1:4
  present_IDs <- names(loaded_data)
  
  get_max_node <- function(sparse_data) {
    max_node <- 0
    for (wk in sparse_data) {
      if (nrow(wk) > 0) max_node <- max(max_node, max(wk$Node_Row, wk$Node_Col))
    }
    max(max_node, 1)
  }
  
  max_node_lookup <- setNames(
    sapply(present_IDs, function(ID) get_max_node(loaded_data[[ID]]$sparse_data)),
    present_IDs
  )
  
  # ------------------------------------------------------------------------
  # Build plot_data: tile centers at (2k-1)/(2*max_node), tile_size=1/max_node
  # so tile edges span exactly [0,1]. Node 1 -> center at 1/(2n), 
  # node n -> center at (2n-1)/(2n). First tile left edge = 0, last tile
  # right edge = 1.
  # ------------------------------------------------------------------------
  plot_data_list <- list()
  
  for (ID in present_IDs) {
    max_node    <- max_node_lookup[[ID]]
    tile_size   <- 1 / max_node
    row_content <- loaded_data[[ID]]$sparse_data
    for (c_idx in seq_along(row_content)) {
      df_coords <- row_content[[c_idx]]
      if (is.null(df_coords) || nrow(df_coords) == 0) next
      original          <- df_coords
      mirrored          <- original
      mirrored$Node_Row <- original$Node_Col
      mirrored$Node_Col <- original$Node_Row
      # unique before rescaling so integer dedup is exact
      combined_df <- unique(rbind(original, mirrored))
      # tile center: (2k - 1) / (2 * max_node)
      combined_df$Node_Row  <- (2 * combined_df$Node_Row - 1) / (2 * max_node)
      combined_df$Node_Col  <- (2 * combined_df$Node_Col - 1) / (2 * max_node)
      combined_df$tile_size <- tile_size
      combined_df$Row_ID    <- ID
      combined_df$Col_ID    <- as.integer(names(row_content)[c_idx])
      plot_data_list[[length(plot_data_list) + 1]] <- combined_df
    }
  }
  
  plot_data        <- do.call(rbind, plot_data_list)
  plot_data$Row_ID <- factor(plot_data$Row_ID, levels = present_IDs)
  plot_data$Col_ID <- factor(plot_data$Col_ID, levels = all_weeks)
  
  # ------------------------------------------------------------------------
  # Build bg_gray_data: one row per (mouse, absent week) — geom_rect with
  # -Inf/Inf floods the entire facet panel with gray for that cell
  # ------------------------------------------------------------------------
  bg_gray_list <- list()
  for (ID in present_IDs) {
    absent_weeks <- loaded_data[[ID]]$absent_weeks
    if (length(absent_weeks) > 0)
      bg_gray_list[[ID]] <- data.frame(Row_ID = ID, Col_ID = absent_weeks)
  }
  
  if (length(bg_gray_list) > 0) {
    bg_gray_data        <- do.call(rbind, bg_gray_list)
    bg_gray_data$Row_ID <- factor(bg_gray_data$Row_ID, levels = present_IDs)
    bg_gray_data$Col_ID <- factor(bg_gray_data$Col_ID, levels = all_weeks)
  } else {
    bg_gray_data <- NULL
  }
  
  # ------------------------------------------------------------------------
  # Boundary data: boundaries are midpoints between regions (e.g. 4.5).
  # Same rescaling as tile centers: (2b - 1) / (2 * max_node)
  # ------------------------------------------------------------------------
  boundary_data_list <- list()
  for (ID in present_IDs) {
    boundaries_i <- loaded_data[[ID]]$boundaries
    max_node     <- max_node_lookup[[ID]]
    if (length(boundaries_i) == 0) next
    if (max_node == 1) next
    bd <- expand.grid(
      Row_ID   = ID,
      Col_ID   = all_weeks,
      boundary = (2 * boundaries_i - 1) / (2 * max_node),
      stringsAsFactors = FALSE
    )
    boundary_data_list[[ID]] <- bd
  }
  
  if (length(boundary_data_list) > 0) {
    boundary_data        <- do.call(rbind, boundary_data_list)
    boundary_data$Row_ID <- factor(boundary_data$Row_ID, levels = present_IDs)
    boundary_data$Col_ID <- factor(boundary_data$Col_ID, levels = all_weeks)
  } else {
    boundary_data <- NULL
  }
  
  # ------------------------------------------------------------------------
  # Build the plot — no anchor layer; panel range fixed via scale limits
  # ------------------------------------------------------------------------
  g <- ggplot() +
    
    # Layer 1: gray background for absent weeks — floods entire panel
    { if (!is.null(bg_gray_data) && nrow(bg_gray_data) > 0)
      geom_rect(data    = bg_gray_data,
                mapping = aes(group = interaction(Row_ID, Col_ID)),
                xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                fill = "gray80", alpha = 0.8)
    } +
    
    # Layer 2: edge tiles
    geom_tile(data = plot_data,
              aes(x = Node_Col, y = Node_Row,
                  width = tile_size, height = tile_size),
              fill = "red") +
    
    # Layer 3: boundary lines
    { if (!is.null(boundary_data) && nrow(boundary_data) > 0)
      list(
        geom_vline(data = boundary_data,
                   aes(xintercept = boundary),
                   color = "gray40", size = 0.5),
        geom_hline(data = boundary_data,
                   aes(yintercept = boundary),
                   color = "gray40", size = 0.5)
      )
    } +
    
    # Hard-set panel range to [0, 1] x [0, 1] with no expansion padding.
    # expand = c(0,0) ensures ggplot doesn't add extra space beyond limits.
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_reverse(limits = c(1, 0), expand = c(0, 0)) +
    
    # scales = "fixed": all panels share [0,1] coordinate space
    # coord_fixed(ratio = 1): enforces square panels
    facet_grid(Row_ID ~ Col_ID, drop = FALSE, scales = "fixed") +
    coord_fixed(ratio = 1) +
    
    theme_minimal(base_size = 15) +
    theme(
      axis.text        = element_blank(),
      axis.title       = element_blank(),
      axis.ticks       = element_blank(),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background  = element_rect(fill = "transparent", color = NA),
      strip.background = element_rect(fill = "gray95"),
      strip.text       = element_text(face = "bold", size = rel(1.5)),
      legend.position  = "none"
    )
  
  g
}
# ------------------------------------------------------------
# Run and display
# ------------------------------------------------------------
g <- plot_toy_adj(toy_data)
print(g)

png_name <- 'visualization_scrap.png'
png(png_name, width = 45, height = 10, units = "in", res = 100)
print(g)
dev.off()
