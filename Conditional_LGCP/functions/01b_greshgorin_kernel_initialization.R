# Gershgorin-Based Robust Precision Matrix Construction
# Implementation with improved bounds for high-dimensional applications

library(Matrix)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(viridis)

# =============================================================================
# CORE IMPLEMENTATION WITH GERSHGORIN BOUNDS
# =============================================================================

#' Construct precision matrix using Gershgorin-based bounds
#' @param p Matrix dimension
#' @param y_c Covariate vector
#' @param alpha_functions Named list of functions alpha_ij(y_c)
#' @param structure_type One of "banded", "exponential_decay", "sparse"
#' @param structure_params List with structure-specific parameters
#' @param epsilon Minimum absolute value threshold
#' @param delta Regularization parameter for spectral adjustment
#' @return List with precision matrix and diagnostics
construct_gershgorin_precision_matrix <- function(p, y_c, alpha_functions, 
                                                  structure_type = "exponential_decay",
                                                  structure_params = list(),
                                                  epsilon = 0.05, delta = 1e-6) {
  
  # Determine optimal rho_max based on structure type
  rho_max <- calculate_optimal_rho_max(p, structure_type, structure_params)
  
  # cat("Using structure:", structure_type, "\n")
  # cat("Optimal rho_max:", round(rho_max, 4), "\n")
  
  # Construct correlation matrix
  rho_values <- matrix(0, p, p)
  diag(rho_values) <- 1
  
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      key <- paste(i, j, sep = "_")
      if (key %in% names(alpha_functions)) {
        alpha_val <- alpha_functions[[key]](y_c)
        rho_val <- rho_max * tanh(alpha_val)
        
        # Apply structure-specific constraints
        rho_val <- apply_structure_constraints(rho_val, i, j, structure_type, structure_params)
        
        # Apply minimum threshold
        if (abs(rho_val) >= epsilon) {
          rho_values[i, j] <- rho_val
          rho_values[j, i] <- rho_val
        }
      }
    }
  }
  
  # Calculate Gershgorin row sums
  row_sums <- calculate_row_sums(rho_values)
  r_max <- max(row_sums)
  
  # cat("Gershgorin r_max:", round(r_max, 4), "\n")
  # cat("Gershgorin condition r_max < 1:", r_max < 1, "\n")
  
  # Spectral adjustment if needed
  R0 <- rho_values
  if (r_max < 1) {
    # No spectral adjustment needed
    P <- R0
    eigenvalues_original <- eigen(R0, symmetric = TRUE, only.values = TRUE)$values
    eigenvalues_adjusted <- eigenvalues_original
    spectral_adjustment_needed <- FALSE
    theoretical_min_eigenvalue <- 1 - r_max
  } else {
    # Apply spectral adjustment
    eigen_decomp <- eigen(R0, symmetric = TRUE)
    Q <- eigen_decomp$vectors
    eigenvalues_original <- eigen_decomp$values
    eigenvalues_adjusted <- pmax(eigenvalues_original, delta)
    P <- Q %*% diag(eigenvalues_adjusted) %*% t(Q)
    diag(P) <- 1  # Ensure exact unit diagonal
    spectral_adjustment_needed <- TRUE
    theoretical_min_eigenvalue <- 1 - r_max
  }
  
  return(list(
    precision_matrix = P,
    original_matrix = R0,
    eigenvalues_original = eigenvalues_original,
    eigenvalues_adjusted = eigenvalues_adjusted,
    spectral_adjustment_needed = spectral_adjustment_needed,
    min_eigenvalue_original = min(eigenvalues_original),
    min_eigenvalue_adjusted = min(eigenvalues_adjusted),
    gershgorin_r_max = r_max,
    gershgorin_bound_satisfied = r_max < 1,
    theoretical_min_eigenvalue = theoretical_min_eigenvalue,
    rho_max = rho_max,
    structure_type = structure_type,
    structure_params = structure_params,
    y_c = y_c,
    row_sums = row_sums
  ))
}

# =============================================================================
# STRUCTURE-SPECIFIC PARAMETER CALCULATION
# =============================================================================

#' Calculate optimal rho_max based on structure type and Gershgorin bounds
#' @param p Matrix dimension
#' @param structure_type Structure type
#' @param structure_params Structure parameters
calculate_optimal_rho_max <- function(p, structure_type, structure_params) {
  
  if (structure_type == "banded") {
    # Banded: rho_max < 1/(2k)
    k <- structure_params$k %||% 5  # Default bandwidth
    rho_max <- 0.9 / (2 * k)  # 90% of theoretical maximum
    
  } else if (structure_type == "exponential_decay") {
    # Exponential decay: rho_max < (e^gamma - 1)/2
    gamma <- structure_params$gamma %||% 0.5  # Default decay rate
    rho_max <- 0.9 * (exp(gamma) - 1) / 2  # 90% of theoretical maximum
    
  } else if (structure_type == "sparse") {
    # Sparse: rho_max < 1/s
    s <- structure_params$s %||% 10  # Default sparsity level
    rho_max <- 0.9 / s  # 90% of theoretical maximum
    
  } else {
    # Fallback to conservative bound
    rho_max <- 0.1
    warning("Unknown structure type. Using conservative rho_max = 0.1")
  }
  
  return(rho_max)
}

#' Apply structure-specific constraints to correlation values
#' @param rho_val Current correlation value
#' @param i Row index
#' @param j Column index
#' @param structure_type Structure type
#' @param structure_params Structure parameters
apply_structure_constraints <- function(rho_val, i, j, structure_type, structure_params) {
  
  if (structure_type == "banded") {
    k <- structure_params$k %||% 5
    if (abs(i - j) > k) {
      return(0)  # Outside band
    }
    
  } else if (structure_type == "exponential_decay") {
    gamma <- structure_params$gamma %||% 0.5
    decay_factor <- exp(-gamma * abs(i - j))
    rho_val <- rho_val * decay_factor
    
  } else if (structure_type == "sparse") {
    # Sparse structure handled in alpha functions
    # No additional constraints here
  }
  
  return(rho_val)
}

#' Calculate row sums for Gershgorin bound
#' @param rho_matrix Correlation matrix
calculate_row_sums <- function(rho_matrix) {
  p <- nrow(rho_matrix)
  row_sums <- numeric(p)
  
  for (i in 1:p) {
    row_sums[i] <- sum(abs(rho_matrix[i, -i]))  # Sum of absolute off-diagonal elements
  }
  
  return(row_sums)
}

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# =============================================================================
# SPECIALIZED ALPHA FUNCTION CREATORS
# =============================================================================

#' Create alpha functions for banded structure
#' @param p Matrix dimension
#' @param k Bandwidth parameter
#' @param covariate_strength Strength of covariate modulation
#' 
create_banded_alpha <- function(p, k = 5, covariate_strength = 0.5) {
  
  # alpha_ij = 1.0 + cs * sum(y_c)
  
  alpha_functions <- list()
  
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      if (abs(i - j) <= k) {  # Only within band
        key <- paste(i, j, sep = "_")
        
        alpha_functions[[key]] <- local({
          i_local <- i
          j_local <- j
          strength <- covariate_strength
          function(y_c) {
            base_strength <- 1.0  # Base correlation strength
            covariate_effect <- strength * sum(y_c)
            return(base_strength + covariate_effect)
          }
        })
      }
    }
  }
  
  return(alpha_functions)
}

#' Create alpha functions for exponential decay structure
#' @param p Matrix dimension
#' @param gamma Decay rate parameter
#' @param covariate_strength Strength of covariate modulation
create_exponential_decay_alpha <- function(p, gamma = 0.5, covariate_strength = 0.3) {
  alpha_functions <- list()
  
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      key <- paste(i, j, sep = "_")
      
      alpha_functions[[key]] <- local({
        i_local <- i
        j_local <- j
        gamma_param <- gamma
        strength <- covariate_strength
        function(y_c) {
          # Base decay handled in constraint application
          base_strength <- 1.0
          covariate_effect <- strength * sum(y_c)
          return(base_strength + covariate_effect)
        }
      })
    }
  }
  
  return(alpha_functions)
}

#' Create alpha functions for sparse structure
#' @param p Matrix dimension
#' @param s Maximum number of connections per node
#' @param connection_prob Probability of connection
#' @param covariate_strength Strength of covariate modulation
create_sparse_alpha <- function(p, s, connection_prob, covariate_strength = 0.4) {
  alpha_functions <- list()
  
  # Generate sparse connectivity pattern
  set.seed(123)  # For reproducibility
  connections <- 0
  
  for (i in 1:(p-1)) {
    connections_for_i <- 0
    for (j in (i+1):p) {
      if (connections_for_i < s && runif(1) < connection_prob) {
        key <- paste(i, j, sep = "_")
        connections_for_i <- connections_for_i + 1
        
        alpha_functions[[key]] <- local({
          i_local <- i
          j_local <- j
          strength <- covariate_strength
          function(y_c) {
            base_strength <- 1.0
            covariate_effect <- strength * sum(y_c)
            return(base_strength + covariate_effect)
          }
        })
      }
    }
  }
  
  cat("Created", length(alpha_functions), "sparse connections\n")
  return(alpha_functions)
}

# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

#' Visualize precision matrix as heatmap
#' @param P Precision matrix
#' @param title Plot title
#' @param show_values Whether to show correlation values as text
plot_precision_heatmap <- function(P, title = "Precision Matrix", show_values = FALSE) {
  p <- nrow(P)
  
  # Convert to long format
  df <- expand.grid(Row = 1:p, Col = 1:p)
  df$Value <- as.vector(P)
  
  # Create base plot
  plot <- ggplot(df, aes(x = Col, y = p + 1 - Row, fill = Value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, name = "Correlation") +
    labs(title = title, x = "Column", y = "Row") +
    theme_minimal() +
    theme(aspect.ratio = 1,
          axis.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5)) +
    coord_fixed()
  
  # Add text values if requested and matrix is small enough
  if (show_values && p <= 15) {
    plot <- plot + geom_text(aes(label = round(Value, 2)), size = 3)
  }
  
  return(plot)
}

#' Plot eigenvalue spectrum
#' @param eigenvalues Vector of eigenvalues
#' @param title Plot title
#' @param theoretical_bound Theoretical minimum eigenvalue bound
plot_eigenvalue_spectrum <- function(eigenvalues, title = "Eigenvalue Spectrum", 
                                     theoretical_bound = NULL) {
  df <- data.frame(
    Index = 1:length(eigenvalues),
    Eigenvalue = sort(eigenvalues, decreasing = TRUE)
  )
  
  plot <- ggplot(df, aes(x = Index, y = Eigenvalue)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_line(alpha = 0.5) +
    labs(title = title, 
         x = "Eigenvalue Index (sorted)", 
         y = "Eigenvalue") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Add theoretical bound line if provided
  if (!is.null(theoretical_bound)) {
    plot <- plot + 
      geom_hline(yintercept = theoretical_bound, color = "red", 
                 linetype = "dashed", size = 1) +
      annotate("text", x = length(eigenvalues) * 0.7, y = theoretical_bound + 0.1, 
               label = paste("Theoretical bound:", round(theoretical_bound, 3)), 
               color = "red")
  }
  
  return(plot)
}

#' Plot Gershgorin row sums
#' @param row_sums Vector of row sums
#' @param title Plot title
plot_gershgorin_bounds <- function(row_sums, title = "Gershgorin Row Sums") {
  df <- data.frame(
    Row = 1:length(row_sums),
    RowSum = row_sums
  )
  
  ggplot(df, aes(x = Row, y = RowSum)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed", size = 1) +
    labs(title = title, 
         x = "Matrix Row", 
         y = "Sum of |off-diagonal elements|") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate("text", x = length(row_sums) * 0.8, y = 1.05, 
             label = "Gershgorin bound (r_max < 1)", color = "red")
}

#' Create comprehensive visualization for precision matrix results
#' @param result Output from construct_gershgorin_precision_matrix
#' @param show_matrix_values Whether to show values in heatmap
create_comprehensive_plots <- function(result, show_matrix_values = FALSE) {
  
  # Matrix heatmap
  p1 <- plot_precision_heatmap(result$precision_matrix, 
                               paste(result$structure_type, "Precision Matrix"),
                               show_matrix_values)
  
  # Eigenvalue spectrum
  p2 <- plot_eigenvalue_spectrum(result$eigenvalues_adjusted,
                                 paste(result$structure_type, "Eigenvalue Spectrum"),
                                 result$theoretical_min_eigenvalue)
  
  # Gershgorin bounds
  p3 <- plot_gershgorin_bounds(result$row_sums,
                               paste(result$structure_type, "Gershgorin Bounds"))
  
  # Sparsity pattern (for visualization of structure)
  sparsity_matrix <- abs(result$precision_matrix) > 1e-10
  sparsity_df <- expand.grid(Row = 1:nrow(sparsity_matrix), Col = 1:ncol(sparsity_matrix))
  sparsity_df$NonZero <- as.vector(sparsity_matrix)
  
  p4 <- ggplot(sparsity_df, aes(x = Col, y = nrow(sparsity_matrix) + 1 - Row, fill = NonZero)) +
    geom_tile() +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "black"), 
                      name = "Non-zero") +
    labs(title = paste(result$structure_type, "Sparsity Pattern"), 
         x = "Column", y = "Row") +
    theme_minimal() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    coord_fixed()
  
  return(list(heatmap = p1, eigenvalues = p2, gershgorin = p3, sparsity = p4))
}

# =============================================================================
# VALIDATION AND TESTING
# =============================================================================

#' Validate Gershgorin-based constraints
#' @param result Output from construct_gershgorin_precision_matrix
#' @param epsilon Minimum threshold parameter
validate_gershgorin_constraints <- function(result, epsilon) {
  P <- result$precision_matrix
  p <- nrow(P)
  
  constraints <- list()
  
  # Constraint 1: Gershgorin bound satisfied
  constraints$gershgorin_bound <- list(
    satisfied = result$gershgorin_bound_satisfied,
    value = result$gershgorin_r_max,
    description = "Gershgorin r_max < 1"
  )
  
  # Constraint 2: Positive semidefinite
  min_eig <- min(result$eigenvalues_adjusted)
  constraints$positive_semidefinite <- list(
    satisfied = min_eig >= -1e-10,
    value = min_eig,
    description = "Positive semidefinite"
  )
  
  # Constraint 3: Theoretical vs actual minimum eigenvalue
  theoretical_gap <- abs(result$min_eigenvalue_original - result$theoretical_min_eigenvalue)
  constraints$theoretical_bound <- list(
    satisfied = theoretical_gap < 0.1,  # Allow some numerical error
    value = paste("Theory:", round(result$theoretical_min_eigenvalue, 4), 
                  "Actual:", round(result$min_eigenvalue_original, 4)),
    description = "Theoretical bound accuracy"
  )
  
  # Standard constraints
  diagonal_values <- diag(P)
  diagonal_ok <- all(abs(diagonal_values - 1) < 1e-10)
  constraints$unit_diagonal <- list(
    satisfied = diagonal_ok,
    value = paste(round(range(diagonal_values), 6), collapse = " to "),
    description = "Unit diagonal elements"
  )
  
  off_diag <- P[upper.tri(P) | lower.tri(P)]
  bounds_ok <- all(off_diag >= -1 - 1e-10 & off_diag <= 1 + 1e-10)
  constraints$bounded_elements <- list(
    satisfied = bounds_ok,
    value = paste(round(range(off_diag), 4), collapse = " to "),
    description = "Off-diagonal elements ∈ [-1, 1]"
  )
  
  nonzero_off_diag <- off_diag[abs(off_diag) > 1e-10]
  if (length(nonzero_off_diag) > 0) {
    min_nonzero <- min(abs(nonzero_off_diag))
    threshold_ok <- min_nonzero >= epsilon - 1e-10
  } else {
    min_nonzero <- NA
    threshold_ok <- TRUE
  }
  
  constraints$minimum_threshold <- list(
    satisfied = threshold_ok,
    value = ifelse(is.na(min_nonzero), "No nonzero elements", round(min_nonzero, 4)),
    description = paste("Minimum |value| ≥", epsilon)
  )
  
  all_satisfied <- all(sapply(constraints, function(x) x$satisfied))
  
  return(list(
    all_satisfied = all_satisfied,
    individual_constraints = constraints,
    summary = paste("Satisfied", sum(sapply(constraints, function(x) x$satisfied)), 
                    "out of", length(constraints), "constraints")
  ))
}

# =============================================================================
# DEMONSTRATION AND TESTING WITH VISUALIZATION
# =============================================================================

# -10 = 0
# -7.5 = 0.055
# -5.2 = 0.1004
# -1.95 = 0.015
# 4.2 = 0.20

# p <- 10  # Higher dimension to show improvement
# y_c <- c(4.2)
# k <- 2
# cs <- 0.1
# epsilon <- 0.05
# 
# # Banded structure with k=3
# banded_params <- list(k = k)
# alpha_funcs_banded <- create_banded_alpha(p, k = k, covariate_strength = cs)
# 
# result_banded <- construct_gershgorin_precision_matrix(
#   p, y_c, alpha_funcs_banded,
#   structure_type = "banded",
#   structure_params = banded_params,
#   epsilon = epsilon
# )
# 
# result_banded$precision_matrix
# 
# eigen(result_banded$precision_matrix)$values

# cat("=== Gershgorin-Based Robust Precision Matrix Construction ===\n\n")
# 
# # Example 1: Banded Structure
# cat("Example 1: Banded Structure\n")
# cat("---------------------------\n")
# 
# p <- 300  # Higher dimension to show improvement
# y_c <- c(0.3, -0.2)
# 
# # Banded structure with k=3
# banded_params <- list(k = 3)
# alpha_funcs_banded <- create_banded_alpha(p, k = 3, covariate_strength = 0.4)
# 
# result_banded <- construct_gershgorin_precision_matrix(
#   p, y_c, alpha_funcs_banded, 
#   structure_type = "banded",
#   structure_params = banded_params,
#   epsilon = 0.05
# )
# 
# validation_banded <- validate_gershgorin_constraints(result_banded, epsilon = 0.05)
# cat("Validation:", validation_banded$summary, "\n")
# cat("Eigenvalue range:", round(range(result_banded$eigenvalues_adjusted), 4), "\n\n")
# 
# # Create visualizations for banded structure
# plots_banded <- create_comprehensive_plots(result_banded, show_matrix_values = TRUE)
# 
# cat("Example 2: Exponential Decay Structure\n")
# cat("--------------------------------------\n")
# 
# p <- 300  # Medium dimension for visualization
# y_c <- c(0.5)
# 
# # Exponential decay with gamma=0.7
# exp_params <- list(gamma = 0.7)
# alpha_funcs_exp <- create_exponential_decay_alpha(p, gamma = 0.7, covariate_strength = 0.3)
# 
# result_exp <- construct_gershgorin_precision_matrix(
#   p, y_c, alpha_funcs_exp,
#   structure_type = "exponential_decay", 
#   structure_params = exp_params,
#   epsilon = 0.03
# )
# 
# validation_exp <- validate_gershgorin_constraints(result_exp, epsilon = 0.03)
# cat("Validation:", validation_exp$summary, "\n")
# cat("Eigenvalue range:", round(range(result_exp$eigenvalues_adjusted), 4), "\n\n")
# 
# # Create visualizations for exponential decay
# plots_exp <- create_comprehensive_plots(result_exp)
# 
# cat("Example 3: Sparse Structure\n")
# cat("---------------------------\n")
# 
# p <- 300
# y_c <- c(0.2, 0.1)
# 
# # Sparse structure with s=10 connections per node
# sparse_params <- list(s = 10)
# alpha_funcs_sparse <- create_sparse_alpha(p, s = 10, connection_prob = 0.01, covariate_strength = 0.5)
# 
# result_sparse <- construct_gershgorin_precision_matrix(
#   p, y_c, alpha_funcs_sparse,
#   structure_type = "sparse",
#   structure_params = sparse_params, 
#   epsilon = 0.04
# )
# 
# validation_sparse <- validate_gershgorin_constraints(result_sparse, epsilon = 0.04)
# cat("Validation:", validation_sparse$summary, "\n")
# cat("Eigenvalue range:", round(range(result_sparse$eigenvalues_adjusted), 4), "\n\n")
# 
# # Create visualizations for sparse structure
# plots_sparse <- create_comprehensive_plots(result_sparse)
# 
# # Display plots
# cat("=== VISUALIZATIONS ===\n")
# cat("Displaying comprehensive plots for each structure type...\n\n")
# 
# # Display banded structure plots
# print("Banded Structure Visualizations:")
# grid.arrange(plots_banded$heatmap, plots_banded$sparsity, 
#              plots_banded$eigenvalues, plots_banded$gershgorin, 
#              ncol = 2, top = "Banded Structure (k=3)")
# 
# # Display exponential decay plots  
# print("Exponential Decay Structure Visualizations:")
# grid.arrange(plots_exp$heatmap, plots_exp$sparsity,
#              plots_exp$eigenvalues, plots_exp$gershgorin,
#              ncol = 2, top = "Exponential Decay Structure (γ=0.7)")
# 
# # Display sparse structure plots
# print("Sparse Structure Visualizations:")
# grid.arrange(plots_sparse$heatmap, plots_sparse$sparsity,
#              plots_sparse$eigenvalues, plots_sparse$gershgorin, 
#              ncol = 2, top = "Sparse Structure (s=8)")
# 
# # Summary comparison
# cat("=== SUMMARY COMPARISON ===\n")
# structures <- c("Banded", "Exponential Decay", "Sparse")
# results <- list(result_banded, result_exp, result_sparse)
# 
# summary_df <- data.frame(
#   Structure = structures,
#   Dimension = sapply(results, function(x) nrow(x$precision_matrix)),
#   rho_max = sapply(results, function(x) round(x$rho_max, 4)),
#   r_max = sapply(results, function(x) round(x$gershgorin_r_max, 4)),
#   Min_Eigenvalue = sapply(results, function(x) round(x$min_eigenvalue_adjusted, 4)),
#   Spectral_Adjustment = sapply(results, function(x) x$spectral_adjustment_needed)
# )
# 
# print(summary_df)