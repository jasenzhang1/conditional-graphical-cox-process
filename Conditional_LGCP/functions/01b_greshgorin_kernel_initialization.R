# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x


calculate_optimal_rho_max <- function(p, structure_type, structure_params) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Calculate optimal rho_max based on structure type and Gershgorin bounds
  # 
  #
  # inputs:
  #
  # - p                 (integer)   Matrix dimension
  # - structure_type    (string)    Structure type, one of "banded", 'exponential_decay', or 'sparse'
  # - structure_params  (list)      Structure parameters  
  #
  # outputs:
  #
  # - rho_max           (number)    optimal rho_max by Gershgorin bounds
  # 
  # ----------------------------------------------------------------------------
  
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

create_sparse_alpha <- function(p, s, connection_prob, base_strength, covariate_strength) {
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: Create alpha functions for sparse structure
  #
  #
  # inputs:
  # 
  # - p                    (integer)       Matrix dimension
  # - s                    (integer)       Maximum number of connections per node
  # - connection_prob      (number)        Probability of connection
  # - base_strength        (number)        Intercept coefficient
  # - covariate_strength   (number)        Strength of covariate modulation
  #
  # outputs:
  #
  # - alpha_functions      (list of alpha functions)
  #
  # ----------------------------------------------------------------------------
  
  
  alpha_functions <- list()
  
  # Generate sparse connectivity pattern
  set.seed(123)  # For reproducibility
  connections <- rep(0, p)
  
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      if (connections[i] < s && connections[j] < s && runif(1) < connection_prob) {
        key <- paste(i, j, sep = "_")
        connections[i] <- connections[i] + 1
        connections[j] <- connections[j] + 1
        
        alpha_functions[[key]] <- function(y_c){
          base_strength <- base_strength
          covariate_effect <- covariate_strength * sum(y_c)
          return(base_strength + covariate_effect)
        }
      }
    }
  }
  
  cat("Created", length(alpha_functions), "sparse connections\n")
  return(alpha_functions)
}


calculate_row_sums <- function(rho_matrix) {
  
  
  # ----------------------------------------------------------------------------
  #
  # Calculate row sums for Gershgorin bound
  # 
  # inputs:
  # 
  # - rho_matrix  (p x p matrix)
  #
  # outputs:
  #
  # row_sums     (p-dim vector)   rowsum values
  #
  # ----------------------------------------------------------------------------
  
  p <- nrow(rho_matrix)
  row_sums <- numeric(p)
  
  for (i in 1:p) {
    row_sums[i] <- sum(abs(rho_matrix[i, -i]))  # Sum of absolute off-diagonal elements
  }
  
  return(row_sums)
}

apply_structure_constraints <- function(rho_val, i, j, structure_type, structure_params) {
  
  #
  # GOAL: after calculating rho = rho_max * tanh(alpha_val), we may need more massaging
  #
  #' Apply structure-specific constraints to correlation values
  #' @param rho_val Current correlation value
  #' @param i Row index
  #' @param j Column index
  #' @param structure_type Structure type
  #' @param structure_params Structure parameters  
  
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

# =============================================================================
# CORE IMPLEMENTATION WITH GERSHGORIN BOUNDS
# =============================================================================


construct_gershgorin_precision_matrix <- function(p, y_c, alpha_functions, 
                                                  structure_type = "exponential_decay",
                                                  structure_params = list(),
                                                  epsilon = 0.05, delta = 1e-6) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Construct precision matrix using Gershgorin-based bounds
  # 
  # inputs:
  # 
  # - p                 (integer)           Matrix dimension
  # - y_c               (q_c-dim vector)    Covariate vector
  # - alpha_functions   (Named list)        each item is a function alpha_ij(y_c)
  # - structure_type    (string)            One of "banded", "exponential_decay", "sparse"
  # - structure_params  (list)              structure-specific parameters
  #   - k               (value)             'banded' parameter
  #   - gamma           (value)             'exponential_decay' parameter
  #   - s               (value)             'sparse' parameter
  #
  # - epsilon           (number)            Minimum absolute value threshold
  # - delta             (number)           Regularization parameter for spectral adjustment
  #
  #
  # outputs:
  #
  # - List with precision matrix and diagnostics  
  #   - precision_matrix   (p x p matrix)   precision matrix with 1's on the diagonal
  #
  #
  # ----------------------------------------------------------------------------
  
  
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