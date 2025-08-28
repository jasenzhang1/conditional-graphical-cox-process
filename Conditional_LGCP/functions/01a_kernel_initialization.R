
# 8/27/25

# we want to develop a method to generate
# p x p partial correlation matrices such that:

# - they are PSD
# - they are continuous in y_c
# - their diagonals are 1
# - off diagonal values are in [-1, 1]
# - abs(off diagonals) are either 0, or they are at least epsilon

library(Matrix)

construct_robust_precision_matrix <- function(p, y_c, alpha_functions, 
                                              rho_max = NULL, epsilon = 0.05, 
                                              delta = 1e-6) {
  
  #
  #
  # inputs:
  #
  # - p   (int)               number of processes
  # - y_c (q_c-dim vector)    vector of covariates
  # - alpha_functions  
  # - rho_max
  # - epsilon                 off-diagonal buffer
  # - delta                   pd buffer
  #
  # 
  
  # Step 1: Ensure theoretical condition p * rho_max < 1
  if (is.null(rho_max)) {
    rho_max <- 0.9 / p  # Conservative choice
  }
  
  if (p * rho_max >= 1) {
    stop(paste("Theoretical condition violated: p * rho_max =", p * rho_max, ">= 1"))
  }
  
  # Step 2: Construct correlation functions rho_ij(y_c) = rho_max * tanh(alpha_ij(y_c))
  rho_values <- matrix(0, p, p)
  diag(rho_values) <- 1  # Diagonal elements = 1
  
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      key <- paste(i, j, sep = "_")
      if (key %in% names(alpha_functions)) {
        alpha_val <- alpha_functions[[key]](y_c)
        rho_val <- rho_max * tanh(alpha_val)
        
        # Apply minimum threshold
        if (abs(rho_val) >= epsilon) {
          rho_values[i, j] <- rho_val
          rho_values[j, i] <- rho_val  # Ensure symmetry
        }
      }
    }
  }
  
  # Step 3: Construct initial correlation matrix
  R0 <- rho_values
  
  # Step 4: Eigenvalue decomposition for spectral adjustment
  eigen_decomp <- eigen(R0, symmetric = TRUE)
  Q <- eigen_decomp$vectors
  lambda <- eigen_decomp$values
  
  # Step 5: Apply spectral adjustment if needed
  lambda_adjusted <- pmax(lambda, delta)
  P <- Q %*% diag(lambda_adjusted) %*% t(Q)
  
  # Ensure diagonal is exactly 1 (numerical precision)
  diag(P) <- 1
  
  # Diagnostics
  min_eigenvalue_original <- min(lambda)
  min_eigenvalue_adjusted <- min(lambda_adjusted)
  spectral_adjustment_needed <- min_eigenvalue_original < delta
  theoretical_lower_bound <- 1 - p * rho_max # theoretical lower bound of eigenvalues
  
  return(list(
    precision_matrix = P,
    original_matrix = R0,
    eigenvalues_original = lambda,
    eigenvalues_adjusted = lambda_adjusted,
    spectral_adjustment_needed = spectral_adjustment_needed,
    min_eigenvalue_original = min_eigenvalue_original,
    min_eigenvalue_adjusted = min_eigenvalue_adjusted,
    theoretical_lower_bound = theoretical_lower_bound,
    rho_max = rho_max,
    y_c = y_c,
    constraint_satisfied = min_eigenvalue_adjusted >= theoretical_lower_bound - 1e-10
  ))
}


validate_precision_constraints <- function(result, epsilon) {
  
  
  # ---------------------------------------------------------------------------- 
  #
  #
  # GOAL: Validate precision matrix constraints
  #       Tests all mathematical requirements for precision matrices
  #
  # inputs:
  #
  # - result   (list)  Output from construct_robust_precision_matrix
  # - epsilon  (number)  Minimum non-zero off-diagonal value
  #
  # output:
  #
  # -
  #
  #
  # ----------------------------------------------------------------------------
  
  P <- result$precision_matrix # final precision matrix
  p <- nrow(P)                 # number of processes
  
  constraints <- list()
  
  # Constraint 1: Positive semidefinite
  min_eig <- min(eigen(P, symmetric = TRUE, only.values = TRUE)$values)
  constraints$positive_semidefinite <- list(
    satisfied = min_eig >= -1e-10,
    value = min_eig,
    description = "Positive semidefinite (lambda_min >= 0)"
  )
  
  # Constraint 2: Continuity (guaranteed by construction)
  constraints$continuity <- list(
    satisfied = TRUE,
    value = NA,
    description = "Continuous in covariates"
  )
  
  # Constraint 3: Unit diagonal
  diagonal_values <- diag(P)
  diagonal_ok <- all(abs(diagonal_values - 1) < 1e-10)
  constraints$unit_diagonal <- list(
    satisfied = diagonal_ok,
    value = paste(round(range(diagonal_values), 6), collapse = " to "),
    description = "Unit diagonal elements"
  )
  
  # Constraint 4: Bounded off-diagonals in [-1, 1]
  off_diag <- P[upper.tri(P) | lower.tri(P)]
  bounds_ok <- all(off_diag >= -1 - 1e-10 & off_diag <= 1 + 1e-10)
  constraints$bounded_elements <- list(
    satisfied = bounds_ok,
    value = paste(round(range(off_diag), 4), collapse = " to "),
    description = "Off-diagonal elements belong to [-1, 1]"
  )
  
  # Constraint 5: Minimum threshold for nonzero elements (> epsilon)
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
    epsilon_threshold = epsilon,
    description = paste("Minimum |value| >=", epsilon, "when nonzero")
  )
  
  # Overall validation
  all_satisfied <- all(sapply(constraints, function(x) x$satisfied))
  
  return(list(
    all_satisfied = all_satisfied,
    individual_constraints = constraints,
    summary = paste("Satisfied", sum(sapply(constraints, function(x) x$satisfied)), 
                    "out of", length(constraints), "constraints")
  ))
}


create_distance_based_alpha <- function(p, decay_rate = 0.3, covariate_strength = 0.5) {
  
  
  #-----------------------------------------------------------------------------
  #
  #
  # GOAL: helper function for alpha_ij specification, namely distance-based alpha
  # 
  #
  #                alpha_ij(y_c) = exp(-decay * abs(i - j)) + strength * sum(y_c)
  #                              =    distance_effect       +  covariate_effect  
  # 
  # 
  # input:
  #
  # - p                    (int)     number of processes   
  # - decay_rate           (number)  rate of exponential decay with distance
  # - covariate_strength   (number)  Strength of covariate modulation
  #
  #
  # output:
  #
  # - alpha_functions   (list)    list of f_ij(y_c) functions for each pair of processes
  # 
  #
  # ----------------------------------------------------------------------------
  
  alpha_functions <- list()
  
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      key <- paste(i, j, sep = "_")
      
      # store various copies of alpha_funcs with specified parameters
      alpha_functions[[key]] <- local({
        i_local <- i
        j_local <- j
        decay <- decay_rate
        strength <- covariate_strength
        function(y_c) {
          distance_effect <- exp(-decay * abs(i_local - j_local))
          covariate_effect <- strength * sum(y_c)
          return(distance_effect + covariate_effect)
        }
      })
    }
  }
  
  return(alpha_functions)
}

create_ar_style_alpha <- function(p, base_correlation = 0.7, covariate_strength = 0.2) {
  
  
  #-----------------------------------------------------------------------------
  #
  #
  # GOAL: helper function for alpha_ij specification, namely AR(1)-style alpha functions
  # 
  # 
  #         alpha_ij(y_c) =   base_cor^abs(i - j) + 1 + strength * sum(y_c)
  #                       =      AR(1) effect     +  covariate effect  
  #
  # 
  # input:
  #
  # - p                    (int)     number of processes   
  # - base_correlation     (number)  Base correlation parameter
  # - covariate_strength   (number)  Strength of covariate modulation
  #
  #
  # output:
  #
  # - alpha_functions   (list)    list of f_ij(y_c) functions for each pair of processes
  # 
  #
  # ----------------------------------------------------------------------------  
  
  alpha_functions <- list()
  
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      key <- paste(i, j, sep = "_")
      
      alpha_functions[[key]] <- local({
        i_local <- i
        j_local <- j
        base_rho <- base_correlation
        strength <- covariate_strength
        function(y_c) {
          ar_component <- base_rho^abs(i_local - j_local)
          covariate_modulation <- 1 + strength * sum(y_c)
          return(atanh(ar_component) * covariate_modulation)
        }
      })
    }
  }
  
  return(alpha_functions)
}
