estimate_precision_operator <- function(C_conditional, gamma2, p) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: 
  #
  #
  # Input: 
  #
  # - C_conditional (list of length p^2, each element m x m)
  # - gamma2 (scalar)
  # - p (scalar)
  #
  # 
  # Output: 
  #
  # - P_conditional (list of length p^2, each element m x m)
  #
  # ------------------------------------------------------------------------
  
  n_time <- nrow(C_conditional[[1]])  # m
  
  # Assemble block correlation matrix: pm x pm
  C_block <- assemble_block_matrix(C_conditional, p, n_time)
  
  # Add regularization: (pm x pm) + (pm x pm) = (pm x pm)
  C_block_reg <- C_block + gamma2 * diag(p * n_time)
  
  # Compute precision operator: (pm x pm)^{-1} = (pm x pm)
  P_block <- solve(C_block_reg)
  
  # Extract block structure: pm x pm -> list of p^2 blocks (m x m each)
  P_conditional <- extract_block_structure(P_block, p, n_time)
  
  return(P_conditional)
}

estimate_precision_operator_v2 <- function(C_conditional, gamma2, p) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: 
  #
  #
  # Input: 
  #
  # - C_conditional (list of length p^2, each element m x m)
  # - gamma2 (scalar)
  # - p (scalar)
  #
  # 
  # Output: 
  #
  # - P_conditional (list of length p^2, each element m x m)
  #
  # ------------------------------------------------------------------------
  
  n_time <- nrow(C_conditional[[1]])  # m
  
  # Assemble block correlation matrix: pm x pm
  C_block <- assemble_block_matrix_v2(C_conditional, p, n_time)
  
  # Add regularization: (pm x pm) + (pm x pm) = (pm x pm)
  C_block_reg <- C_block + gamma2 * diag(p * n_time)
  
  # Compute precision operator: (pm x pm)^{-1} = (pm x pm)
  P_block <- solve(C_block_reg)
  
  # Extract block structure: pm x pm -> list of p^2 blocks (m x m each)
  P_conditional <- extract_block_structure_v2(P_block, p, n_time)
  
  return(P_conditional)
}


estimate_precision_operator_v3 <- function(C_conditional, p) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: 
  #
  # - v3: use psd_jitter to adaptively massage the matrices instead of a fixed gamma
  #
  # Input: 
  #
  # - C_conditional (list of length p^2, each element m x m)
  # - p (scalar)
  #
  # 
  # Output: 
  #
  # - P_conditional (list of length p^2, each element m x m)
  #
  # ------------------------------------------------------------------------
  
  n_time <- nrow(C_conditional[[1]])  # m
  
  # Assemble block correlation matrix: pm x pm
  C_block <- assemble_block_matrix_v2(C_conditional, p, n_time)
  
  # Add regularization: (pm x pm) + (pm x pm) = (pm x pm)
  C_block_reg <- psd_jitter(C_block)
  
  # Compute precision operator: (pm x pm)^{-1} = (pm x pm)
  P_block <- solve(C_block_reg)
  
  # Extract block structure: pm x pm -> list of p^2 blocks (m x m each)
  P_conditional <- extract_block_structure_v2(P_block, p, n_time)
  
  return(P_conditional)
}
