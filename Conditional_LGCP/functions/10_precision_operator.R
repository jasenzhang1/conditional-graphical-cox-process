estimate_precision_operator_v3 <- function(C_cond, p, block = FALSE, MP = FALSE) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: 
  #
  # - v3: use psd_jitter to adaptively massage the matrices instead of a fixed gamma
  #
  # Input: 
  #
  # - C_cond        (pm x pm matrix)
  # - p             (scalar)
  # - block         (boolean) are we taking the inverse of each m x m matrix individually?
  # - MP            (boolean)  are we using moore penrose?
  #
  # 
  # Output: 
  #
  # - P_cond (pm x pm matrix)
  #
  # ----------------------------------------------------------------------------
  
  n_time <- dim(C_cond)[1] / p
  
  if(block){
    
    C_cond_block <- extract_block_structure_v2(C_cond, p, n_time)
    
    if(MP){
      P_cond_block <- lapply(C_cond_block, function(x){ginv(x)})
    } else{
      P_cond_block <- lapply(C_cond_block, function(x){solve_sym(x)})
    }
    
    P_cond <- assemble_block_matrix_v2(P_cond_block, p, block_size)
    
  } else{
    

    if(MP){
      P_cond <- ginv(C_cond)
    } else{
      # Add regularization: (pm x pm) + (pm x pm) = (pm x pm)
      C_cond_reg <- psd_jitter(C_cond)
      
      # Compute precision operator: (pm x pm)^{-1} = (pm x pm)
      P_cond <- solve_sym(C_cond_reg)
    }
    

  }
  
  return(P_cond)
}