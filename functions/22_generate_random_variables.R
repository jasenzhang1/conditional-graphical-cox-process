# all functions to generate y_c

# 1) generate continuous and discrete RV's
# 2) generate only continuous RV's because we assume 1 strata
# 3) generate continuous RV's that reflect week; assume 1 strata
# 4) simulation_seed: deterministic seeds for data generation



simulation_seed <- function(adj_type, rep_i, stage, base_seed = NULL) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: deterministic RNG seed for one piece of the simulation, so every
  #       replication (and every group within it) can be regenerated on its own
  #
  # inputs:
  #
  # - adj_type    (string)    simulation setting, e.g. 'hub_block_v2'
  # - rep_i       (integer)   replication index
  # - stage       (integer)   0 = conditioning variables (part 0),
  #                           g = event generation for group g (parts 1 and 2)
  # - base_seed   (integer)   defaults to the CGCP_SEED environment variable,
  #                           or 2025 if it is unset
  #
  # output:
  #
  # - integer seed in [1, 2^31 - 2]
  #
  # ----------------------------------------------------------------------------
  
  if (is.null(base_seed)) base_seed <- as.numeric(Sys.getenv('CGCP_SEED', '2025'))
  
  adj_hash <- sum(utf8ToInt(adj_type) * seq_along(utf8ToInt(adj_type)))
  seed <- (base_seed * 7919 + adj_hash * 104729 + rep_i * 1009 + stage) %% 2147483646 + 1
  
  as.integer(seed)
}




generate_y_c_adj_type <- function(n, adj_type, params, seed = NULL){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: depending on the adj_type, generate Y_continuous
  #
  # - for now, assume q_c = 1, so Y_continuous is just a scalar for each subject
  #
  # inputs:
  #
  # - n               (integer)          sample size
  # - adj_type        (string)           precision matrix adjacency type (e.g. single_c2)
  # - params          (vector)           associated vector of parameters of the adjacency type
  # - seed            (integer)          reproducibility seed 
  #
  # 
  # output:
  #
  # - Y_continuous (n x q_c matrix)
  #
  #
  # ----------------------------------------------------------------------------
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # extract last character
  
  y_c_type <- substring(adj_type, nchar(adj_type))
  
  if(y_c_type == '0'){
    
    # params[1] = value
    # all covariates are the same
    
    Y_continuous <- matrix(rep(params[1], n), nrow = n, ncol = 1)
    
    return(Y_continuous)
    
  }   
  
  if(y_c_type == '1'){
    
    # params[1] = min
    # params[2] = max
    # y_c ~ uniform[min, max]
    
    # the premise is that regardless of the value of Y, the underlying graph is the same
    # randomly sample between min and max
    
    Y_continuous <- matrix(runif(n, params[1], params[2]), nrow = n, ncol = 1)
    
    return(Y_continuous)
    
  }
  
  if(y_c_type == '2'){
    
    # params[1] = min
    # params[2] = max
    # equal spacing between all y_c_k's
    
    Y_continuous <- matrix(seq(params[1], params[2], length.out = n), nrow = n, ncol = 1)
    
    return(Y_continuous)
    
  }  
  
  if(y_c_type == '3'){
    
    # params[1] = min
    # params[2] = max
    # params[3] = num_groups
    # mimic weeks. have y_c take on repeated values spread throughtout [min, max]
    m <- params[3]
    values <- seq(params[1], params[2], length.out = m)
    
    counts <- rep(floor(n / m), m)
    counts[1:(n %% m)] <- counts[1:(n %% m)] + 1
    
    # repeat each value accordingly
    Y_continuous <- matrix(rep(values, counts), nrow = n, ncol = 1)
    
  }
  

}  
    
