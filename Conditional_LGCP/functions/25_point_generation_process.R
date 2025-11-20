generate_cox_process_events <- function(X_functions, time_grid, T_max, 
                                        max_intensity = 5, seed = NULL) {
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: Generate point process events from log-intensity functions
  #
  # 
  # Input: 
  #
  # - X_functions    (p x m)
  # - time_grid      (m x 1)
  # - T_max          (scalar)
  # - max_intensity  (scalar)
  # - seed           (scalar)
  #
  # Output: 
  #
  # - list with event_times and event_counts
  #
  #   - [[1]]
  #     - event_times (unnamed list of p vectors)
  #     - event_counts (unnamed p-dim vector)
  #
  # 
  # ----------------------------------------------------------------------------
  
  if (!is.null(seed)) set.seed(seed)
  
  p <- nrow(X_functions)
  event_times <- list()
  event_counts <- numeric(p)
  
  for (i in 1:p) {
    # Convert to intensity: Lambda_i(t) = exp(X_i(t))
    log_intensities <- X_functions[i, ]
    intensities <- pmin(exp(log_intensities), max_intensity)  # Cap intensity
    
    # Generate events via thinning algorithm
    events_i <- c()
    
    # Upper bound for intensity (for thinning)
    max_int <- max(intensities)
    if (max_int <= 0) {
      event_times[[i]] <- numeric(0)
      event_counts[i] <- 0
      next
    }
    
    # Thinning algorithm implementation
    t <- 0  # Current time
    while (t < T_max) {
      # Generate next potential event time
      # Inter-arrival time ~ Exp(max_int)
      t <- t + rexp(1, max_int)
      
      if (t >= T_max) break
      
      # Compute actual intensity at time t via interpolation
      intensity_at_t <- approx(time_grid, intensities, t, rule = 2)$y
      
      # Accept with probability intensity_at_t / max_int
      if (runif(1) < intensity_at_t / max_int) {
        events_i <- c(events_i, t)
      }
    }
    
    if(length(events_i) == 0){
      event_times[[i]] <- numeric(0)
      event_counts[i] <- 0
    } else{
      event_times[[i]] <- sort(events_i)
      event_counts[i] <- length(events_i)      
    }

  }
  
  return(list(
    event_times = event_times,
    event_counts = event_counts
  ))
}
