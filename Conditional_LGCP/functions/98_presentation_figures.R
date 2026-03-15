library(MASS)
library(tidyverse)
source('functions/28_Simulation_Visualization.R')
source('functions/25_point_generation_process.R')



GP_kernel <- function(timepoints, kernel, kernel_params) {
  
  
  # ----------------------------------------------------------------------------
  # 
  # - timepoints     (vector)    vector of timepoints
  # - kernel         (string)    'RBF', 'Matern'
  # - kernel_params  (vector)    vector of parameters
  #
  # ----------------------------------------------------------------------------
  

  m <- length(timepoints)
  
  # Initialize empty matrix
  K <- matrix(0, nrow = m, ncol = m)
  
  # Compute pairwise distances
  dists <- as.matrix(dist(timepoints))
  
  if (kernel == "RBF") {
    # RBF / squared exponential kernel
    
    variance <- kernel_params[1]
    lengthscale <- kernel_params[2]
    
    K <- variance * exp(-0.5 * (dists / lengthscale)^2)
    
  } else if (kernel == "Matern") {
    # Matern kernel using fields::Matern
    
    variance <- kernel_params[1]
    lengthscale <- kernel_params[2]
    nu <- kernel_params[3]
    
    K <- variance * Matern(d = dists, range = lengthscale, smoothness = nu)
  }
  
  return(K)
}


present_GP_funcs <- function(p, timepoints, kernel, kernel_params){
  
  # ----------------------------------------------------------------------------
  # 
  # - p              (integer)         number of random function draws
  # - timepoints     (m-dim vector)    vector of timepoints
  # - kernel         (string)          'RBF', 'Matern'
  # - kernel_params  (vector)          vector of parameters
  #
  # output:
  #
  # - gp_funcs       (p x m matrix)    p realizations of random functions
  #
  # ----------------------------------------------------------------------------
  
  m <- length(timepoints)
  
  # Compute the kernel/covariance matrix
  K <- GP_kernel(timepoints, kernel, kernel_params)
  
  # Sample p realizations from multivariate normal with mean zero
  # MASS::mvrnorm returns p x m if we transpose
  gp_funcs <- t(MASS::mvrnorm(n = p, mu = rep(2, m), Sigma = K))
  
  # Output: p rows = p realizations, m columns = timepoints
  return(gp_funcs)
  
  
  
}

thinning <- function(time_vec, lambda_vec){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: thinning algorithm to draw events from an inhomogeneous point process
  #
  # inputs:
  #
  # - time_vec   (vector) times
  # - lambda_vec (vector) values of the intensity at certain times
  #
  # 
  # outputs:
  #
  # - events   (vector) vector of drawn events
  #
  # ----------------------------------------------------------------------------
  
  # 1. Clean data: Ensure intensity is non-negative
  # (Using pmax because intensities cannot be negative for thinning)
  lambda_t <- pmax(0, lambda_vec) 
  
  # 2. Get bounds from the data provided
  t_min <- min(time_vec)
  t_max <- max(time_vec)
  max_lambda <- max(lambda_t)
  
  # If max_lambda is 0, no events can occur
  if (max_lambda <= 0) return(numeric(0))
  
  # 3. Generate potential events from a Homogeneous Poisson Process
  # Intensity * duration
  n_potential <- rpois(1, lambda = max_lambda * (t_max - t_min))
  
  if (n_potential == 0) return(numeric(0))
  
  potential_events <- sort(runif(n_potential, t_min, t_max))
  
  # 4. Thinning: accept points with probability lambda(t) / max_lambda
  # Use linear interpolation to find lambda at the exact potential event times
  event_intensities <- approx(x = time_vec, y = lambda_t, xout = potential_events)$y
  
  acceptance_prob <- event_intensities / max_lambda
  keep <- runif(n_potential) <= acceptance_prob
  
  events <- potential_events[keep]
  
  return(events)
  
}

sine_random_function_with_points_generate <- function(tmin, tmax, delta_t, mu, sigma, seed, sin_prop = NULL){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: generate a random function + point process points from the following:
  #
  #        x(t) = mu(t) + a1 * sin(100t) + a2 * sin(189t) + a3 * sin(71t) 
  # 
  #        where a1, a2, a3 are drawn uniformly from [0, 1] and sum to 1. So use stick breaking to find them. 
  # 
  # inputs:
  #
  # - tmin     (scalar) 
  # - tmax     (scalar)
  # - delta_t  (scalar)  discretization of time
  # - mu       (scalar)  constant value for constant mu(t) function
  # - sigma    (scalar)  how much to amplify each sine function
  # - seed     (scalar)  seed for replicability
  # 
  # 
  # outputs:
  # 
  # - list of the following:
  #   - t_vec      (vector)  time discretization vector
  #   - x_t        (vector)  realized random function starting from tmin to tmax at step sizes of delta_t
  #   - events     (vector)  realized timestamps of the events
  #   - sin_prop   (vector)  [a1, a2, a3] vector of allocated proportions to each sine function
  # 
  # ----------------------------------------------------------------------------
  
  set.seed(seed)
  
  # 1. Generate random phase shifts for each sine wave
  # Range [0, 2*pi]
  phases <- runif(3, 0, 2 * pi)
  

  
  # 2. Stick-breaking to find a1, a2, a3 ONLY if not provided
  if (is.null(sin_prop)) {
    breaks <- sort(runif(2))
    a1 <- breaks[1]
    a2 <- breaks[2] - breaks[1]
    a3 <- 1 - breaks[2]
    sin_prop <- c(a1, a2, a3)
  }
  
  # 3. Define time vector
  t_vec <- seq(tmin, tmax, by = delta_t)
  
  # 4. Generate the random function x(t)
  # x(t) = mu + sigma * (a1*sin(100t) + a2*sin(189t) + a3*sin(71t))
  x_t <- mu + sigma * (
    sin_prop[1] * sin(0.40 * t_vec + phases[1]) + 
    sin_prop[2] * sin(1.19 * t_vec + phases[2]) + 
    sin_prop[3] * sin(2.12 * t_vec + phases[3])
  )
  
  # 5. Generate Point Process (Inhomogeneous Poisson Process) via Thinning
  # Ensure intensity lambda(t) is non-negative
  lambda_t <- pmax(0, x_t)
  max_lambda <- max(lambda_t)
  
  # Generate potential events from a Homogeneous Poisson Process
  n_potential <- rpois(1, lambda = max_lambda * (tmax - tmin))
  potential_events <- sort(runif(n_potential, tmin, tmax))
  
  # Thinning: accept points with probability lambda(t) / max_lambda
  if (n_potential > 0) {
    # Interpolate to find intensity at specific random potential_event times
    event_intensities <- approx(t_vec, lambda_t, xout = potential_events)$y
    acceptance_prob <- event_intensities / max_lambda
    keep <- runif(n_potential) <= acceptance_prob
    events <- potential_events[keep]
  } else {
    events <- numeric(0)
  }
  
  # 6. Return results
  return(list(
    t_vec = t_vec,
    x_t = x_t,
    events = events,
    sin_prop = sin_prop
  ))
  
}



sine_random_function_with_points_plot <- function(t_vec, x_t, events, event_y, event_alpha = 0.5) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: continuation of sine_random_function_with_points_generate, 
  #       plot the realized random function with its points
  #
  # inputs:
  # 
  # - t_vec       (vector)  time discretization vector
  # - x_t         (vector)  intensity function vector
  # - events      (vector)  realized timestamps of events
  # - event_y     (scalar)  y-value for which all the events are lying on
  # - event_alpha (scalar)  opacity of the event points (0 to 1)
  #
  # outputs:  
  # 
  # - a list of:
  #   - graph_events_only: point processes timestamps alone
  #   - graph_with_intensity: point processes timestamps with intensity function 
  #
  # ----------------------------------------------------------------------------
  
  # Data frames
  df_func <- data.frame(time = t_vec, intensity = x_t)
  df_events <- data.frame(time = events, y = event_y)
  
  # Define the minimal theme inside the function to avoid scope errors
  minimal_y_only_theme <- theme_minimal() +
    theme(
      # Remove x-axis entirely
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      # Add y-axis line and keep labels
      axis.line.y = element_line(color = "black"),
      # Remove all grid lines for a clean background
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # 1. Graph with point processes timestamps alone
  g1 <- ggplot(df_events, aes(x = time, y = y)) +
    # Events as black dots with custom opacity
    geom_point(color = "black", alpha = event_alpha) +
    labs(title = "Point Process Events", y = "") +
    # Vertical windowing around the event_y line
    ylim(event_y - 0.5, event_y + 0.5) +
    minimal_y_only_theme +
    theme(
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )
  
  # 2. Graph with point processes timestamps with intensity function
  g2 <- ggplot(df_func, aes(x = time, y = intensity)) +
    # The intensity function line
    geom_line(color = "royalblue", linewidth = 0.8) +
    # Events as black dots at the specified event_y level
    geom_point(data = df_events, 
               aes(x = time, y = y), 
               color = "black", 
               alpha = event_alpha, 
               inherit.aes = FALSE) +
    labs(title = "Intensity Function x(t) with Events", y = "x(t)") +
    minimal_y_only_theme
  
  # Return as a list of ggplot objects
  return(list(
    graph_events_only = g1,
    graph_with_intensity = g2
  ))
}
