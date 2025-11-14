library(MASS)
source('functions/28_Simulation_Visualization.R')
source('functions/25_point_generation_process.R')

GP_kernel <- function(timepoints, kernel, kernel_params) {
  
  #
  # - timepoints     (vector)    vector of timepoints
  # - kernel         (string)    'RBF', 'Matern'
  # - kernel_params  (vector)    vector of parameters
  #
  

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


kernel <- 'RBF'
kernel_params <- c(2, 3)
p <- 10
Tmax <- 10
timepoints <- seq(0, Tmax, length.out = 200)

mat <- present_GP_funcs(p, timepoints, kernel, kernel_params)

visualize_log_intensity(t(mat), timepoints, '10 GP Draws')


# simulate a drawing

events <- generate_cox_process_events(matrix(mat[,9], nrow = 1), timepoints, Tmax, max_intensity = Inf)

g_point_process <- visualize_intensity_with_points(mat[,9], timepoints, events$event_times[[1]])
