# goal:
#        see if bivariate intensities are accurately reconstructed



# algorithm
#
# 1) generate a ground truth GP(mean, var) process
# 2) draw a log-intensity, exponentiate to intensity
# 2a) get the ground truth bivariate intensity 
# 3) draw events from the intensity
# 4) estiate the intensity
# 5) estimate the bivariate intensity 


# ------------------------------------------------------------------------------

seed_num <- 1
set.seed(seed_num)


source('functions/20_simulation_function_wrapper.R')
source('functions/00_function_wrapper.R')

ncores <- parallel::detectCores() - 1


# simulation params -------------------------------------------

m <- 50
T_max <- 1

time_grid <- seq(0, T_max, length.out = m)

GP_mean <- 6
GP_mean_vec <- rep(GP_mean, length(time_grid))

GP_var <- 1

GP_cov <- generate_covariance_matrix(time_grid, kernel = 'rbf', gamma = log(2 * m) * m, variance = GP_var)


# get our ground truths ---------------------------------

# log-intensity
# intensity
# bivariate intensity 

X_k <- rmvnorm(1, 
               mean = GP_mean_vec, 
               sigma = GP_cov)       # drawn log-intensity



rho_i_truth <- exp(X_k)

rho_ii_truth <- 

# draw events from thinning algorithm ---------------------------------

events_k <- generate_cox_process_events(X_k, time_grid, T_max, max_intensity = Inf, seed = seed_num)




