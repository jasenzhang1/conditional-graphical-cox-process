

# What do dependent intensities look like?
# And can we reconstruct the precision matrix from ground truth dependent intensities?

seed_num <- 1
set.seed(seed_num)


source('functions/20_simulation_function_wrapper.R')
source('functions/00_function_wrapper.R')

ncores <- parallel::detectCores() - 1


# 1) define simulation params -------------------------------------------

# time discretization
p <- 10
m <- 50
T_max <- 1
time_grid <- seq(0, T_max, length.out = m)


# m x m GP parameters
GP_mean <- 6
GP_mean_vec <- rep(GP_mean, length(time_grid))
GP_var <- 1
gamma <- 20

GP_cov <- generate_covariance_matrix(time_grid, kernel = 'rbf', gamma = gamma, variance = GP_var)

 
# 2) generate adjacency matrix ----------------------------------



corr_mat_truth <- 0.6 * matrix(1, p,p) + 0.4 * diag(p)
prec_mat_truth <- solve(corr_mat_truth)
p_corr_mat_truth <- -prec_mat_truth / sqrt(outer(diag(prec_mat_truth), diag(prec_mat_truth)))
diag(p_corr_mat_truth) <- 1

corr_mat_indep <- diag(p)
prec_mat_indep <- solve(corr_mat_indep)
p_corr_mat_indep <- -prec_mat_indep / sqrt(outer(diag(prec_mat_indep), diag(prec_mat_indep)))
diag(p_corr_mat_indep) <- 1

# 2.3) simulate dependent GP's


GP_mean_full <- rep(GP_mean, p*m)

GP_cov_full <- kronecker(corr_mat_truth, GP_cov)
GP_cov_full_indep <- kronecker(corr_mat_indep, GP_cov)



X_k <-        generate_log_intensity_functions_full_mat(GP_cov_full, time_grid, 
                                                        baseline_mean = GP_mean_vec, 
                                                        sample_mode = 'multi', seed = seed_num)

X_k_indep <-  generate_log_intensity_functions_full_mat(GP_cov_full_indep, time_grid, 
                                                        baseline_mean = GP_mean_vec, 
                                                        sample_mode = 'multi', seed = seed_num)

# 2.4) visualize the dependent and independent GP's

g <- visualize_log_intensity(X_k, time_grid)
g_indep <- visualize_log_intensity(X_k_indep, time_grid)

grid.arrange(g, g_indep, nrow = 1)

g_intensities <-  visualize_log_intensity(exp(X_k), time_grid)
g_intensities_indep <-  visualize_log_intensity(exp(X_k_indep), time_grid)

grid.arrange(g_intensities, g_intensities_indep, nrow = 1)

# 3) we skip the data generation and reconstruction of log-intensities. We immediately go to estimation.


