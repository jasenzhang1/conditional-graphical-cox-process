
# 1) generate a log-intensity from a GP
# 2) do eigendecomposition on GP_cov
# 3) estimate the kl coefficients to get a linear combination of the eigenfunctions to match the log-intensity


# ------------------------------------------------------------------------------

seed_num <- 1
set.seed(seed_num)


source('functions/20_simulation_function_wrapper.R')
source('functions/00_function_wrapper.R')

ncores <- parallel::detectCores() - 1


# simulation params -------------------------------------------

n <- 10
m <- 50
p <- 1
T_max <- 1

time_grid <- seq(0, T_max, length.out = m)
  
GP_mean <- 6
GP_mean_vec <- rep(GP_mean, length(time_grid))
GP_var <- 1
GP_gamma <- 20

GP_cov <- generate_covariance_matrix(time_grid, kernel = 'rbf', gamma = GP_gamma, variance = GP_var)


# get our ground truths ---------------------------------

# log-intensity
# intensity
# bivariate intensity 

X_k <- rmvnorm(n, 
               mean = GP_mean_vec, 
               sigma = GP_cov)       # drawn log-intensity

X_k_centered <- X_k - GP_mean_vec

visualize_log_intensity(X_k_centered, time_grid)
  
# eigendecomposition

eigen_decomp_truth <- compute_eigendecomposition_ii(array(GP_cov, dim = c(1, m, m)))

eigenfunctions <- t(eigen_decomp_truth$eigenfunctions[[1]])

visualize_log_intensity(eigenfunctions, time_grid) # the eigenfunctions look like sine waves, good

orthogonality_check <- eigenfunctions %*% t(eigenfunctions) # the eigenfunctions are orthonormal yay


# kl expansion

X_k_centered_prep <- array(t(X_k_centered), dim = c(p, m, n))

kl_coeffs_v2 <- estimate_kl_coefficients_parallel_v2(X_k_centered_prep, 
                                                     eigen_decomp_truth$eigenfunctions, 
                                                     time_grid, 
                                                     ncores)



X_k_reconstruct <- validate_kl_coeffs(eigen_decomp_truth$eigenfunctions, kl_coeffs_v2) # obtain X_k from the kl coeffs and eigenfunctions

df1 <- reshape2::melt(X_k_centered)
df2 <- reshape2::melt(t(X_k_reconstruct[1,,]))
df1$Var2 <- time_grid[df1$Var2]
df2$Var2 <- time_grid[df2$Var2]

df1$class <- 'truth'
df2$class <- 'reconstructed'
df_graph <- rbind(df1, df2)
df_graph$class <- factor(df_graph$class)

ggplot() + geom_line(data = df_graph, aes(x = Var2, y = value, group = class, color = class)) +
  facet_wrap(~ Var1 )



