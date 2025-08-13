# goal:
#        see if bivariate intensities are accurately reconstructed



# algorithm
#
# 1) define a ground truth GP(mean, var) process
# 2) draw n log-intensities, exponentiate to intensity
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

n <- 50000
m <- 500
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



rho_i_truth <- exp(GP_mean + 0.5 * GP_var)


rho_i_truth_vec <- exp(GP_mean_vec + 0.5 * diag(GP_cov))

rho_ii_truth <- (rho_i_truth_vec %o% rho_i_truth_vec) * exp(GP_cov)

# skip drawing events, assume we get full reconstruction of intensities --------


process_i_means <- apply(exp(X_k), 1, mean)

rho_i_est_cum <- cumsum(process_i_means) / seq_along(process_i_means)


# since we have constant mean, we can take the average across all timepoints and subjects
ggplot() + geom_point(aes(x = 1:n, y = rho_i_est_cum)) + 
  geom_hline(yintercept = rho_i_truth) +
  ylim(rho_i_truth - 50, rho_i_truth + 50)



rho_i_est <- apply(exp(X_k), 2, mean)

ggplot() + geom_line(aes(x = time_grid, y = rho_i_truth_vec, color = 'truth')) + 
  geom_line(aes(x = time_grid, y = rho_i_est, color = 'estimate')) + 
  ylim(rho_i_truth - 20, rho_i_truth + 20)

# bivariate intensity

rho_ii_est <- crossprod(exp(X_k)) / n

summary(as.numeric(1 - rho_ii_est/rho_ii_truth))  # 3% to 5% error


# ratio

g_ij_st_truth1 <- log(rho_ii_truth / (rho_i_truth_vec %o% rho_i_truth_vec))
g_ij_st_truth <- GP_cov

g_ij_st_est <- log(rho_ii_est / (rho_i_est %o% rho_i_est))

summary(as.numeric(1 - g_ij_st_est/g_ij_st_truth))  # very far apart


df1 <- reshape2::melt(g_ij_st_est)
df2 <- reshape2::melt(g_ij_st_truth)

# Heatmap plots
p1 <- ggplot(df1, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed()

p2 <- ggplot(df2, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed()

# Arrange side-by-side
grid.arrange(p1, p2, ncol = 2) # est left, truth right

## DONE 

## WE HAVE MEAN INTENSITY AND MEAN BIVARIATE INTENSITY ESTIMATES


  
  