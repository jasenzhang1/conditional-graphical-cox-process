# want to show that KDE recovers the underlying density

# 1) thinning algorithm to sample from a density
# 2) KDE to estimate the density from data
# 3) see if KDE is working/converging

seed_num <- 2
set.seed(seed_num)



source('functions/20_simulation_function_wrapper.R')
source('functions/00_function_wrapper.R')

ncores <- parallel::detectCores() - 1

m <- 50
T_max <- 1

time_grid <- seq(0, T_max, length.out = m)

GP_mean <- 6
GP_mean_vec <- rep(GP_mean, length(time_grid))

GP_var <- 1

GP_cov <- generate_covariance_matrix(time_grid, kernel = 'rbf', gamma = log(2 * m) * m, variance = GP_var)

# SANITY CHECK 1: 
# 
# check that the number of events generated is consistent with theory
# since our GP has constant mean and variance, the expected number of events is exp(mu + 0.5 * K(t,t))

num_events <- c()
n_simu <- 100

for(i in 1:n_simu){
  
  if(i %% 10 == 0){
    print(i)
  }
  
  X_k <- rmvnorm(1, 
                 mean = GP_mean_vec, 
                 sigma = GP_cov)
  


  events_k <- generate_cox_process_events(X_k, time_grid, T_max, max_intensity = Inf, seed = seed_num)
  num_events <- c(num_events, events_k$event_counts)
  
}



theoretical_event_num <- exp(GP_mean + 0.5 * GP_var)
empirical_event_num <- mean(num_events)

theoretical_event_num
empirical_event_num


# SANITY CHECK 2: 
# 
# do the events look like the drawn density?


# store ground truth log(intensity)

df_truth <- data.frame(time = time_grid, log_intensity = as.numeric(X_k))

# plot events on top of intensity

g_events <- ggplot() + geom_histogram(aes(x = events_k$event_times[[1]], y = after_stat(density)), bins = 50) + 
  geom_line(data = df_truth, aes(x = time, y = exp(log_intensity) / theoretical_event_num))
g_events


# SANITY CHECK 3:
#
# does the KDE estimate look like the ground truth intenstiy?

df_estimate <- data.frame(feature_id = 1, 
                          time = events_k$event_times[[1]], 
                          subject_num = 1) %>% as.data.table()

Tseq <- 0:200/200
patient_sel <- 1
feature_sel <- 1

# estimation procedure
rho_list <- estimate_intensities_stratum_parallel_v3(df_estimate, patient_sel, feature_sel, Tseq, ncores)

# store estimation
empirical_intensity <- rho_list[[1]]$rho_i

df_empirical <- data.frame(time = Tseq, intensity = empirical_intensity)


g_compare <- ggplot() + 
  geom_histogram(aes(x = events_k$event_times[[1]], y = after_stat(density)), bins = 50) + 
  geom_line(data = df_empirical, aes(x = time, y = intensity / theoretical_event_num, color = 'empirical')) + 
  geom_line(data = df_truth, aes(x = time, y = exp(log_intensity) / theoretical_event_num, color = 'truth')) 
  
g_compare
