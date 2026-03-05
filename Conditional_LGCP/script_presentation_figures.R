source('functions/98_presentation_figures.R')


# kernel <- 'RBF'
# kernel_params <- c(2, 3)
# p <- 10
# Tmax <- 10
# timepoints <- seq(0, Tmax, length.out = 200)
# 
# mat <- present_GP_funcs(p, timepoints, kernel, kernel_params)
# 
# visualize_log_intensity(t(mat), timepoints, '10 GP Draws')
# 
# 
# # simulate a drawing
# 
# events <- generate_cox_process_events(matrix(mat[,9], nrow = 1), timepoints, Tmax, max_intensity = Inf)
# 
# g_point_process <- visualize_intensity_with_points(mat[,9], timepoints, events$event_times[[1]])




# 3) draw events from a sine random function -----------------------------------

tmin <- -1
tmax <- 5
delta_t <- 0.001
mu <- 20
sigma <- 35
seed <- 3
event_y <- 2

l1 <- sine_random_function_with_points_generate(tmin, tmax, delta_t, mu, sigma, seed)

g1 <- sine_random_function_with_points_plot(l1$t_vec, l1$x_t, l1$events, event_y)

g1$graph_with_intensity

# 3b) 

tmin <- -1
tmax <- 5
delta_t <- 0.001
mu <- 10
sigma <- 6
seed <- 17
event_y <- 2

l2 <- sine_random_function_with_points_generate(tmin, tmax, delta_t, mu, sigma, seed)

g2 <- sine_random_function_with_points_plot(l2$t_vec, l2$x_t, l2$events, event_y)

g2$graph_with_intensity

# 3b) 

tmin <- -1
tmax <- 5
delta_t <- 0.001
mu <- 5
sigma <- 2
seed <- 63
event_y <- 2

l3 <- sine_random_function_with_points_generate(tmin, tmax, delta_t, mu, sigma, seed)

g3 <- sine_random_function_with_points_plot(l3$t_vec, l3$x_t, l3$events, event_y)

g3$graph_with_intensity

# 3b) 

tmin <- -1
tmax <- 5
delta_t <- 0.001
mu <- 10
sigma <- 8
seed <- 2
event_y <- 4
sin_prop <- c(0, 1, 0)

l2 <- sine_random_function_with_points_generate(tmin, tmax, delta_t, mu, sigma, seed, sin_prop)

g2 <- sine_random_function_with_points_plot(l2$t_vec, l2$x_t, l2$events, event_y)

g2$graph_with_intensity
