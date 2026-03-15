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

# 4) generate events from an intensity vector

lambda_1 <- read.csv('data/lambda_1.csv', header = T) %>%
  # Replace 'V1' with whatever the actual column name is (check with names(lambda_1))
  separate(1, into = c("time", "intensity"), sep = ",") %>%
  mutate(across(everything(), as.numeric))

lambda_2 <- read.csv('data/lambda_2.csv', header = T) %>%
  # Replace 'V1' with whatever the actual column name is (check with names(lambda_1))
  separate(1, into = c("time", "intensity"), sep = ",") %>%
  mutate(across(everything(), as.numeric))

lambda_3 <- read.csv('data/lambda_p.csv', header = T) %>%
  # Replace 'V1' with whatever the actual column name is (check with names(lambda_1))
  separate(1, into = c("time", "intensity"), sep = ",") %>%
  mutate(across(everything(), as.numeric))

m <- 4
shift <- 0.8

lambda_df <- data.frame(lambda_1$time, 
                        m*(lambda_1$intensity + shift), 
                        m*(lambda_2$intensity + shift), 
                        m*(lambda_3$intensity + shift))

colnames(lambda_df) <- c('time', 'l1', 'l2', 'l3')


ggplot() + geom_point(data = lambda_df, aes(x = time, y = l1), color = 'red') + 
  geom_point(data = lambda_df, aes(x = time, y = l2), color = 'blue') + 
  geom_point(data = lambda_df, aes(x = time, y = l3), color = 'green')


events_1 <- thinning(lambda_df$time, lambda_df$l1)
events_2 <- thinning(lambda_df$time, lambda_df$l2)
events_3 <- thinning(lambda_df$time, lambda_df$l3)
