# 1) generate simulated data
# 2) estimate ground truth
# 3) calculate accuracy metrics

source('functions/20_simulation_function_wrapper.R')

n = 100                     # Sample size (n)
p = 10                      # Number of processes (p)  
T_max = 10                  # Time horizon (T)
q_c = 2                     # Continuous conditioning dimension (q_c)
K = 4                       # Discrete combinations (K)
sparsity = 0.2              # Graph sparsity (s)
theta = 1.0                 # Signal strength (theta)
graph_type = "random"       # Graph topology
dependence_type = "linear"  # Conditional dependence type
time_grid_size = 50         # Time discretization (m)
seed = 1

dataset <- simulate_conditional_cox_data(n, p, T_max, q_c, K,
                                         sparsity,
                                         theta,
                                         graph_type,
                                         dependence_type,
                                         time_grid_size,
                                         seed)

# fitting

query_yd <- unique(dataset$Y_discrete)
query_yc <- rep(0, q_c)


