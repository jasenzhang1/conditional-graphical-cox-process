

inner_1_truth_reconstruct <- t(efuncs_1_truth) %*% efuncs_1_truth
inner_1_est_reconstruct <- t(efuncs_1_est) %*% efuncs_1_est

outer_1_truth_reconstruct <- efuncs_1_truth %*% t(efuncs_1_truth)
outer_1_est_reconstruct <- efuncs_1_est %*% t(efuncs_1_est)

inner_1_truth_reconstruct_norm <- t(efuncs_1_truth / sqrt(50)) %*% efuncs_1_truth / sqrt(50)
inner_1_est_reconstruct_norm <- t(efuncs_1_est / sqrt(20)) %*% efuncs_1_est / sqrt(20)

outer_1_truth_reconstruct_norm <- (efuncs_1_truth / sqrt(50)) %*% t(efuncs_1_truth / sqrt(50))
outer_1_est_reconstruct_norm <- (efuncs_1_est / sqrt(20)) %*% t(efuncs_1_est / sqrt(20))

# using data

visualize_matrix_heatmap(outer_1_truth_reconstruct, zmid = 0)

ef1 <- trig_basis(5)
ef2 <- trig_basis(10)
m1 <- 100
m2 <- 200
time_grid_1 <- make_time_grid(m1)
time_grid_2 <- make_time_grid(m2)

mat1 <- trig_basis_realization(ef1, time_grid_1)
mat2 <- trig_basis_realization(ef2, time_grid_2)

outer_1 <- mat1 %*% t(mat1)
outer_2 <- mat2 %*% t(mat2)


visualize_matrix_heatmap(outer_1, zmid = 0)
visualize_matrix_heatmap(outer_2, zmid = 0)
hilbert_schmidt_norm_diff(outer_1, outer_2, time_grid_1, time_grid_2)


outer_1_norm <- (mat1 / sqrt(m1)) %*% t(mat1 / sqrt(m1))
outer_2_norm <- (mat2 / sqrt(m2)) %*% t(mat2 / sqrt(m2))

grid.arrange(visualize_matrix_heatmap(outer_1_norm, zmid = 0),
             visualize_matrix_heatmap(outer_2_norm, zmid = 0), nrow = 1)
hilbert_schmidt_norm_diff(outer_1_norm, outer_2_norm, time_grid_1, time_grid_2)


