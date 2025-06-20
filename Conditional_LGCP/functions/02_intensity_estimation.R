estimate_intensities_stratum <- function(data_all, patient_sel, feature_sel, 
                                         t_seq) {
  
  #  
  # estimate the first and second order intensities
  # - rho_i  (p x m matrix)
  # - rho_{ij} (m x m matrix) for all (p x p) pairs
  #  
  #
  # Input: 
  #
  #
  # - data_all   (data.table):  'feature_id', 'time', 'subject_num' 
  # - patient_sel (vector)    ID's of non-discarded subjects
  # - feature_sel (vector)    ID's of non-discarded features
  # - t_seq   (vector of length m), 
  #
  #
  #
  # Output: 
  #
  # - rho_hat (p x m matrix), 
  # - rho_hat_pairs (p x p x m x m array)
  #
  # ---------------------------------------------------------------------------
  
  NN = length(patient_sel)
  p <- length(feature_sel)
  n_time <- length(t_seq)
  
  # First-order intensities: p x m matrix
  rho_hat <- matrix(0, nrow=p, ncol=n_time)
  
  # Second-order intensities: p x p x m x m array  
  rho_hat_pairs <- array(0, dim=c(p, p, n_time, n_time))
  
  for (i in 1:p) {
    
    data_i = data_all[feature_id==i,]
    data_i_count = data_i[,.(count=.N),by="subject_num"] # xi_i's
    
    if(nrow(data_i) == 0){
      rho_i <- rep(0, length(t_seq))
    } else{
      
    
      h_i <- get_gamma(data_i$time)
        
      Gamma_i <- data_i[,estimate_density(time, t_seq, h_i),by=c("subject_num")]   # vertically stacked \Gamma_i^k's
      
      
      patient_now = unique(Gamma_i$subject_num)                            # patients left after filtering for process i
      data_i_count = data_i_count[match(patient_now, subject_num),]        # filter out subjects that are not in "patient_now". 
      
      density_mat = matrix(Gamma_i$V1,nrow=length(t_seq))                    # rearrange to be a 19 x 380 matrix
      
      xi_i_mat <- matrix(data_i_count$count, nrow = 1)
      
      density_xi_mat <- sweep(density_mat, 2, xi_i_mat, `*`)         # multiply by xi_i (19 x 89 times 1 x 89)
      
      rho_i = apply(density_xi_mat, 1, sum) / NN                             # rho_i(t) for mark i 
      
      
    }
    
    rho_hat[i,] <- rho_i
    
    # Estimate second-order intensities (double loop over time grid)
    for (j in 1:p) {
      
      print(paste0(i, ',', j))
      
      data_j = data_all[feature_id==j,]
      data_j_count = data_j[,.(count=.N),by="subject_num"] # xi_j's  
      
      if(nrow(data_j) == 0 || nrow(data_i) == 0){ # moot case
        rho_ij_mat <- matrix(0, nrow = length(t_seq), ncol = length(t_seq))
      } else{
        times_i <- data_i[, .(event_times_i = list(time)), by = subject_num]
        times_j <- data_j[, .(event_times_j = list(time)), by = subject_num]
        
        times_ij <- merge(times_i, times_j, by = "subject_num")
        
        Gamma_ij <- times_ij[, estimate_bivariate_density(event_times_i[[1]], 
                                                          event_times_j[[1]], 
                                                          eval_grid_s, 
                                                          eval_grid_t
                                                          ),  by = 'subject_num']
        
        bivariate_density_mat <- matrix(Gamma_ij$V1,nrow=length(t_seq)^2)    # rarrange to be a 361 x 380 matrix
        
        data_ij_count <- data_i_count %>% full_join(data_j_count, by = 'subject_num') %>% na.omit()
        
        colnames(data_ij_count) <- c('subjet_num', 'count_i', 'count_j')
        
        xi_ij_mat <- matrix(data_ij_count$count_i * data_ij_count$count_j, nrow = 1)
        
        bivariate_density_xi_mat <- sweep(bivariate_density_mat, 2, xi_ij_mat, `*`) 
        
        rho_ij = apply(bivariate_density_xi_mat, 1, sum) / NN               # rho_ij(s, t) = 19 x 19  
        rho_ij_mat <- matrix(rho_ij, nrow = length(t_seq), ncol = length(t_seq))
      }
      
      rho_hat_pairs[i, j, ,] <- rho_ij_mat
      
      
    } # end of j
  } # end of i
  
  return(list(rho_hat = rho_hat, rho_hat_pairs = rho_hat_pairs))
}
