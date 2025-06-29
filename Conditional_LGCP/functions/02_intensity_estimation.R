library(future.apply)
library(dplyr)
library(data.table)
library(pbmcapply)
estimate_intensities_stratum <- function(data_all, patient_sel, feature_sel, 
                                         t_seq) {
  
  # ----------------------------------------------------------------------------  
  #
  #
  # GOAL: estimate the first and second order intensities
  #
  # - rho_i    (p x m matrix)
  # - rho_{ij} (m x m matrix) for all (p x p) pairs
  #  
  #
  # Input: 
  #
  #
  # - data_all       (data.table):  'feature_id', 'time', 'subject_num' 
  # - patient_sel    (vector)        subject ID's of non-discarded subjects
  # - feature_sel    (vector)        neuron ID's of non-discarded features
  # - t_seq          (vector of length m) 
  #
  #
  #
  # Output: 
  #
  # - rho_hat       (p x m matrix) 
  # - rho_hat_pairs (p x p x m x m array)
  #
  # ----------------------------------------------------------------------------
  

  
  NN = length(patient_sel) # 89
  p <- length(feature_sel) # 249
  n_time <- length(t_seq)  # 19
  
  # First-order intensities: p x m matrix
  rho_hat <- matrix(0, nrow=p, ncol=n_time)
  
  # Second-order intensities: p x p x m x m array  
  rho_hat_pairs <- array(0, dim=c(p, p, n_time, n_time))
  
  for (i in 1:p) {
    
    data_i = data_all[feature_id==i,]
    data_i_count = data_i[,.(count=.N),by="subject_num"] # xi_i's
    
    # 1a) Estimate first-order intensity (rho_i = 19 dim vector)
    
    if(nrow(data_i) == 0){ # trivial case where process i has no events, it's intensity (rho) is zero
      rho_i <- rep(0, length(t_seq))
    } else{
      
    
      h_i <- get_gamma(data_i$time)
        
      Gamma_i <- data_i[,estimate_density(time, t_seq),by=c("subject_num")]   # vertically stacked \Gamma_i^k's (89 * 19 length matrix)
      
      
      # patient_now = unique(Gamma_i$subject_num)                            # patients left after filtering for process i
      # data_i_count = data_i_count[match(patient_now, subject_num),]        # filter out subjects that are not in "patient_now". 
      
      # density_mat = matrix(Gamma_i$gamma_hat,nrow=length(t_seq))                    # rearrange to be a 19 x 89 matrix
      # 
      # xi_i_mat <- matrix(data_i_count$count, nrow = 1) # (1 x 89 mat) xi's in neuron number order
      # 
      # density_xi_mat <- sweep(density_mat, 2, xi_i_mat, `*`)         # multiply by xi_i (19 x 89 times 1 x 89)
      # 
      # rho_i = apply(density_xi_mat, 1, sum) / NN                             # rho_i(t) for mark i (length 19 vector)
      
      rho_mat <- matrix(Gamma_i$rho_hat, nrow = n_time)
      rho_i <- apply(rho_mat, 1, sum) / NN
      
      
    }
    
    rho_hat[i,] <- rho_i
    
    # 1b) Estimate second-order intensities (double loop over time grid)
    for (j in i:p) {
      
      print(paste0(i, ',', j))
      
      data_j = data_all[feature_id==j,]
      data_j_count = data_j[,.(count=.N),by="subject_num"] # xi_j's  
      
      if(nrow(data_j) == 0 || nrow(data_i) == 0){ # moot case where either process has no events
        rho_ij_mat <- matrix(0, nrow = length(t_seq), ncol = length(t_seq))
      } else{
        times_i <- data_i[, .(event_times_i = list(time)), by = subject_num]    # concatenate all the event times of rep k into a vector, then make dataframe
        times_j <- data_j[, .(event_times_j = list(time)), by = subject_num]
        
        times_ij <- merge(times_i, times_j, by = "subject_num", all = TRUE)
        
        Gamma_ij <- times_ij[, estimate_bivariate_density(event_times_i[[1]], 
                                                          event_times_j[[1]], 
                                                          t_seq, 
                                                          t_seq,
                                                          'i'
                                                          ),  by = 'subject_num']
        
        # Gamma_ij_troubleshoot <- list()
        # for(k in 1:nrow(times_ij)){
        #   Gamma_ij_k <- estimate_bivariate_density(times_ij$event_times_i[[k]],
        #                                            times_ij$event_times_j[[k]],
        #                                            t_seq,
        #                                            t_seq, 
        #                                            'i')
        # 
        #   Gamma_ij_troubleshoot[[k]] <- Gamma_ij_k
        # }
        
        
        # bivariate_density_mat <- matrix(Gamma_ij$V1,nrow=length(t_seq)^2)    # rarrange to be a [361 x 89] mat
        # 
        # # calculate rho_ij from Gamma_ij
        # 
        # data_ij_count <- data_i_count %>% full_join(data_j_count, by = 'subject_num')
        # 
        # colnames(data_ij_count) <- c('subjet_num', 'count_i', 'count_j')
        # 
        # if(i == j){
        #   # (xi) (xi - 1) correction
        #   xi_ij_mat <- matrix(data_ij_count$count_i * (data_ij_count$count_j - 1), nrow = 1)  # [1 x 89] mat
        #   xi_ij_mat[is.na(xi_ij_mat)] <- 0
        # } else{
        #   # xi_i * xi_j (full_join created some NA's so we need to set them equal to zero)
        #   xi_ij_mat <- matrix(data_ij_count$count_i * data_ij_count$count_j, nrow = 1)  # [1 x 89] mat
        #   xi_ij_mat[is.na(xi_ij_mat)] <- 0
        # }
        # 
        # 
        # bivariate_density_xi_mat <- sweep(bivariate_density_mat, 2, xi_ij_mat, `*`)   # [361 x 89 mat]
        
        bivariate_intensity <- matrix(Gamma_ij$V1,nrow=length(t_seq)^2)
        
        rho_ij = apply(bivariate_intensity, 1, sum) / NN               # rho_ij(s, t) = 19 x 19  
        rho_ij_mat <- matrix(rho_ij, nrow = length(t_seq), ncol = length(t_seq))
      }
      
      rho_hat_pairs[i, j, ,] <- rho_ij_mat
      rho_hat_pairs[j, i, ,] <- t(rho_ij_mat)
      
    } # end of j
  } # end of i
  
  return(list(rho_hat = rho_hat, rho_hat_pairs = rho_hat_pairs))
}

estimate_intensities_stratum_parallel <- function(data_all, patient_sel, feature_sel, 
                                         t_seq, ncores) {
  
  # ----------------------------------------------------------------------------
  #
  #  
  # GOAL: estimate the first and second order intensities
  #
  # - using future.apply
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
  
  Sys.setenv(
    OMP_NUM_THREADS = 1,
    OPENBLAS_NUM_THREADS = 1,
    MKL_NUM_THREADS = 1,
    NUMEXPR_NUM_THREADS = 1
  )
  
  plan(multisession, workers = ncores)
  
  NN = length(patient_sel) # 89
  p <- length(feature_sel) # 249
  n_time <- length(t_seq)  # 19
  
  

  

  
  rho_list <- future_lapply(1:p, function(i){
    
    invisible(suppressWarnings(suppressMessages(library(dplyr))))
    invisible(suppressWarnings(suppressMessages(library(data.table))))

    
    data_i = data_all[feature_id==i,]
    data_i_count = data_i[,.(count=.N),by="subject_num"] # xi_i's
    
    # 1a) Estimate first-order intensity (rho_i = 19 dim vector)
    
    if(nrow(data_i) == 0){ # trivial case where process i has no events, it's intensity (rho) is zero
      rho_i <- rep(0, length(t_seq))
    } else{
      
      
      h_i <- get_gamma(data_i$time)
      
      Gamma_i <- data_i[,estimate_density(time, t_seq),by=c("subject_num")]   # vertically stacked \Gamma_i^k's (89 * 19 length matrix)
      
      
      # patient_now = unique(Gamma_i$subject_num)                            # patients left after filtering for process i
      # data_i_count = data_i_count[match(patient_now, subject_num),]        # filter out subjects that are not in "patient_now". 
      
      # density_mat = matrix(Gamma_i$gamma_hat,nrow=length(t_seq))                    # rearrange to be a 19 x 89 matrix
      # 
      # xi_i_mat <- matrix(data_i_count$count, nrow = 1) # (1 x 89 mat) xi's in neuron number order
      # 
      # density_xi_mat <- sweep(density_mat, 2, xi_i_mat, `*`)         # multiply by xi_i (19 x 89 times 1 x 89)
      # 
      # rho_i = apply(density_xi_mat, 1, sum) / NN                             # rho_i(t) for mark i (length 19 vector)
      
      rho_mat <- matrix(Gamma_i$rho_hat, nrow = n_time)
      rho_i <- apply(rho_mat, 1, sum) / NN
      
      
    }

    
    # 1b) Estimate second-order intensities (double loop over time grid)
    
    rho_ij_list <- vector('list', p)
    
    for (j in i:p) {
      
      
      data_j = data_all[feature_id==j,]
      data_j_count = data_j[,.(count=.N),by="subject_num"] # xi_j's  
      
      if(nrow(data_j) == 0 || nrow(data_i) == 0){ # moot case where either process has no events
        rho_ij_mat <- matrix(0, nrow = n_time, ncol = n_time)
      } else{
        times_i <- data_i[, .(event_times_i = list(time)), by = subject_num]    # concatenate all the event times of rep k into a vector, then make dataframe
        times_j <- data_j[, .(event_times_j = list(time)), by = subject_num]
        
        times_ij <- merge(times_i, times_j, by = "subject_num", all = TRUE)
        
        Gamma_ij <- times_ij[, estimate_bivariate_density(event_times_i[[1]], 
                                                          event_times_j[[1]], 
                                                          t_seq, 
                                                          t_seq,
                                                          'i'
        ),  by = 'subject_num']
        
        bivariate_intensity <- matrix(Gamma_ij$V1,nrow=length(t_seq)^2)
        
        rho_ij = apply(bivariate_intensity, 1, sum) / NN               # rho_ij(s, t) = 19 x 19  
        rho_ij_mat <- matrix(rho_ij, nrow = length(t_seq), ncol = length(t_seq))
      }
      
      rho_ij_list[[j]] <- rho_ij_mat
      
    } # end of j
    
    list(rho_i = rho_i, rho_ij_list = rho_ij_list)
    
  }, future.seed = TRUE) 
  
  # reconstruct output --------------------------------------
  
  # First-order intensities: p x m matrix
  rho_hat <- matrix(0, nrow=p, ncol=n_time)
  
  # Second-order intensities: p x p x m x m array  
  rho_hat_pairs <- array(0, dim=c(p, p, n_time, n_time))  
  
  for (i in 1:p) {
    rho_hat[i, ] <- rho_list[[i]]$rho_i
    for (j in i:p) {
      rho_ij_mat <- rho_list[[i]]$rho_ij_list[[j]]
      rho_hat_pairs[i, j, , ] <- rho_ij_mat
      rho_hat_pairs[j, i, , ] <- t(rho_ij_mat)
    }
  }
  
  
  return(list(rho_hat = rho_hat, rho_hat_pairs = rho_hat_pairs))
}

estimate_intensities_stratum_parallel_v2 <- function(data_all, patient_sel, feature_sel, 
                                                     t_seq, ncores) {
  
  # ----------------------------------------------------------------------------
  #
  #  
  # GOAL: estimate the first and second order intensities
  #
  # - using pbmclapply
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
  
  
  NN = length(patient_sel) # 89
  p <- length(feature_sel) # 249
  n_time <- length(t_seq)  # 19
  
  # first order
  rho_i_list <- pbmclapply(1:p, function(i) {
    data_i <- data_all[feature_id == i, ]
    
    if (nrow(data_i) == 0) {
      rho_i <- rep(0, n_time)
    } else {
      h_i <- get_gamma(data_i$time)
      Gamma_i <- data_i[, estimate_density(time, t_seq), by = "subject_num"]
      rho_mat <- matrix(Gamma_i$rho_hat, nrow = n_time)
      rho_i <- apply(rho_mat, 1, sum) / NN
    }
    
    rho_i
  }, mc.cores = ncores)
  
  # second order
  ijs <- data.frame(rep(1:p, each = p),
                    rep(1:p, p))
  colnames(ijs) <- c('i', 'j')
  
  ijs <- ijs %>% filter(j >= i)
  
  rho_ij_list <- pbmclapply(1:nrow(ijs), function(k){
    i <- ijs[k,1]
    j <- ijs[k,2]
    
    data_i <- data_all[feature_id == i, ]
    data_j <- data_all[feature_id == j, ]
    
  
    # fitting
    if (nrow(data_j) == 0 || nrow(data_i) == 0) {
      rho_ij_mat <- matrix(0, nrow = n_time, ncol = n_time)
    } else {
      times_i <- data_i[, .(event_times_i = list(time)), by = subject_num]
      times_j <- data_j[, .(event_times_j = list(time)), by = subject_num]
      times_ij <- merge(times_i, times_j, by = "subject_num", all = TRUE)
      
      Gamma_ij <- times_ij[, estimate_bivariate_density(
        event_times_i[[1]], event_times_j[[1]],
        t_seq, t_seq, 'i'), by = 'subject_num']
      
      bivariate_intensity <- matrix(Gamma_ij$V1, nrow = n_time^2)
      rho_ij <- apply(bivariate_intensity, 1, sum) / NN
      rho_ij_mat <- matrix(rho_ij, nrow = n_time)
    }
    rho_ij_mat # return this
  }, mc.cores = ncores)
  
  # Since they are all lists, we need to assemble them now
  
  rho_hat <- matrix(0, nrow = p, ncol = n_time)
  rho_hat_pairs <- array(0, dim = c(p, p, n_time, n_time))
  
  ij_index <- 1
  
  for (i in 1:p) {
    rho_hat[i, ] <- rho_i_list[[i]]
    for (j in i:p) {
      rho_ij <- rho_ij_list[[ij_index]]
      ij_index <- ij_index + 1
      rho_hat_pairs[i, j, , ] <- rho_ij
      rho_hat_pairs[j, i, , ] <- t(rho_ij)
    }
  }
  
  return(list(rho_hat = rho_hat, rho_hat_pairs = rho_hat_pairs))  
  
}

estimate_intensities_stratum_parallel_v3 <- function(data_all, patient_sel, feature_sel, 
                                                     t_seq, ncores) {
  
  # ----------------------------------------------------------------------------
  #
  #  
  # GOAL: estimate the first and second order intensities
  #
  # - using pbmclapply
  # - rho_i  (p x m matrix)
  # - rho_{ij} (m x m matrix) for all (p x p) pairs
  #
  # - in v3, we only care about i = j entries for bivariate estimates
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
  # - rho_i_list (list of p entries)
  #   - each entry is a list of 2 matrices
  #     - rho_i      (m-dim vector)            univariate intensity
  #     - rho_ii_mat (m x m dim matrix)        bivariate intensity
  #
  # ---------------------------------------------------------------------------
  
  
  NN = length(patient_sel) # 89
  p <- length(feature_sel) # 249
  n_time <- length(t_seq)  # 19
  

  rho_i_list <- pbmclapply(1:p, function(i) {
    data_i <- data_all[feature_id == feature_sel[i], ]
    
    if (nrow(data_i) == 0) {
      rho_i <- rep(0, n_time)
      rho_ii_mat <- matrix(0, nrow = n_time, ncol = n_time)
    } else {
      
      # univariate
      Gamma_i <- data_i[, estimate_density(time, t_seq), by = "subject_num"]
      rho_mat <- matrix(Gamma_i$rho_hat, nrow = n_time)
      rho_i <- apply(rho_mat, 1, sum) / NN
      
      # bivariate
      times_i <- data_i[, .(event_times_i = list(time)), by = subject_num]
      Gamma_ii <- times_i[, estimate_bivariate_density(
        event_times_i[[1]], event_times_i[[1]],
        t_seq, t_seq, 'i'), by = 'subject_num']
      
      bivariate_intensity <- matrix(Gamma_ii$V1, nrow = n_time^2)
      rho_ii <- apply(bivariate_intensity, 1, sum) / NN
      rho_ii_mat <- matrix(rho_ii, nrow = n_time)      
    }
    
    list(rho_i=rho_i,
         rho_ii_mat=rho_ii_mat)
  }, mc.cores = ncores)
  
  
  return(rho_i_list)  
  
}
