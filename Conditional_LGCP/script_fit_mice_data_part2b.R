t0 <- Sys.time()
library(RhpcBLASctl)

# limit threads in BLAS/LAPACK
blas_set_num_threads(1)   # limit BLAS
omp_set_num_threads(1)    # limit OpenMP


source('functions/00_function_wrapper.R')

# args 
args <- commandArgs(trailingOnly = TRUE)

ID <- args[1]                       # ID <- 'Tau1'
y_c_structure <- args[2]            # y_c_structure <- "week_only" or "time_and_week"
time_scale <- as.numeric(args[3])   # time_scale <- 10   (each replicate is 5 seconds)
method <- args[4]                   # method <- 'CPGM'
movement <- as.numeric(args[5])     # movement <- 0
VR <- as.numeric(args[6])           # VR <- 0
cont_ind <- as.numeric(args[7])

discrete_level <- paste0('m', movement, 'vr', VR)

setting_info_list <- list(ID = ID,
                          y_c_structure = y_c_structure,
                          time_scale = time_scale,
                          method = method,
                          discrete_level = discrete_level)

temp_file_dir <- 'temp_data'

# ---------------------------
# estimation
# ---------------------------

full_conditional_estimation_with_no_truth_part2b(temp_file_dir, setting_info_list, cont_ind)



t1 <- Sys.time()

print(paste0('Time to finish: ', round(as.numeric(t1 - t0, units = "mins"), 2), ' minutes'))  
print(strrep("-", 50))


# ALL WARNINGS
print('all warnings below:')
print(warnings())

