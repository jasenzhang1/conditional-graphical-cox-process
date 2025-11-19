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
n_keys_univariate <- as.numeric(args[7])          
n_keys_bivariate  <- as.numeric(args[8])

discrete_level <- paste0('m', movement, 'vr', VR)

setting_info_list <- list(ID = ID,
                          y_c_structure = y_c_structure,
                          time_scale = time_scale,
                          method = method,
                          discrete_level = discrete_level)



temp_file_dir <- 'temp_data'
if (!dir.exists(temp_file_dir)) dir.create(temp_file_dir)  # /temp_data


# ---------------------------
# estimation
# ---------------------------

if(method == 'CPGM'){
  estimate_intensities_stratum_parallel_with_yc_part3(temp_file_dir, setting_info_list, n_keys_univariate, n_keys_bivariate)
} else{
  stop('Invalid method. Must be CPGM')
}







