# ------------------------------------------------------------------------------
# 
# GOAL: pre-process the data in data/with_ts to make it look just like "dataset"
#
#
# inputs:
#
# - ID               (string)    mouse name such as 'Tau1'
# - y_c_structure    (string)    "week_only" or "time_and_week"
# - time_scale       (integer)   how many seconds is each replicate? Values may be 1, 2, 5, 10
# - method           (string)    estimation method, "CPGM"
# - m                (integer)   time_grid spacing
# - movement         (integer)   0 (resting) or 1 (moving)
# - VR               (integer)   0 (off) or 1 (on)
# - region           (string)    'HIP', 'EHC', or 'HIP_EHC'
# - min_events       (integer)   what is the minimum number of events in subject's neuron to be included?
# - max_processes    (integer)   how many processes should we truncate? 
#
# 
# outputs:
#
# 
# ------------------------------------------------------------------------------

library(RhpcBLASctl)

# limit threads in BLAS/LAPACK
blas_set_num_threads(1)   # limit BLAS
omp_set_num_threads(1)    # limit OpenMP


source('functions/00_function_wrapper.R')

# 1) load args 
args <- commandArgs(trailingOnly = TRUE)

ID             <- args[1]               # ID <- 'Tau3'
y_c_structure  <- args[2]               # y_c_structure <- "week_only" or "time_and_week"
time_scale     <- as.numeric(args[3])   # time_scale <- 10   (each replicate is 5 seconds)
method         <- args[4]               # method <- 'CPGM'
m              <- as.numeric(args[5])   # m <- 20
movement       <- as.numeric(args[6])
VR             <- as.numeric(args[7])
region         <- args[8]
min_events     <- as.numeric(args[9])
max_processes  <- as.numeric(args[10])
n_weeks        <- as.numeric(args[11])

# ID <- 'Tau3'
# y_c_structure <- 'week_only'
# time_scale <- 10
# method <- 'CPGM'
# m <- 30
# movement <- 0
# VR <- 0
# region <- 'HIP'
# min_events <- 5
# max_processes <- 12
# n_weeks <- 6

# 2) create folder and print settings

new_data_folder <- 'mice_data'
if (!dir.exists(new_data_folder)) dir.create(new_data_folder)  # /mice_data

new_data_folder <- paste0(new_data_folder, "/", y_c_structure)
if (!dir.exists(new_data_folder)) dir.create(new_data_folder)   # /mice_data/week_only OR /mice_data/time_and_week


print(paste0("mouse: ", ID))
print(paste0("y_c_structure: ", y_c_structure))
print(paste0("time_scale: ", time_scale))
print(paste0("estimation method: ", method))
print(paste0("num timepoints: ", m))

# 3) load `dataset`
old_data_folder <- 'data/with_ts'
load(paste0(old_data_folder, '/', ID, '_t', time_scale, '_data.rda'))

# load brain region info
load('data/Brain_Region.RData')

# ------------------------------------------------------------
# 4) pre-processing to make it look just like "dataset"
# ------------------------------------------------------------


time_grid_est <- make_time_grid(m)

dataset_k <- convert_data_for_storage(LGCP_data, df_brain_region, ID, y_c_structure, movement, VR, region, time_scale,
                                      time_grid_est, min_events, n_weeks, max_processes = max_processes, seed = NULL) # 00e

# 5) store
discrete_name <- paste0('m', movement, 'vr', VR)
file_dir <- paste0(new_data_folder, '/', ID, '_', discrete_name, '_t', time_scale, '.RData')
save(dataset_k, file = file_dir)


print('saved dataset')


# ALL WARNINGS
print('all warnings below:')
print(warnings())

