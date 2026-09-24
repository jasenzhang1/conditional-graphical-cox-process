

packages <- c("dplyr", 'tidyr',
              "data.table", "pROC", "grid", "abind", "future.apply", "pbmcapply",
              "Matrix", "ggplot2", 
              "ggridges", # for ridgeline plot 
              "reshape2", "gridExtra", "viridis", "MASS", 'igraph',
              'patchwork')



for (package in packages) {
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}

source('functions/00a_matrix_massaging.R')
source('functions/00b_matrix_norms.R')
source('functions/00c_block_matrix_arrange.R')
source('functions/00e_preprocessing.R')
source('functions/01_kernel_estimation.R')
source('functions/01b_greshgorin_kernel_initialization.R')
source('functions/02_intensity_estimation.R')
source('functions/03_covariance_function_estimation.R')
source('functions/03a_subject_covariance_function.R')
source('functions/04_eigendecomposition.R')
source('functions/05_KL_coefficients.R')
source('functions/05a_KL_cross_coeffs.R')
source('functions/06_kernel_matrix.R')
source('functions/09_conditional_correlation.R')
source('functions/10b_GIC.R')
source('functions/10c_GIC_parallel.R')
source('functions/11_graph_estimation.R')
source('functions/12b_full_estimation_procedure_mice.R')
source('functions/12z_full_estimation_function_blocks.R')
source('functions/13_estimation_validation.R')



