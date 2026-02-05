packages <- c("dplyr", "data.table", "pROC", "grid", "abind", "future.apply", "pbmcapply",
              "Matrix", "ggplot2", "reshape2", "gridExtra", "viridis", "MASS", 'igraph',
              'knitr', 'tidyr', 'rmarkdown', 'patchwork')

for (package in packages) {
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}

source('functions/00a_matrix_massaging.R')
source('functions/00b_matrix_norms.R')
source('functions/00c_block_matrix_arrange.R')
source('functions/00e_preprocessing.R')
source('functions/01_kernel_estimation.R')
source('functions/01a_kernel_initialization.R')
source('functions/01b_greshgorin_kernel_initialization.R')
source('functions/02_intensity_estimation.R')
source('functions/03_covariance_function_estimation.R')
source('functions/03a_subject_covariance_function.R')
source('functions/04_eigendecomposition.R')
source('functions/05_KL_coefficients.R')
source('functions/05a_KL_cross_coeffs.R')
source('functions/06_kernel_matrix.R')
source('functions/07_regression_operator.R')
# source('functions/07a_regression_operator_speedups.R')
source('functions/07b_regression_operator_G_ij.R')
source('functions/08_conditional_covariance.R')
source('functions/08a_conditional_covariance_G_ij.R')
source('functions/08z_cross_covariance_bundle.R')
source('functions/09_conditional_correlation.R')
source('functions/10_precision_operator.R')
source('functions/10b_GIC.R')
source('functions/11_graph_estimation.R')
source('functions/12_full_estimation_procedure.R')
source('functions/12a_full_estimation_procedure_G_ij.R')
source('functions/12b_full_estimation_procedure_mice.R')
source('functions/12c_full_estimation_procedure_JASA.R')
source('functions/12z_full_estimation_function_blocks.R')
source('functions/13_estimation_validation.R')



