

packages <- c("dplyr", "tidyr", "data.table", "gridExtra", "grid", "patchwork",
              "Matrix", "mvtnorm", "igraph", "MASS", "abind", "reshape2",
              "ggplot2", "parallel", "purrr", "pracma")

for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

source('functions/21_graph_structure_generation.R')
source('functions/22_generate_random_variables.R')
source('functions/23_conditional_dependence_functions.R')
source('functions/24_log_intensity_generation.R')
source('functions/24a_precision_operator_from_E_yc_yd.R')
source('functions/25_point_generation_process.R')
source('functions/26_complete_data_generation_pipeline.R')
source('functions/27_performance_evaluation_functions.R')
source('functions/28_Simulation_Visualization.R')
source('functions/28b_Simulation_Visualization_2.R')
source('functions/00a_matrix_massaging.R')
source('functions/00b_matrix_norms.R')



