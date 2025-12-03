source('functions/28z_Visualization_helpers.R')
source('functions/28y_Visualization_Blocks.R')
source('functions/28_Simulation_Visualization.R')
source('functions/28c_Grand_Visualizations.R')

visualize_finite_basis <- function(truth_file_name, estimates_file_name, graph_ids){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: visualize intermediate plots of the finite basis simulation setting 
  #
  # - just merging `truths` and `estimates` from the finite basis setting
  # - then using `visualize over time` from 28_c
  #
  # inputs:
  # 
  # - truth_file_name       (string)              'simu_data/block_banded_v2_n_100_truths.RData'
  # - estimates_file_name   (string)              'simu_results/block_banded_v2/CPGM/CPGM_n_100.RData'
  # - graph_ids             (vector of strings)   which graphs do we want to see
  # - i
  # - j
  #
  #
  # outputs:
  #
  # - g_comparisons    (list)  list of graphs that display heatmaps (estimate, truth) over queries
  #
  # ----------------------------------------------------------------------------
  
  # 1) merge truths and estimates
  
  merged <- convergence_metrics_part1(truth_file_name, estimates_file_name)
  

  # 2) now, utilize the code in 28c

  g_comparisons <- visualize_over_time(merged, graph_ids, full)
  
  return(g_comparisons)
  

  
  
}
