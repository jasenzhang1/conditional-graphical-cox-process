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
  
  truths <- load_file(truth_file_name)
  estimates <- load_file(estimates_file_name)
  
  
  full <- F
  
  # estimates has the hierarchy of step_X --> [[i]] --> 'est'
  # truths    has the hierarchy of step_X --> [[i]] --> 'truth'
  
  # merged    has the same hierarchy and ests and truths are put together
  
  merged <- lapply(names(estimates), function(step_name) {
    est_step <- estimates[[step_name]]
    tru_step <- truths[[step_name]]
    
    if (!is.null(tru_step)) {
      # If both are lists of lists (e.g., step_2), merge elementwise
      if (is.list(est_step[[1]]) && is.list(tru_step[[1]])) {
        mapply(function(e, t) c(e, t), est_step, tru_step, SIMPLIFY = FALSE)
      } else {
        # Otherwise, just combine their contents directly (e.g., step_1)
        c(est_step, tru_step)
      }
    } else {
      # No matching truth: keep as-is
      est_step
    }
  })
  
  names(merged) <- names(estimates)
  
  # ----------------------------------------------------------------------------
  # now, utilize the code in 28c
  # ----------------------------------------------------------------------------
  
  
  g_comparisons <- visualize_over_time(merged, graph_ids, full)
  
  return(g_comparisons)
  

  
  
}
