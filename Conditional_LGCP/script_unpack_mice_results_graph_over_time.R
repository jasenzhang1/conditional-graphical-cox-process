# only interested in plotting estimated graphs over time


source('functions/00_function_wrapper.R')
source('functions/20_simulation_function_wrapper.R')



# ------------------------------------------------------------------------------
# for all 4 discrete settings, group them and plot adjacency matrices over time


y_c_structure <- "week_only_bw_0003"
method <- 'CPGM'
region <- 'BOTH_100'  # HIP, EHC, BOTH_100
time_scale <- 10

base_folder <- 'mice_results'




results_folder <- paste0(base_folder, '/', y_c_structure, '/', method, '_', region)


results_folder_2 <- paste0(results_folder, '/export')

if (!dir.exists(results_folder_2)) {
  dir.create(results_folder_2)
}


IDs <- c('Tau1', 'Tau2', 'Tau3', 'WT1', 'WT2', 'WT3')


discrete_levels <- paste0('m', c(0, 0, 1, 1), 'vr', c(0, 1, 0, 1))

discrete_levels_list <- list(discrete_levels,
                             discrete_levels,
                             discrete_levels,
                             discrete_levels,
                             discrete_levels,
                             discrete_levels)



for(i in 1:length(IDs)){
  print(IDs[i])
  png_name <- paste0(results_folder_2, '/', IDs[i], '_t', time_scale, '_edge_sets.png')
  
  png(png_name, width = 45, height = 15, units = "in", res = 100)
  
  g <- visualize_discrete_comparison(results_folder, IDs[i], time_scale, discrete_levels_list[[i]], 'adj')
  print(g)
  dev.off()
}


  
