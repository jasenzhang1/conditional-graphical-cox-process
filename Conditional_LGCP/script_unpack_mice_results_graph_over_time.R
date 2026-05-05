# only interested in plotting estimated graphs over time


source('functions/00_function_wrapper.R')
source('functions/20_simulation_function_wrapper.R')



# ------------------------------------------------------------------------------
# for all 4 discrete settings, group them and plot adjacency matrices over time


y_c_structure <- "week_only_bw_003_min_01_before_fix"
method <- 'CPGM'
region <- 'BOTH_150'  # HIP, EHC, BOTH_100, BOTH_150
time_scale <- 10

base_folder <- 'mice_results'

# base_folder <- '../../../project-biostat-chair/mice_results'


results_folder <- paste0(base_folder, '/', y_c_structure, '/', method, '_', region)


results_folder_2 <- paste0(results_folder, '/export')

if (!dir.exists(results_folder_2)) {
  dir.create(results_folder_2)
}


IDs <- c('Tau1', 'Tau2', 'Tau3', 'WT3')


discrete_levels <- paste0('m', c(0, 0, 1, 1), 'vr', c(0, 1, 0, 1))

discrete_levels_list <- list(discrete_levels,
                             discrete_levels,
                             discrete_levels,
                             discrete_levels)

region_border <- T

# ------------------------------------------------------------------------------
# plot graphs for a single mice, facet wrapping all 4 discrete strata


# for(i in 1:length(IDs)){
#   print(IDs[i])
#   png_name <- paste0(results_folder_2, '/', IDs[i], '_t', time_scale, '_edge_sets.png')
#   
#   png(png_name, width = 45, height = 15, units = "in", res = 100)
#   
#   g <- visualize_discrete_comparison(results_folder, IDs[i], time_scale, discrete_levels_list[[i]], 'adj', region_border)
#   print(g)
#   dev.off()
# }

# ------------------------------------------------------------------------------
# plot graphs comparing the same discrete strata across all mice

# strata_results <- visualize_strata_all_mice(
#   results_folder  = results_folder,
#   time_scale      = time_scale,
#   discrete_levels = discrete_levels,
#   output          = 'adj',
#   region_border   = region_border
# )
# 
# for (d_level in names(strata_results$plots)) {
#   png_name <- paste0(results_folder_2, '/all_mice_', d_level, '_t', time_scale, '_edge_sets.png')
#   n_mice   <- strata_results$n_mice[[d_level]]
#   png(png_name, width = 45, height = 4 * n_mice, units = "in", res = 100)
#   print(strata_results$plots[[d_level]])
#   dev.off()
# }

# ------------------------------------------------------------------------------
# plot proportion of connected edges over time with a line graph

# prop_plots <- plot_edge_proportion_all_mice(
#   results_folder  = results_folder,
#   time_scale      = time_scale,
#   discrete_levels = discrete_levels
# )
# 
# for (d_level in names(prop_plots)) {
#   png_name <- paste0(results_folder_2, '/all_mice_', d_level, '_t', time_scale, '_edge_proportion.png')
#   png(png_name, width = 10, height = 6, units = "in", res = 100)
#   print(prop_plots[[d_level]])
#   dev.off()
# }


# ------------------------------------------------------------------------------
# plot edge instability (if they were still present in the next week) over time with a line graph  


instability_plots <- plot_edge_instability_all_mice(
  results_folder  = results_folder,
  time_scale      = time_scale,
  discrete_levels = discrete_levels
)

for (d_level in names(instability_plots)) {
  png_name <- paste0(results_folder_2, '/all_mice_', d_level, '_t', time_scale, '_edge_instability.png')
  png(png_name, width = 14, height = 4, units = "in", res = 100)
  print(instability_plots[[d_level]])
  dev.off()
}

