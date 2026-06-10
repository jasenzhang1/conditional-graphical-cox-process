# only interested in plotting estimated graphs over time


source('functions/00_function_wrapper.R')
source('functions/20_simulation_function_wrapper.R')



# ------------------------------------------------------------------------------
# for all 4 discrete settings, group them and plot adjacency matrices over time


args <- commandArgs(trailingOnly = TRUE)
experiment_folder <- args[1]  # "mice_results/experiment_12"
results_folder    <- args[2]  # "mice_results/experiment_12/week_only_bw_001_min_05" etc.

print(experiment_folder)
print(results_folder)

method <- 'CPGM'
region <- 'BOTH_100_NORMALIZED'  # HIP, EHC, BOTH_100, BOTH_150
time_scale <- 10



# base_folder <- '../../../project-biostat-chair/mice_results'


results_folder <- paste0(results_folder, '/', method, '_', region)  # mice_results/experiment_9/week_only_bw_default_min_1/CPGM_BOTH_150

# place "_results" before second "/"
graph_folder <- sub("^([^/]*/[^/]*)", "\\1_results", args[2]) # mice_results/experiment_9_results/week_only_bw_default_min_01/


if (!dir.exists(graph_folder)) {
  dir.create(graph_folder, recursive = TRUE)
}




IDs <- c('Tau1', 'Tau2', 'Tau3', 'WT1', 'WT2', 'WT3')


discrete_levels <- paste0('m', c(0, 1), 'vr', c(1, 1))

discrete_levels_list <- list(discrete_levels,
                             discrete_levels,
                             discrete_levels,
                             discrete_levels)

region_border <- T

# ------------------------------------------------------------------------------
# plot graphs for a single mice, facet wrapping all 4 discrete strata


# for(i in 1:length(IDs)){
#   print(IDs[i])
#   png_name <- paste0(graph_folder, '/', IDs[i], '_t', time_scale, '_edge_sets.png')
#   
#   png(png_name, width = 45, height = 15, units = "in", res = 100)
#   
#   g <- visualize_discrete_comparison(results_folder, IDs[i], time_scale, discrete_levels_list[[i]], 'adj', region_border)
#   print(g)
#   dev.off()
# }

# ------------------------------------------------------------------------------
# Figure 1a: Tau1, weeks 17-22, both strata

result <- visualize_edge_set_one_mouse_all_strata(
  results_folder  = results_folder,
  time_scale      = time_scale,
  discrete_levels = discrete_levels,
  output          = "adj",
  region_border   = TRUE,
  mouse_ID        = "Tau1"
)

# --- save ---

png(file.path(graph_folder, paste0("Fig1a_edge_set_strata.png")),
    width = 4800, height = 600, res = 150)
print(result$plot)
dev.off()


# ------------------------------------------------------------------------------
# plot adjacency matrices for edges, comparing the same discrete strata across all mice

# strata_results <- visualize_strata_all_mice_v2(
#   results_folder  = results_folder,
#   all_weeks       = 17:22,
#   time_scale      = time_scale,
#   discrete_levels = discrete_levels,
#   output          = 'adj',
#   region_border   = region_border
# )
# 
# for (d_level in names(strata_results$plots)) {
#   png_name <- paste0(graph_folder, '/all_mice_', d_level, '_t', time_scale, '_edge_sets.png')
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
#   png_name <- paste0(graph_folder, '/all_mice_', d_level, '_t', time_scale, '_edge_proportion.png')
#   png(png_name, width = 10, height = 6, units = "in", res = 100)
#   print(prop_plots[[d_level]])
#   dev.off()
# }


# ------------------------------------------------------------------------------
# plot jaccard similarity temporally - loess

# instability_plots <- plot_edge_instability_all_mice_v2(
#   results_folder  = results_folder,
#   time_scale      = time_scale,
#   discrete_levels = discrete_levels
# )
# 
# for (scale_type in c("linear", "sqrt")) {
#   png_name <- paste0(graph_folder, '/all_mice_t', time_scale, '_temporal_jaccard_loess_', scale_type, '.png')
#   png(png_name, width = 7, height = 10, units = "in", res = 300)
#   print(instability_plots[[scale_type]])
#   dev.off()
# }

# ------------------------------------------------------------------------------
# plot jaccard similarity across strata - loess

# discrete_levels <- paste0('m', c(0, 1), 'vr', c(1, 1))
# 
# d_level_A <- discrete_levels[1]
# d_level_B <- discrete_levels[2]
# 
# strata_instability_plots <- plot_strata_instability_all_mice_v2(
#   results_folder  = results_folder,
#   time_scale      = time_scale,
#   discrete_levels = discrete_levels
# )
# 
# 
# for (plot_id in names(strata_instability_plots$linear)) {
#   for (scale_type in c("linear", "sqrt")) {
#     png_name <- paste0(graph_folder, '/all_mice_', d_level_A, '_vs_', d_level_B, '_t', time_scale, '_strata_jaccard_loess_', plot_id, '_', scale_type, '.png')
#     png(png_name, width = 14, height = 4, units = "in", res = 300)
#     print(strata_instability_plots[[scale_type]][[plot_id]])
#     dev.off()
#   }
# }

# ------------------------------------------------------------------------------
# plot jaccard similarity temporally - dots and lines

# instability_plots <- plot_edge_instability_all_mice(
#   results_folder  = results_folder,
#   time_scale      = time_scale,
#   discrete_levels = discrete_levels
# )
# for (scale_type in c("linear", "sqrt")) {
#   png_name <- paste0(graph_folder, '/all_mice_t', time_scale, '_temporal_jaccard_lines_', scale_type, '.png')
#   png(png_name, width = 7, height = 10, units = "in", res = 300)
#   print(instability_plots[[scale_type]])
#   dev.off()
# }

# ------------------------------------------------------------------------------
# plot jaccard similarity across strata - dots and lines

# discrete_levels <- paste0('m', c(0, 1), 'vr', c(1, 1))
# 
# d_level_A <- discrete_levels[1]
# d_level_B <- discrete_levels[2]
# 
# strata_instability_plots <- plot_strata_instability_all_mice(
#   results_folder  = results_folder,
#   time_scale      = time_scale,
#   discrete_levels = discrete_levels
# )
# 
# 
# for (plot_id in names(strata_instability_plots$linear)) {
#   for (scale_type in c("linear", "sqrt")) {
#     png_name <- paste0(graph_folder, '/all_mice_', d_level_A, '_vs_', d_level_B, '_t', time_scale, '_strata_jaccard_lines_', plot_id, '_', scale_type, '.png')
#     png(png_name, width = 14, height = 4, units = "in", res = 300)
#     print(strata_instability_plots[[scale_type]][[plot_id]])
#     dev.off()
#   }
# }

# ------------------------------------------------------------------------------
# plot temporal and cross-stratum Jaccard similarity (combined faceted plot)

# discrete_levels_vr1 <- paste0('m', c(0, 1), 'vr', c(1, 1))   # c('m0vr1', 'm1vr1')
# 
# stability_plots <- plot_edge_stability_combined(
#   results_folder  = results_folder,
#   time_scale      = time_scale,
#   discrete_levels = discrete_levels_vr1
# )
# 
# for (scale_type in c("linear", "sqrt")) {
#   png_name <- paste0(graph_folder, '/all_mice_t', time_scale, '_edge_stability_combined_', scale_type, '.png')
#   png(png_name, width = 7, height = 7, units = "in", res = 300)
#   print(stability_plots[[scale_type]])
#   dev.off()
# }

# ------------------------------------------------------------------------------
# plot edge regional proportions across strata

# discrete_levels <- paste0('m', c(0, 1), 'vr', c(1, 1))
# 
# edge_regional_plots <- plot_edge_regional_proportion_all_mice(
#   results_folder  = results_folder,
#   time_scale      = time_scale,
#   discrete_levels = discrete_levels
# )
# 
# for (plot_id in names(edge_regional_plots)) {
#   png_name <- file.path(graph_folder,
#                         paste0('all_mice_', plot_id, '_t', time_scale, '_edge_regional_proportion.png'))
#   png(png_name, width = 14, height = 4, units = "in", res = 300)
#   print(edge_regional_plots[[plot_id]])
#   dev.off()
# }

# v2

# discrete_levels <- paste0('m', c(0, 1), 'vr', c(1, 1))
# edge_regional_plot <- plot_edge_regional_proportion_all_mice_v2(
#   results_folder  = results_folder,
#   time_scale      = time_scale,
#   discrete_levels = discrete_levels
# )
# png_name <- file.path(graph_folder,
#                       paste0('all_mice_t', time_scale, '_edge_regional_proportion_v2.png'))
# png(png_name, width = 14, height = 7, units = "in", res = 300)
# print(edge_regional_plot)
# dev.off()

# ------------------------------------------------------------------------------
# plot median non-zero degree across strata

# discrete_levels <- paste0('m', c(0, 1), 'vr', c(1, 1))
# 
# degree_plots <- plot_median_nonzero_degree_all_mice(
#   results_folder  = results_folder,
#   time_scale      = time_scale,
#   discrete_levels = discrete_levels
# )
# 
# for (plot_id in names(degree_plots)) {
#   png_name <- file.path(graph_folder,
#                         paste0('all_mice_', plot_id, '_t', time_scale, '_median_nonzero_degree.png'))
#   png(png_name, width = 14, height = 4, units = "in", res = 300)
#   print(degree_plots[[plot_id]])
#   dev.off()
# }

# v2

# discrete_levels <- paste0('m', c(0, 1), 'vr', c(1, 1))
# degree_plot <- plot_median_nonzero_degree_all_mice_v2(
#   results_folder  = results_folder,
#   time_scale      = time_scale,
#   discrete_levels = discrete_levels
# )
# png_name <- file.path(graph_folder,
#                       paste0('all_mice_t', time_scale, '_median_nonzero_degree.png'))
# png(png_name, width = 7, height = 5, units = "in", res = 300)
# print(degree_plot)
# dev.off()

# v3 - median degree (including zeros)

# discrete_levels <- paste0('m', c(0, 1), 'vr', c(1, 1))
# degree_plot <- plot_median_degree_all_mice_v3(
#   results_folder  = results_folder,
#   time_scale      = time_scale,
#   discrete_levels = discrete_levels
# )
# png_name <- file.path(graph_folder,
#                       paste0('all_mice_t', time_scale, '_median_degree.png'))
# png(png_name, width = 7, height = 5, units = "in", res = 300)
# print(degree_plot)
# dev.off()

# ------------------------------------------------------------------------------
# ridgeline plot of normalized degree

# discrete_levels <- paste0('m', c(0, 1), 'vr', c(1, 1))
# ridgeline_plot <- plot_degree_ridgeline_all_mice(
#   results_folder  = results_folder,
#   time_scale      = time_scale,
#   discrete_levels = discrete_levels
# )
# png_name <- file.path(graph_folder,
#                       paste0('all_mice_t', time_scale, '_degree_ridgeline.png'))
# png(png_name, width = 9, height = 20, units = "in", res = 300)
# print(ridgeline_plot)
# dev.off()



# ------------------------------------------------------------------------------
# edge set for a single mouse across strata

# result <- visualize_edge_set_one_mouse_all_strata_v2(
#   results_folder  = results_folder,
#   time_scale      = time_scale,
#   discrete_levels = discrete_levels,
#   output          = "adj",
#   region_border   = TRUE,
#   mouse_ID        = "Tau2"
# )
# 
# 
# png(file.path(graph_folder, paste0("Tau2_t", time_scale, "_edge_set_strata.png")),
#     width = 4800, height = 600, res = 150)
# print(result$plot)
# dev.off()

