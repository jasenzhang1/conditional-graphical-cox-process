# ------------------------------------------------------------------------------
#
# Data analysis figures (Section 6 and Supplementary Figures S3-S6)
#
# Compares mice and strata once scripts_exec/run_data_analysis.sh has fitted
# every mouse in both strata.
#
#   Rscript scripts_figures/data_analysis_figures.R [results folder] [time_scale]
#
# default results folder: mice_results/paper/week_only_bw_001_min_05/CPGM_BOTH_100_NORMALIZED
#                         (h_Y = 0.001, 5% connectivity floor, top 50 neurons per region)
#
# inputs:  <results folder>/<ID>_<stratum>_t<time_scale>.RData
# outputs: figures/data_analysis/*.png and *.pdf
#
# ------------------------------------------------------------------------------

source('functions/00_function_wrapper.R')
source('functions/20_simulation_function_wrapper.R')

args <- commandArgs(trailingOnly = TRUE)

results_folder <- if (length(args) >= 1) args[1] else 'mice_results/paper/week_only_bw_001_min_05/CPGM_BOTH_100_NORMALIZED'
time_scale     <- if (length(args) >= 2) as.numeric(args[2]) else 10

discrete_levels <- c('m0vr1', 'm1vr1')   # resting / moving, VR on
out_folder      <- 'figures/data_analysis'
dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)

result_files <- list.files(results_folder, pattern = paste0('_t', time_scale, '\\.RData$'), full.names = TRUE)
if (length(result_files) == 0) stop('No results in ', results_folder, '. Run `make analysis` first.')
print(paste0('Found ', length(result_files), ' mouse x stratum fits in ', results_folder))

# every week between the first and last observed week, across all fits
observed_weeks <- unlist(lapply(result_files, function(f) as.numeric(load_file(f)$y_c_query[, 1])))
all_weeks <- seq(floor(min(observed_weeks)), ceiling(max(observed_weeks)))
print(paste0('Weeks: ', min(all_weeks), '-', max(all_weeks)))


# save a ggplot / patchwork / gtable, or a (nested) list of them, as png and pdf
save_figure <- function(obj, name, width, height) {

  if (is.null(obj)) return(invisible(NULL))

  if (inherits(obj, c('ggplot', 'patchwork', 'gtable', 'grob'))) {
    for (ext in c('png', 'pdf')) {
      ggsave(file.path(out_folder, paste0(name, '.', ext)), plot = obj, width = width, height = height, dpi = 300, limitsize = FALSE)
    }
  } else if (is.list(obj)) {
    for (nm in names(obj)) save_figure(obj[[nm]], paste0(name, '_', nm), width, height)
  }
  invisible(NULL)
}

# draw one figure; report and carry on if it fails
make_figure <- function(name, width, height, expr) {
  print(paste0('Drawing ', name))
  tryCatch(save_figure(expr, name, width, height),
           error = function(e) message('  failed: ', name, ': ', conditionMessage(e)))
}


# ------------------------------------------------------------------------------
# estimated edge sets over weeks, all mice, one panel per stratum

make_figure('edge_sets_all_mice', 30, 15,
            visualize_strata_all_mice_v2(results_folder, all_weeks, time_scale, discrete_levels, 'adj', TRUE)$plots)

# edge stability over weeks (within stratum) and between strata
make_figure('edge_stability', 12, 14,
            plot_edge_stability_combined(results_folder, time_scale, discrete_levels))

# proportion of edges that are HIP-HIP, HIP-EHC, EHC-EHC over weeks
make_figure('edge_regional_proportion', 14, 8,
            plot_edge_regional_proportion_all_mice_v2(results_folder, time_scale, discrete_levels))

# median normalized degree over weeks
make_figure('median_degree', 14, 8,
            plot_median_degree_all_mice_v3(results_folder, time_scale, discrete_levels))

# degree distribution over weeks (ridgeline)
make_figure('degree_ridgeline', 14, 12,
            plot_degree_ridgeline_all_mice(results_folder, time_scale, discrete_levels))

# per-mouse edge sets in both strata
for (ID in unique(sub(paste0('_m[01]vr[01]_t', time_scale, '\\.RData$'), '', basename(result_files)))) {
  make_figure(paste0('edge_sets_', ID), 30, 10,
              visualize_edge_set_one_mouse_all_strata(results_folder, all_weeks, time_scale, discrete_levels, 'adj', TRUE, ID)$plot)
}

print('Done')
