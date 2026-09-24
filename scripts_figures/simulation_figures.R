# ------------------------------------------------------------------------------
#
# Simulation figures (Section 5 and Supplementary Figures S1-S2)
#
# Compares all six topology-trend settings once every setting has finished.
# Called at the end of scripts_exec/run_simulations.sh, or on its own:
#
#   Rscript scripts_figures/simulation_figures.R <n_large> <max_rep> [method]
#
# inputs (from scripts_exec/run_simulations.sh):
#
# - simu_results/<adj_type>/<method>/rep<i>/<adj_type>_n_<n>_rep_<i>.RData
# - simu_data/<adj_type>_n_<n_large>_rep_<i>/<adj_type>_n_<n_large>_rep_<i>_truths.RData
#
# outputs:
#
# - figures/simulations/   copies of every figure below, one folder to look in
# - simu_results/          the same figures where the plotting functions write them
#
# ------------------------------------------------------------------------------

source('functions/00_function_wrapper.R')
source('functions/20_simulation_function_wrapper.R')

args <- commandArgs(trailingOnly = TRUE)

n_large <- as.numeric(args[1])            # n_large <- 4000
n_reps  <- as.numeric(args[2])            # largest replication index, n_reps <- 50
method  <- if (length(args) >= 3) args[3] else 'CPGM'

base_folder <- 'simu_results'
data_root   <- 'simu_data'
out_folder  <- 'figures/simulations'
dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)

# the six settings, laid out as rows (trend) x columns (topology)
row_names <- c('Linear', 'Jump')
col_names <- c('Banded', 'Hub', 'Complete')

adj_grid <- list(c('flexible_block_banded_v2', 'hub_block_v2', 'complete_block_v2'),
                 c('flexible_block_banded_j2', 'hub_block_j2', 'complete_block_j2'))
vert_dashed_line <- list(c(F, F, F),
                         c(T, T, T))  # the jump settings get a dashed line at the jump

adj_types <- unlist(adj_grid)
verts     <- unlist(vert_dashed_line)

settings_done <- dir.exists(file.path(base_folder, adj_types, method))
if (!any(settings_done)) stop('No results found in ', base_folder, '/<adj_type>/', method)
if (!all(settings_done)) message('Missing settings (skipped): ', paste(adj_types[!settings_done], collapse = ', '))

# draw one figure; report and carry on if it fails (e.g. a metric that is
# undefined in every replication)
try_figure <- function(name, expr) {
  tryCatch({ expr; TRUE }, error = function(e) { message('  failed: ', name, ': ', conditionMessage(e)); FALSE })
}

# first replication with saved results for a setting
first_rep <- function(adj_type) {
  rep_dirs <- list.dirs(file.path(base_folder, adj_type, method), full.names = FALSE, recursive = FALSE)
  rep_nums <- as.numeric(sub('^rep', '', grep('^rep[0-9]+$', rep_dirs, value = TRUE)))
  if (length(rep_nums) == 0) return(NA)
  min(rep_nums)
}


# ------------------------------------------------------------------------------
# 1) accuracy / sensitivity / specificity / PPV / NPV / F1 across y_c,
#    loess curves with bands over replications, 2 x 3 grid of settings

print('Faceted metric curves across y_c')

metrics      <- c('accuracy', 'sensitivity', 'specificity', 'ppv', 'npv', 'f1_score')
metric_names <- c('Accuracy', 'Sensitivity', 'Specificity', 'Positive Predictive Value (PPV)', 'Negative Predictive Value (NPV)', 'F1 Score')

results_folders <- lapply(adj_grid, function(x) paste0(base_folder, '/', x, '/', method))

visualize_retrieve_metrics(base_folder, results_folders, metrics, n_reps, row_names, col_names)

for (i in seq_along(metrics)) {
  ok <- try_figure(metrics[i], visualize_metric_CI_across_yc_faceted(metrics[i], metric_names[i], 'local', base_folder, fig_title = NULL))
  if (ok) file.copy(file.path(base_folder, paste0('local_', metrics[i], '_faceted_2x3.png')), out_folder, overwrite = TRUE)
}


# ------------------------------------------------------------------------------
# 2) per setting: metric confidence intervals for each n, and accuracy across
#    y_c for each n

print('Per-setting metric intervals')

for (i in seq_along(adj_types)[settings_done]) {
  results_folder <- file.path(base_folder, adj_types[i], method)

  try_figure(paste(adj_types[i], 'metric intervals'), visualize_metrics_CI(results_folder, n_reps, adj_types[i]))
  try_figure(paste(adj_types[i], 'accuracy across y_c'), visualize_accuracy_CI_across_yc(results_folder, 'local', n_reps, adj_types[i], verts[i]))

  pngs <- list.files(results_folder, pattern = paste0('^', adj_types[i], '_.*\\.(png|pdf)$'), full.names = TRUE)
  file.copy(pngs, out_folder, overwrite = TRUE)
}


# ------------------------------------------------------------------------------
# 3) estimated vs true graph across y_c and n (heatmaps), one replication each

print('Accuracy heatmaps')

for (adj_type in adj_types[settings_done]) {
  rep_i <- first_rep(adj_type)
  if (is.na(rep_i)) next

  results_folder <- file.path(base_folder, adj_type, method, paste0('rep', rep_i))
  truth_file <- paste0(data_root, '/', adj_type, '_n_', n_large, '_rep_', rep_i, '/', adj_type, '_n_', n_large, '_rep_', rep_i, '_truths.RData')
  if (!file.exists(truth_file)) next

  try_figure(paste(adj_type, 'heatmap'), visualize_accuracy_heatmap_across_yc(results_folder, truth_file, adj_type))

  heatmaps <- list.files(file.path(results_folder, 'export'), pattern = '_heatmap_combined_n\\.png$', full.names = TRUE)
  file.copy(heatmaps, out_folder, overwrite = TRUE)
}

print('Done')
