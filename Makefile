# ------------------------------------------------------------------------------
# One-command entry points. Run `make` to list them.
#
# Long runs are best started in the background, e.g.
#   nohup make simulations > simulations.out 2>&1 &
#
# Set MAX_JOBS to cap the number of concurrent R processes, e.g.
#   make simulations MAX_JOBS=32
# ------------------------------------------------------------------------------

SIM_CONFIG      ?= config/simulation_paper.sh
ANALYSIS_CONFIG ?= config/analysis_paper.sh
MAX_JOBS        ?=
export MAX_JOBS

.PHONY: help install quick-test simulations simulation-figures analysis analysis-figures all clean-temp

help:
	@echo "make install             install R packages (renv.lock)"
	@echo "make quick-test          smoke test: one simulation setting, tiny n (~5-10 min on 4 cores)"
	@echo "make simulations         Section 5: all settings, fits, and Figure 1 / S1-S2"
	@echo "make simulation-figures  redraw the Section 5 figures from saved results"
	@echo "make analysis            Section 6: fit every mouse in both strata (needs data/, see data/README.md)"
	@echo "make analysis-figures    draw the Section 6 figures (Figures 2-3, S3-S6) from saved results"
	@echo "make all                 simulations + analysis + all figures"
	@echo ""
	@echo "Settings: SIM_CONFIG=$(SIM_CONFIG)  ANALYSIS_CONFIG=$(ANALYSIS_CONFIG)"

install:
	Rscript install/install_packages.R

quick-test:
	./scripts_exec/run_simulations.sh config/simulation_quick.sh

simulations:
	./scripts_exec/run_simulations.sh $(SIM_CONFIG)

simulation-figures:
	bash -c 'source $(SIM_CONFIG) && Rscript scripts_figures/simulation_figures.R $$n_large $${rep_ids[$${#rep_ids[@]}-1]} $$method'

analysis:
	./scripts_exec/run_data_analysis.sh $(ANALYSIS_CONFIG)

analysis-figures:
	bash -c 'source $(ANALYSIS_CONFIG) && Rscript scripts_figures/data_analysis_figures.R mice_results/$$experiment_folder/$${y_c_structure}_bw_$${y_c_bandwidth#*.}_min_$$(printf "%.2f" $${min_connect_pcts[0]} | cut -d. -f2)/$${method}_$$region $$time_scale'

all: simulations analysis analysis-figures

clean-temp:
	rm -rf temp_data
