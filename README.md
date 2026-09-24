# Conditional Graphical Models for Multivariate Cox Processes

Code to reproduce the simulation study (Section 5) and the analysis of mouse
hippocampal / entorhinal spike trains (Section 6) of

> **TODO: authors, title, journal / arXiv link**

The method estimates a graph among point processes (neurons) whose edges change
with a continuous covariate (here, age in weeks), within strata of a discrete
covariate (moving / resting), and selects thresholds with the locally weighted
GIC (lwGIC).

---

## Quick start

```bash
git clone https://github.com/jasenzhang1/conditional-graphical-cox-process.git
cd conditional-graphical-cox-process

make install       # R packages, pinned in renv.lock
make quick-test    # tiny end-to-end simulation run (about 10 minutes on 4 cores)
```

`make quick-test` generates one small dataset, fits the model, selects thresholds,
saves results to `simu_results/` and draws figures to `figures/simulations/`.
If it finishes without errors, the full runs below will work.

## Reproducing the manuscript

| Command | What it does | Output |
|---|---|---|
| `make simulations` | Section 5: 6 settings × n ∈ {500, 1000, 2000, 4000} × 50 replications, then the simulation figures | `simu_results/`, `figures/simulations/` |
| `make simulation-figures` | redraw the simulation figures from saved results | `figures/simulations/` |
| `make analysis` | Section 6: every mouse in both strata, at every observed week (needs the data, see [Data](#data)) | `mice_results/paper/` |
| `make analysis-figures` | draw the data analysis figures from saved results | `figures/data_analysis/` |
| `make all` | all of the above | |

The full runs take a long time (see [Runtime](#runtime-and-hardware)); start them in the background:

```bash
nohup make simulations > simulations.out 2>&1 &
```

`MAX_JOBS` caps the number of R processes running at once (default 80 for
simulations, 75 for the analysis), e.g. `make simulations MAX_JOBS=32`.
Progress is logged per dataset in `script_outputs/`.

### Figures and tables

<!-- TODO: confirm each row against the manuscript -->

| Manuscript | Produced by | File |
|---|---|---|
| Figure 1 | `scripts_figures/simulation_figures.R` | `figures/simulations/local_accuracy_faceted_2x3.png` **(to confirm)** |
| Figures S1–S2 | `scripts_figures/simulation_figures.R` | `figures/simulations/local_<metric>_faceted_2x3.png`, `<setting>_metrics_CI_plots.pdf`, `<setting>_local_heatmap_combined_n.png` **(to confirm)** |
| Figures 2–3 | `scripts_figures/data_analysis_figures.R` | `figures/data_analysis/edge_sets_all_mice_*.png`, `edge_stability_*.png` **(to confirm)** |
| Figures S3–S6 | `scripts_figures/data_analysis_figures.R` | `figures/data_analysis/edge_regional_proportion.png`, `median_degree.png`, `degree_ridgeline.png`, `edge_sets_<mouse>.png` **(to confirm)** |

Per-setting and per-mouse diagnostics (intensity estimates ρ<sub>i</sub>(t) and
ρ<sub>ij</sub>(s, t), covariance functions, eigenfunctions, Hilbert–Schmidt norms,
selected thresholds) are also written next to each result, in
`simu_results/<setting>/CPGM/rep<i>/export/` and
`mice_results/paper/.../export/`.

## Settings

All settings live in `config/`. Each run script takes a config file as its
only argument.

| Config | Used by | Settings |
|---|---|---|
| `simulation_paper.sh` | `make simulations` | Section 5: six settings, n<sub>large</sub> = 4000, n ∈ {500, 1000, 2000, 4000}, 50 replications, 6 query points, p = 16 |
| `simulation_quick.sh` | `make quick-test` | one setting, n = 40, 1 replication, 3 query points |
| `simulation_small.sh` | | n<sub>large</sub> = 500, n ∈ {100, 250, 500}, replications 11–50 |
| `analysis_paper.sh` | `make analysis` | Section 6: six mice, strata (moving, VR on) and (resting, VR on), 10 s windows, 5-event replicate filter, top 50 neurons per region, h<sub>Y</sub> = 0.001, 5% lwGIC connectivity floor |
| `analysis_sensitivity.sh` | | h<sub>Y</sub> = 0.0003, connectivity floors 1%, 2%, …, 15% |

The six simulation settings cross three topologies (hub, complete,
flexible block-banded) with two trends in y<sub>c</sub> (linear `_v2`,
jump `_j2`); their parameters are listed in `config/simulation_paper.sh`.

### Random seeds

Only the simulated data are random; estimation and the data analysis are
deterministic. Every dataset is generated from a seed derived from one base
seed (`CGCP_SEED=2025` in the config), the setting, the replication, and the
generation group (`simulation_seed()` in `functions/22_generate_random_variables.R`).
Each seed is printed to that dataset's log in `script_outputs/simu/`, so any
single replication can be regenerated on its own.

## Repository layout

```
Makefile                    one-command entry points (run `make` for a list)
config/                     settings for each run (see Settings)
scripts_exec/
  run_simulations.sh        Section 5 pipeline: generate -> fit -> lwGIC -> save -> diagnostics
  run_data_analysis.sh      Section 6 pipeline: preprocess -> fit -> lwGIC -> save -> diagnostics
scripts_middle/             one R script per pipeline step, called by the scripts above
  1_generate_data/          simulated datasets (conditioning variables, log-intensities, events, truths)
  1_preprocess_data/        spike trains -> windows, stratum, replicate and neuron filters
  2_fit_model/              intensities rho_i, rho_ij over all replicates, then kernel weights
  3_threshold_selection/    operator estimation and (lw)GIC threshold selection
  4_save_results/           one result file per sample size / mouse x stratum
  9_unpack/                 per-setting and per-mouse diagnostic figures
scripts_figures/
  simulation_figures.R      manuscript figures for Section 5
  data_analysis_figures.R   manuscript figures for Section 6
functions/                  R functions, loaded by 00_function_wrapper.R and 20_simulation_function_wrapper.R
  00*-13*                   estimation: kernels, intensities, covariance functions, eigendecomposition,
                            KL coefficients, conditional correlation, GIC, graph estimation
  21*-26*                   simulation: graph structures, finite-basis processes, event generation
  28*, 31*-35*, 98*         visualization
data/                       spike-train data (not included, see data/README.md)
install/install_packages.R  package installation (make install)
DESCRIPTION, renv.lock      R dependencies and exact versions
```

Generated during runs (not tracked): `temp_data/` (intermediate files, removed
after each dataset), `simu_data/`, `mice_data/`, `simu_results/`,
`mice_results/`, `script_outputs/` (logs), `figures/`.

### Pipeline

Each step is a separate `Rscript` call, so the steps run in parallel across
processes, neuron pairs, query points and thresholds. For one dataset:

1. **Data.** Simulation: `1_generate_data/` draws the conditioning variables,
   the finite-basis log-intensities and the Cox process events, and computes
   the true graphs at each query point. Data analysis: `1_preprocess_data/`
   applies the filters described in `data/README.md`.
2. **Intensities.** `2_fit_model/script_step2_part1_v5.R` and `part2_v5.R`
   estimate ρ<sub>i</sub> and ρ<sub>ij</sub> for each process and pair once, on
   all replicates; `part3_v5.R` merges them; `part4_v5.R` re-weights them for
   each sample size and query point with the covariate kernel (bandwidth h<sub>Y</sub>).
3. **Operators and lwGIC.** `3_threshold_selection/script_fit_mice_data_part2b_before_GIC.R`
   estimates the conditional covariance and precision operators;
   `script_GIC_local_part1.R` … `part4.R` evaluate lwGIC over the threshold grid
   (with the connectivity floor); `part2b_after_GIC.R` applies the selected thresholds;
   `part2d.R` computes ROC curves and edge sets.
4. **Save.** `4_save_results/script_fit_mice_data_part3.R` writes one `.RData` per
   sample size (simulation) or per mouse × stratum (data analysis).

The `script_GIC_global_*` and `script_GIC_hybrid_*` scripts are alternative
(global / hybrid) threshold rules. They are off in all configs
(`global_thresh_method="neither"`) and not used in the manuscript.

## Requirements

- **R** ≥ 4.2. `renv.lock` pins the package versions of the tested environment
  (R 4.3.3); `make install` restores them with `renv`. Without `renv`, run
  `CGCP_NO_RENV=true make install` to install the packages listed in `DESCRIPTION`
  from CRAN.
- **bash**, **make**, **bc** and **pgrep** (procps). Tested on Linux (Ubuntu 24.04);
  macOS should also work.
- A multi-core machine for the full runs. The scripts start one R process per
  task and keep at most `MAX_JOBS` running.

## Runtime and hardware

<!-- TODO: fill in the machine and times of the manuscript runs -->

| Run | Hardware | Time |
|---|---|---|
| `make quick-test` | 4 vCPU, 16 GB RAM (cloud VM) | about 10 minutes |
| `make simulations` | **TODO** (e.g. N cores, M GB RAM) | **TODO** |
| `make analysis` | **TODO** | **TODO** |

Disk: the simulation results take about **TODO** GB; `temp_data/` holds up to a
few GB per dataset while it is being fitted.

## Data

The spike-train data (Zwang et al., 2025) are not redistributed here.
[`data/README.md`](data/README.md) explains how to obtain them, which files to
place in `data/`, and what the preprocessing does.

## Citation, license and archive

- **License:** MIT, see [`LICENSE`](LICENSE).
- **Citation:** see [`CITATION.cff`](CITATION.cff) (GitHub's "Cite this repository").
- **Archived release:** the version used for the manuscript is archived on Zenodo:
  **TODO: DOI badge / link**.
