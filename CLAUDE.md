# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

R + bash code reproducing the simulation study (Section 5) and mouse spike-train analysis (Section 6) of a manuscript on conditional graphical models for multivariate Cox processes. It estimates a graph among point processes (neurons) whose edges vary with a continuous covariate y_c (age in weeks) within strata of a discrete covariate (moving/resting, VR), with thresholds chosen by locally weighted GIC (lwGIC). It is not an R package: there is no test suite, lint setup or build step. `DESCRIPTION` only lists dependencies.

## Commands

```bash
make install            # restore packages from renv.lock (CGCP_NO_RENV=true make install -> CRAN from DESCRIPTION)
make quick-test         # end-to-end smoke test (config/simulation_quick.sh, ~10 min on 4 cores); the closest thing to a test
make simulations        # Section 5 (config/simulation_paper.sh), then figures
make simulation-figures # redraw Section 5 figures from simu_results/
make analysis           # Section 6 (config/analysis_paper.sh); needs data/ (see data/README.md)
make analysis-figures
make clean-temp         # rm -rf temp_data
make simulations MAX_JOBS=32 SIM_CONFIG=config/simulation_small.sh   # override concurrency / config
```

## Comments

Each function in the /functions folder must have this docstring format:

  # ----------------------------------------------------------------------------
  #
  # GOAL: describe the goal of this function
  #
  #
  # input:
  #
  # - arg1           (p x p matrix)         description1 
  # - arg2           (string)               description2
  # - arg3           (number)               description3
  #
  #
  # output:
  #
  # - output1        (datatype)             desecription
  #
  # ----------------------------------------------------------------------------
