# Data

The data analysis (Section 6) uses the in vivo spike-train recordings of
Zwang et al. (2025). The data are not redistributed in this repository.

## Getting the data

<!-- TODO: replace with the actual source and access terms -->
1. Request / download the recordings from: **TODO: data source URL or contact**
   (access terms: **TODO**).
2. Place the following files directly in this folder (`data/`):

| File | Contents |
|---|---|
| `Tau1_t10_data.rda`, `Tau2_t10_data.rda`, `Tau3_t10_data.rda` | spike times of each tau-model mouse, split into 10-second windows; loads an object `LGCP_data` |
| `WT1_t10_data.rda`, `WT2_t10_data.rda`, `WT3_t10_data.rda` | the same for each wild-type mouse |
| `Brain_Region.RData` | brain region (hippocampus / entorhinal cortex) of every recorded neuron; loads `df_brain_region` |

`make analysis` checks that all seven files are present before it starts.

## What the pipeline does with them

`scripts_middle/1_preprocess_data/script_preprocess_mice_data.R` (functions in
`functions/00e_preprocessing.R`), for each mouse and stratum:

1. keeps the 10-second windows (replicates) of the stratum: VR on, and
   movement = 1 (moving) or 0 (resting);
2. pre-selects the 75 most active neurons per region (hippocampus, entorhinal
   cortex) by total spike count;
3. in each stratum separately, repeats until nothing changes:
   (a) drops a (window, neuron) replicate with fewer than 5 spikes, and
   (b) keeps the largest clique of neurons in which every pair shares at least
   one window, so that every pairwise intensity can be estimated;
4. intersects the surviving neurons across the two strata, so both strata use
   the same neurons;
5. keeps the top 50 neurons per region, ranked by spike count over VR-on
   windows (`region = "BOTH_100_NORMALIZED"`);
6. writes `mice_data/week_only/<ID>_<stratum>_t10.RData`, the input to fitting.
