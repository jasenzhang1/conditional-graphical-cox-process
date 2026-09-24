# ------------------------------------------------------------------------------
# Smaller exploratory run (n_large = 500), used in fall 2026 re-runs.
#
#   ./scripts_exec/run_simulations.sh config/simulation_small.sh
#
#
# Each entry of adj_type_params is "<adj_type> <adj_params...>", i.e. one of the
# six topology-trend combinations (hub / complete / flexible-block-banded
# topology x linear (v2) / jump (j2) trend in y_c).
# ------------------------------------------------------------------------------

adj_type_params=(
                                                        # c1 = 4, c2 = 3 always
                                                        # hub:      c3 < c1 / sqrt(3), c4 < c2 / sqrt(3)
  "hub_block_v2 0 1 4 4 3 0 1.2 0 0.9"                  # linear trend, [0, 1.2] and [0, 0.9]
  "hub_block_j2 0 1 4 0.5 4 3 0 1.2 0 0.9"              # jump at y_c = 0.5
                                                        # complete: c3 < c1 / 3, c4 < c2 / 3
  "complete_block_v2 0 1 4 4 3 0 1.2 0 0.9"
  "complete_block_j2 0 1 4 0.5 4 3 0 1.2 0 0.9"
                                                        # flexible: c3 < c1 / 1.618, c4 < c2 / 1.618
  "flexible_block_banded_v2 0 1 4 4 3 0 1.2 0 0.9"
  "flexible_block_banded_j2 0 1 4 0.5 4 3 0 1.2 0 0.9"
)

n_large=500                 # size of each generated dataset
ns=(100 250 500)            # sample sizes
rep_ids=($(seq 11 50))      # replications 11-50
n_query=6                   # number of y_c query points (the query grid)

p=16                        # number of processes (neurons)
d=2                         # number of basis functions
beta_0=4.7                  # baseline log-intensity
n_group=20                  # subjects generated per parallel job
min_events=10               # regenerate a group if any process has fewer events ...
max_events=5000             # ... or more events than these bounds

y_c_bandwidth=""            # empty = default bandwidth gamma_c
eigen_setting="trig_simple"
min_connect_pcts=(0)        # lwGIC connectivity floor; 0 = no floor
global_thresh_method="neither"
method="CPGM"
beta_truth="F"
X_truth="F"

CGCP_SEED=2025              # base RNG seed; every dataset's seed is derived from it
max_jobs=${MAX_JOBS:-80}    # maximum concurrent R processes
