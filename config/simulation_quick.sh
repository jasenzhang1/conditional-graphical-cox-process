# ------------------------------------------------------------------------------
# Smoke test: one setting, one replication, tiny n (about 5-10 minutes on 4 cores).
#
#   make quick-test
#   ./scripts_exec/run_simulations.sh config/simulation_quick.sh
#
# Each entry of adj_type_params is "<adj_type> <adj_params...>", i.e. one of the
# six topology-trend combinations (hub / complete / flexible-block-banded
# topology x linear (v2) / jump (j2) trend in y_c).
# ------------------------------------------------------------------------------

adj_type_params=(
                                                        # c1 = 4, c2 = 3 always
                                                        # hub:      c3 < c1 / sqrt(3), c4 < c2 / sqrt(3)
  "hub_block_v2 0 1 4 4 3 0 1.2 0 0.9"                  # linear trend, [0, 1.2] and [0, 0.9]
)

n_large=40                  # size of each generated dataset
ns=(40)                     # sample sizes
rep_ids=(1)                 # one replication
n_query=3                   # number of y_c query points

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
max_jobs=${MAX_JOBS:-4}     # maximum concurrent R processes
