# ------------------------------------------------------------------------------
# Dummy run: Section 6 settings for one mouse (Tau1) at 5 evenly spaced weeks.
#
#   make analysis ANALYSIS_CONFIG=config/analysis_dummy.sh
#   ./scripts_exec/run_data_analysis.sh config/analysis_dummy.sh
#
# Needs the spike-train files in data/ (see data/README.md).
# ------------------------------------------------------------------------------

IDs=("Tau1")                                   # one mouse
movement=(1 0)                                  # strata, paired with VR below:
VR=(1 1)                                        #   m1vr1 = moving, VR on; m0vr1 = resting, VR on

time_scale=10                   # window (replicate) length in seconds
m=30                            # time grid points per window
min_freq=0.5                    # replicate filter: >= time_scale * min_freq = 5 events
n_weeks=5                       # 5 query weeks, evenly spaced over the observed range
deducted_weeks=''               # weeks to exclude, comma separated (none)
max_processes=10000             # no cap on the number of neurons
region="BOTH_100_NORMALIZED"    # largest-clique retention, then top 50 neurons per region

y_c_structure="week_only"       # conditioning variable: age in weeks
y_c_bandwidth=0.001             # h_Y
min_connect_pcts=(0.05)         # lwGIC connectivity floor: 5% of possible edges
eigen_setting="only_joint"
global_thresh_method="neither"
method="CPGM"

experiment_folder="dummy"       # results go to mice_results/<experiment_folder>/...
max_jobs=${MAX_JOBS:-75}        # maximum concurrent R processes
