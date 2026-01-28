#!/bin/bash

# ------------------------------------------------------------------------------
#
# GOAL: Fit data replicates for mice data
# 
# v2:   (1/12/2026)
#       Incorporate parallelized estimation code 
#
# 
# inputs:
#
# 
# - IDs               (vector)   mice ID that we want to make data for
# - y_c_structure     (string)   "week_only" or "time_and_week"
# - method            (string)
# - model_type        (string)
# - time_scale        (vector)   how long is each data replicate? 2, 5, or 10 seconds
# - m                 (integer)  granularity of time_grid
# - movement          (integer)  not moving (0) or moving (1)
# - VR                (integer)  VR off (0) or on (1)
# - min_events        (integer)  what is the minimum number of events in subject's neuron to be included?
# - max_processes     (integer)  do we truncate the number of neurons?
# ------------------------------------------------------------------------------

cd "$(dirname "$0")/.."   # go one level up (from /scripts to /)


IDs=("Tau3")
y_c_structure="week_only"
method="CPGM"
model_type="mice"  # simu or mice
max_jobs=30

time_scale=10 
m=30
movement=(0 0 1 1)
VR=(0 1 0 1)
min_events=5
max_processes=1000
n_weeks=5
eigen_setting="only_joint"  #only_joint, trig_and_joint

# function wait_for_slot {
#     # Wait until the number of background jobs is strictly less than max_jobs
#     while true; do
#         running=$(jobs -rp | wc -l)
#         if (( running < max_jobs )); then
#             break
#         fi
#         sleep 0.5
#     done
# }

function wait_for_slot {
    while true; do
        # Count any process named 'R' or 'Rscript' owned by you
        running=$(pgrep -u jasen -x "R" | wc -l)
        if (( running < max_jobs )); then
            break
        fi
        sleep 1
    done
}

mkdir -p script_outputs
mkdir -p script_outputs/mice

# -------------------------------------------------
# preprocess
# -------------------------------------------------


# Loop over all combinations of movement and VR
for ID in "${IDs[@]}"; do

    # Loop over all combinations of movement and VR
    for i in "${!movement[@]}"; do
    
        wait_for_slot
        (
    
            mov=${movement[i]}
            vr=${VR[i]}
                
            outfile="script_outputs/mice/${ID}_m${mov}vr${vr}_t${time_scale}.log"   #Tau1_m0vr0_t10.log
            rm -f "$outfile"   # delete old log if it exists
            
            
            echo "===================================================" | tee -a "$outfile"
            echo "Starting full pipeline for mouse=$ID, movement=$mov, VR=$vr, time_scale=$time_scale" | tee -a "$outfile"
            echo "Logging to: $outfile" | tee -a "$outfile"
            echo "Start time: $(date)" | tee -a "$outfile"
            echo "===================================================" | tee -a "$outfile"
            echo "" | tee -a "$outfile"
            
            start_time=$(date +%s)
          
            # -------------------
            # Step 1: Preprocess - divide the dataset by discrete covariates
            # -------------------
            
            echo "[STEP 1] Preprocess dataset..." | tee -a "$outfile"
            
            wait_for_slot
            Rscript script_preprocess_mice_data.R "$ID" "$y_c_structure" "$time_scale" "$method" "$m" "$mov" "$vr" "$min_events" "$max_processes" "$n_weeks" >> "$outfile" 2>&1 

            
            end_time=$(date +%s)
            runtime=$((end_time - start_time))
            
            echo "[STEP 1] Finished" | tee -a "$outfile"
            echo "Total elapsed time: ${runtime} seconds (~$((runtime/60)) minutes)." >> "$outfile"
            echo "" | tee -a "$outfile"
            echo "===========================================" >> "$outfile"
            echo "" | tee -a "$outfile"
            
            # -----------------------
            # Step 2: Fit
            # -----------------------
            
            echo "[STEP 2] Fitting dataset..." | tee -a "$outfile"
            echo "" | tee -a "$outfile"
            echo "===================================================" >> "$outfile"
            echo "" | tee -a "$outfile"
        
                
            # ----------------
            # Part 1 - get all outer products before re-weighting by y_c
            # ----------------
            
            echo "[PART 1] Collecting parameters ..." >> "$outfile"
            
            wait_for_slot
            output=$(Rscript script_fit_mice_data_part1.R \
                      "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" \
                      2>&1 | tee -a "$outfile")
            
            
            echo "" | tee -a "$outfile"
            echo "===================================================" >> "$outfile"
            echo "" | tee -a "$outfile"
            
            # Extract n_queries, n_i and n_ij from output
      
            n_queries=$(echo "$output" | grep "n_queries" | awk -F= '{print $2}')
            n_queries=$(echo "$n_queries" | xargs)
      
            n_i=$(echo "$output" | grep "n_processes" | awk -F= '{print $2}')
            n_i=$(echo "$n_i" | xargs)
            
            n_ij=$(echo "$output" | grep "n_bivariate_processes" | awk -F= '{print $2}')
            n_ij=$(echo "$n_ij" | xargs)
            
            
            echo "num queries=$n_queries" >> "$outfile"
            echo "num processes=$n_i" >> "$outfile"
            echo "num pairwise processes=$n_ij" >> "$outfile"
            echo "" | tee -a "$outfile"
            echo "===================================================" >> "$outfile"
            echo "" | tee -a "$outfile"
            
            # ----------------
            # Part 2 - getting raw rho_i, rho_ij values without weighing 
            # ----------------
            
            echo "[PART 2] Calculating Rho_i ..." >> "$outfile"
            echo "" | tee -a "$outfile"
          
          
            for k in $(seq 1 "$n_i"); do
                wait_for_slot
                
                # temp_data/simu/step_2_v5_rho_i...
                Rscript script_step2_part1_v5.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$k" >> "$outfile" 2>&1 &
            done
          
          
            echo "[PART 3] Calculating Rho_ij ..." >> "$outfile"
            echo "" | tee -a "$outfile"
          
            for kl in $(seq 1 "$n_ij"); do
                wait_for_slot
                
                # temp_data/simu/step_2_rho_ij...
                Rscript script_step2_part2_v5.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$kl" >> "$outfile" 2>&1 &
            done
            wait  
            
            echo "[PART 4] Merging all Rho_i and Rho_ij ..." >> "$outfile"
            echo "" | tee -a "$outfile"
            
            
            # temp_data/simu/step_2_v5_raw_rho_list...
            wait_for_slot
            Rscript script_step2_part3_v5.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_i" "$n_ij" >> "$outfile" 2>&1
            
          
            echo "[PART 5] Downstream estimation for all y_c values ..." >> "$outfile"
            echo "" | tee -a "$outfile"
            
            echo "Number of queries: $n_queries" >> "$outfile"
            echo "" | tee -a "$outfile"
            
            # Loop over all subsets of y_c
      
            for j in $(seq 1 "$n_queries"); do
                echo "[START] y_c = $j" >> "$outfile"
      
                wait_for_slot
                (
                    # gets weights and pads rho_i and rho_ii
                    # temp_data/simu/step_2_rho_list...
                    Rscript script_step2_part4_v5.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" >> "$outfile" 2>&1
                    
          
                    # estimation from step3 onwards
                    # temp_data/simu/part3...
                    Rscript script_fit_mice_data_part2b.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" "$eigen_setting" >> "$outfile" 2>&1
                    
                    echo "[END] y_c = $j" >> "$outfile"
                ) &
            
                
            done
            wait 
            
            # ----------------
            # Part 3- when all part 2's are done, do part 3
            # ----------------
            
            wait_for_slot
            
            # simu_results/adj_type/CPGM/... .RData
            Rscript script_fit_mice_data_part3.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_queries" >> "$outfile" 2>&1
            
            echo "" | tee -a "$outfile"
            echo "[DONE] Estimating all y_cs" >> "$outfile"
                
          
            wait
            
            # ============================================================
            # END TIMER
            # ============================================================
            
            echo "" | tee -a "$outfile"
            echo "===================================================" >> "$outfile"
            
            
            end_time=$(date +%s)
            runtime=$((end_time - start_time))
            
            echo "Pipeline finished at: $(date)" >> "$outfile"
            echo "Total runtime: ${runtime} seconds (~$((runtime/60)) minutes)." >> "$outfile"
        ) &
    done  # discrete strata loop
done  # mouse loop

wait


