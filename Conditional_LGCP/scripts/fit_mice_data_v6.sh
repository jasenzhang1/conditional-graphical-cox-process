#!/bin/bash

# ------------------------------------------------------------------------------
#
# GOAL: Fit data replicates for mice data
# 
# v5:   (2/18/2026)
#       mimic gen_fit v5 code
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
# - eigen_setting     (string)   only_joint, or trig_and_joint, or trig_simple
# - global_thresh_method (string)  both, joint, tau_c, or neither. (both = global and hybrid)
# 
# ------------------------------------------------------------------------------

cd "$(dirname "$0")/.."   # go one level up (from /scripts to /)


IDs=("WT3" "Tau3" "WT1" "WT2" "Tau1" "Tau2")
movement=(0 0 1 1)
VR=(0 1 0 1)

IDs=("Tau3")
movement=(1)
VR=(1)

y_c_structure="week_only"
method="CPGM"
model_type="mice"  # simu or mice
max_jobs=60

y_c_bandwidth=0.0001 # usually its 0.3, exp(-gamma * y_c_diff^2)
m=30

time_scale=10 
min_freq=0.5

min_events=$((time_scale * min_freq))
max_processes=10000

region="HIP_EHC"   # "HIP", "EHC", or "HIP_EHC"
n_weeks=50
eigen_setting="only_joint"  #only_joint, trig_and_joint, trig_simple
global_thresh_method="neither" # both, joint, tau_c, or neither



function wait_for_slot {
    while true; do
        running=$(pgrep -u "$USER" -x "R|Rscript" | wc -l)
        if (( running < max_jobs )); then
            # Small sleep to allow the process you're about to launch 
            # to actually show up in the process table for the next check.
            sleep 0.1 
            break
        fi
        sleep 1
    done
}

print_bar() {
    local current=$1 total=$2 width=40
    local filled=$(( current * width / total ))
    local empty=$(( width - filled ))
    printf "\r["
    printf "%${filled}s" | tr ' ' '#'
    printf "%${empty}s" | tr ' ' '-'
    printf "] %3d / %3d" "$current" "$total"
}

mkdir -p script_outputs
mkdir -p script_outputs/mice

# -------------------------------------------------
# preprocess
# -------------------------------------------------


# Loop over all mice
for ID in "${IDs[@]}"; do

    # Loop over all combinations of movement and VR
    for i in "${!movement[@]}"; do
    
        mov=${movement[i]}
        vr=${VR[i]}
            
        outfile="script_outputs/mice/${ID}_m${mov}vr${vr}_t${time_scale}.log"   #Tau1_m0vr0_t10.log
        rm -f "$outfile"   # delete old log if it exists
        
        echo "===================================================" | tee -a "$outfile"
        echo "Starting full pipeline for mouse=$ID, region=$region, movement=$mov, VR=$vr, time_scale=$time_scale" | tee -a "$outfile"
        echo "Logging to: $outfile" | tee -a "$outfile"
        echo "Start time: $(date)" | tee -a "$outfile"
        echo "===================================================" | tee -a "$outfile"
        echo "" | tee -a "$outfile"
        
        start_time=$(date +%s)
      
        # -------------------
        # Step 1: Preprocess - divide the dataset by discrete covariates
        # -------------------
        
        echo "[STEP 1] Preprocess dataset..." | tee -a "$outfile"
        
        rm -rf "mice_data/$y_c_structure"   # delete old datasets
        
        wait_for_slot
        Rscript script_preprocess_mice_data.R "$ID" "$y_c_structure" "$time_scale" "$method" "$m" "$mov" "$vr" "$region" "$min_events" "$max_processes" "$n_weeks" >> "$outfile" 2>&1 

        
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
        
        echo "[STEP 2] Fitting dataset ..." | tee -a "$outfile"
        echo "" | tee -a "$outfile"
        echo "===================================================" >> "$outfile"
        echo "" | tee -a "$outfile"
    
            
        # ----------------
        # Part 1 - get all outer products before re-weighting by y_c
        # ----------------
        
        echo "[PART 1] Collecting parameters ..." >> "$outfile"
        
        wait_for_slot
        
        # temp_data/simu/part1_hub_block_v2_n_100 
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
                # temp_data/simu/part2_...
                Rscript script_step2_part4_v5.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" "$y_c_bandwidth" >> "$outfile" 2>&1
                
        
                # estimation until GIC
                # temp_data/simu/part2b...
                Rscript script_fit_mice_data_part2b_before_GIC.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" "$eigen_setting" "$y_c_bandwidth" >> "$outfile" 2>&1
                
                
                # ---------------------------------------------------------
                # GIC_LOCAL PROCEDURE
                # ---------------------------------------------------------
                
                echo "At GIC local" >> "$outfile"
    
                # PART 1: Precompute tau_c quantiles and get num_k
                output=$(Rscript script_GIC_local_part1.R \
                         "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" \
                         2>&1 | tee -a "$outfile")
                

                
                # retrieve variables
                            
                num_k=$(echo "$output" | grep "max_k" | awk -F= '{print $2}')
                num_k=$(echo "$num_k" | xargs)
                
                num_suffixes=$(echo "$output" | grep "num_suffixes" | awk -F= '{print $2}')
                num_suffixes=$(echo "$num_suffixes" | xargs)

                echo "The Max K is: $num_k" >> "$outfile"
                echo "The Number of Suffixes is: $num_suffixes" >> "$outfile"


                # Loop through the IDs found in the map
                for id_suffix in $(seq 1 "$num_suffixes"); do
                
                    wait_for_slot
                    (
                        echo "[GIC] y_c = $j, suffix ID = $id_suffix" >> "$outfile"
                    
                        for ((k=1; k<=num_k; k++)); do
                        
                          wait_for_slot
                          # Launch one R process per tau_c. 
                          # Inside this R script, it will loop through all tau_p (l=1...num_l)
                          Rscript script_GIC_local_part2and3_serial.R \
                              "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" \
                              "$j" "$id_suffix" "$k" >> "$outfile" 2>&1 &                              

                                
                        done
                        wait
                        echo "[GIC DONE] Processing suffix ID: $id_suffix" >> "$outfile"
                    ) &
                done
                wait
                
                echo "All GIC local tasks complete. Combining." >> "$outfile"
                
                # PART 4: Finalize and combine results
                Rscript script_GIC_local_part4.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" >> "$outfile" 2>&1
    
                # ---------------------------------------------------------
                
                # estimation after GIC
                # temp_data/simu/part3...
                
                echo "" | tee -a "$outfile"
                echo "===================================================" >> "$outfile"
                echo "" | tee -a "$outfile"
                echo "Merging GIC local info with part2b" >> "$outfile"
                
                Rscript script_fit_mice_data_part2b_after_GIC.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" >> "$outfile" 2>&1
                
                  
                
                
                echo "[END] y_c = $j" >> "$outfile"
            ) &
        
            
        done
        wait 
        
        # ----------------
        # Part 2c - joint calculation of tau_c and tau_p across all timepoints                      
        #           
        #           AND/OR
        #
        #         - joint calculation of global tau_c and individual tau_p across all timepoints
        # ----------------
        
        wait_for_slot 
        
        # Check for joint or both
        if [[ "$global_thresh_method" == "joint" || "$global_thresh_method" == "both" || "$global_thresh_method" == "tau_c" ]]; then
        
            echo "" | tee -a "$outfile"
            echo "===================================================" >> "$outfile"
            echo "" | tee -a "$outfile"
            echo "Running global/hybrid script..." >> "$outfile"
            echo "" | tee -a "$outfile"
            echo "===================================================" >> "$outfile"
            echo "" | tee -a "$outfile"
            

            # PART 1: Precompute tau_c quantiles and get num_k
            output=$(Rscript script_GIC_global_part1.R \
                     "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_queries" \
                     2>&1 | tee -a "$outfile")
            

            
            # retrieve variables
                        
            num_k=$(echo "$output" | grep "max_k" | awk -F= '{print $2}')
            num_k=$(echo "$num_k" | xargs)
            
            num_suffixes=$(echo "$output" | grep "num_suffixes" | awk -F= '{print $2}')
            num_suffixes=$(echo "$num_suffixes" | xargs)

            echo "The Max K is: $num_k" >> "$outfile"
            echo "The Number of Suffixes is: $num_suffixes" >> "$outfile"
            
            
            # Loop through the IDs found in the map
            for id_suffix in $(seq 1 "$num_suffixes"); do

                echo "suffix="$id_suffix >> "$outfile"
                
                for ((k=1; k<=num_k; k++)); do
                    
                    # Run JOINT if requested
                    wait_for_slot
                    
                    if [[ "$global_thresh_method" == "joint" || "$global_thresh_method" == "both" ]]; then
                        Rscript script_GIC_global_part2and3_serial.R \
                            "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" \
                            "$id_suffix" "$k" >> "$outfile" 2>&1 &
                    fi
                    
                    # Run HYBRID if requested
                    wait_for_slot
                    
                    if [[ "$global_thresh_method" == "tau_c" || "$global_thresh_method" == "both" ]]; then
                        Rscript script_GIC_hybrid_part2and3_serial.R \
                            "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" \
                            "$id_suffix" "$k" >> "$outfile" 2>&1 &
                    fi

                done

            done
            wait
            
            echo "All GIC global/hybrid tasks complete. Combining." >> "$outfile"
            
            # PART 4: Finalize and combine results
            
            if [[ "$global_thresh_method" == "joint" || "$global_thresh_method" == "both" ]]; then
                Rscript script_GIC_global_part4.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_queries" >> "$outfile" 2>&1
            fi
            
            if [[ "$global_thresh_method" == "tau_c" || "$global_thresh_method" == "both" ]]; then
                Rscript script_GIC_hybrid_part4.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_queries" >> "$outfile" 2>&1
            fi
            
            
            

        fi
        
        
        # ----------------
        # Part 2d - get step_12 and step_12b, ROC and edge set after all w_mats have been calculated                     
        # ----------------
        
        echo "At part 2d" >> "$outfile"
        
        wait_for_slot
        
        Rscript script_fit_mice_data_part2d.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_queries" >> "$outfile" 2>&1
              
        
        # ----------------
        # Part 3- when all part 2's are done, do part 3
        # ----------------
        
        wait_for_slot
        
        echo "At part 3" >> "$outfile"
        
        # mice_results/adj_type/CPGM/... .RData
        Rscript script_fit_mice_data_part3.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_queries" >> "$outfile" 2>&1
        
        echo "" | tee -a "$outfile"
        echo "[DONE] Estimating all y_cs" >> "$outfile"
            
      
        wait
        
        # ----------------
        # Last step - visualize results
        # ----------------
        
        echo "" | tee -a "$outfile"
        echo "===================================================" >> "$outfile"
        echo "" | tee -a "$outfile"
        echo "[STEP 3] Graphing results... " >> "$outfile"
        echo "" | tee -a "$outfile"
        
        wait_for_slot
    


        Rscript script_unpack_mice_results.R "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$eigen_setting"  >> "$outfile" 2>&1
        
        echo "Deleting Files ..." >> "$outfile"
        
        MICE_FOLDER="temp_data/mice/${ID}_m${mov}vr${vr}_t${time_scale}"
        
        rm -rf "$MICE_FOLDER"
        
        echo "Cleanup complete for $ID, m${mov}vr${vr}, t=$time_scale." >> "$outfile"
        
        # ============================================================
        # END TIMER
        # ============================================================
        
        echo "" | tee -a "$outfile"
        echo "===================================================" >> "$outfile"
        
        
        end_time=$(date +%s)
        runtime=$((end_time - start_time))
        
        echo "Pipeline finished at: $(date)" >> "$outfile"
        echo "Total runtime: ${runtime} seconds (~$((runtime/60)) minutes)." >> "$outfile"

    done  # discrete strata loop
    wait
    
done  # mouse loop



