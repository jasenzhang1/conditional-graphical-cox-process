#!/bin/bash

# ------------------------------------------------------------------------------
#
# Data analysis (Section 6): fit the model to the spike-train data of every mouse
# in both strata, select thresholds with lwGIC, save results, and draw per-mouse
# diagnostic figures.
#
# usage:
#
#   ./scripts_exec/run_data_analysis.sh [config file]      (default: config/analysis_paper.sh)
#
# For every mouse and stratum (movement, VR):
#
#   1. pre-process: 10 s windows, stratum filter, 5-event replicate
#      filter, largest-clique neuron retention, top 50 per region    scripts_middle/1_preprocess_data
#   2. estimate rho_i, rho_ij over all windows                       scripts_middle/2_fit_model
#   3. at every observed week: kernel weights (bandwidth h_Y),
#      operator estimation and lwGIC threshold selection             scripts_middle/3_threshold_selection
#   4. save results                                                  scripts_middle/4_save_results
#        -> mice_results/<experiment>/week_only_bw_<h_Y>_min_<floor>/CPGM_<region>/<ID>_<stratum>_t10.RData
#   5. per-mouse diagnostic figures                                  scripts_middle/9_unpack
#
# Afterwards, draw the manuscript figures with scripts_figures/data_analysis_figures.R
# (make analysis-figures).
#
# Needs the spike-train files in data/ (see data/README.md).
# Logs: script_outputs/mice/<ID>_m<movement>vr<VR>_t<time_scale>.log
#
# ------------------------------------------------------------------------------

cd "$(dirname "$0")/.."  # repository root

CONFIG="${1:-config/analysis_paper.sh}"
if [ ! -f "$CONFIG" ]; then
    echo "Config file not found: $CONFIG" >&2
    exit 1
fi
source "$CONFIG"

PREPROCESS_SCRIPTS="scripts_middle/1_preprocess_data"
FIT_SCRIPTS="scripts_middle/2_fit_model"
THRESHOLD_SCRIPTS="scripts_middle/3_threshold_selection"
SAVE_SCRIPTS="scripts_middle/4_save_results"
UNPACK_SCRIPTS="scripts_middle/9_unpack"

model_type="mice"
min_events=$(echo "$time_scale * $min_freq" | bc | xargs printf "%.0f")
USER="${USER:-$(id -un)}"

for f in data/Brain_Region.RData $(for ID in "${IDs[@]}"; do echo "data/${ID}_t${time_scale}_data.rda"; done); do
    if [ ! -f "$f" ]; then
        echo "Missing data file: $f (see data/README.md)" >&2
        exit 1
    fi
done

# results folders written by script_fit_mice_data_part3.R, one per connectivity floor
bw_string="${y_c_bandwidth#*.}"                  # 0.001 -> 001
results_folders=()
for pct in "${min_connect_pcts[@]}"; do
    pct_string=$(printf "%.2f" "$pct")
    case "$pct_string" in
        0.00) min_string="min_00" ;;
        1.00) min_string="min_100" ;;
        *)    min_string="min_${pct_string#*.}" ;;  # 0.05 -> min_05
    esac
    results_folders+=("mice_results/${experiment_folder}/${y_c_structure}_bw_${bw_string}_${min_string}/${method}_${region}")
done

echo "Config: $CONFIG"
echo "Mice: ${IDs[*]}; strata: ${#movement[@]}; h_Y=$y_c_bandwidth; connectivity floor(s): ${min_connect_pcts[*]}"

mkdir -p script_outputs/mice

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


for ID in "${IDs[@]}"; do

    # Loop over all combinations of movement and VR
    for i in "${!movement[@]}"; do
    

        mov=${movement[i]}
        vr=${VR[i]}        
        
        outfile="script_outputs/mice/${ID}_m${mov}vr${vr}_t${time_scale}.log"   #Tau1_m0vr0_t10.log
        rm -f "$outfile"   # delete old log if it exists
        
        echo "===================================================" | tee -a "$outfile"
        echo "Starting full pipeline for mouse=$ID, region=$region, bandwidth=$y_c_bandwidth, movement=$mov, VR=$vr, time_scale=$time_scale" | tee -a "$outfile"
        echo "Logging to: $outfile" | tee -a "$outfile"
        echo "Start time: $(date)" | tee -a "$outfile"
        echo "===================================================" | tee -a "$outfile"
        echo "" | tee -a "$outfile"
        
        step_1_start_time=$(date +%s)
      
        # -------------------
        # Step 1: Preprocess - divide the dataset by discrete covariates
        # -------------------
        
        echo "[STEP 1] Preprocess dataset..." | tee -a "$outfile"
        
        rm -rf "mice_data/$y_c_structure"   # delete old datasets
        
        # creates mice_data/week_only/Data.RData
        Rscript "$PREPROCESS_SCRIPTS/script_preprocess_mice_data.R" "$ID" "$y_c_structure" "$time_scale" "$method" "$m" "$mov" "$vr" "$region" "$min_events" "$max_processes" "$n_weeks" "$deducted_weeks" >> "$outfile" 2>&1 

        
        step_1_end_time=$(date +%s)
        step_1_runtime=$((step_1_end_time - step_1_start_time))
        
        echo "[STEP 1] Finished" | tee -a "$outfile"
        echo "Total elapsed time: ${step_1_runtime} seconds (~$((step_1_runtime/60)) minutes)." >> "$outfile"
        echo "" | tee -a "$outfile"
        echo "===========================================" >> "$outfile"
        echo "" | tee -a "$outfile"
        
        # -----------------------
        # Step 2: Fit
        # -----------------------
        
        step_2_start_time=$(date +%s)
        echo "[STEP 2] Fitting dataset ..." | tee -a "$outfile"
        echo "" | tee -a "$outfile"
        echo "===================================================" >> "$outfile"
        echo "" | tee -a "$outfile"
    
            
        # ----------------
        # Part 1 - get all outer products before re-weighting by y_c
        # ----------------
        
        echo "[PART 1] Collecting parameters ..." >> "$outfile"
        
        
        # temp_data/simu/part1_hub_block_v2_n_100 
        output=$(Rscript "$FIT_SCRIPTS/script_fit_mice_data_part1.R" \
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
            Rscript "$FIT_SCRIPTS/script_step2_part1_v5.R" "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$k" >> "$outfile" 2>&1 &
        done
        wait
        
        echo "[PART 3] Calculating Rho_ij ..." >> "$outfile"
        echo "" | tee -a "$outfile"
      
        for kl in $(seq 1 "$n_ij"); do
            wait_for_slot
            
            # temp_data/simu/step_2_rho_ij...
            Rscript "$FIT_SCRIPTS/script_step2_part2_v5.R" "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$kl" >> "$outfile" 2>&1 &
        done
        wait  
        
        echo "[PART 4] Merging all Rho_i and Rho_ij ..." >> "$outfile"
        echo "" | tee -a "$outfile"
        
        
        # temp_data/simu/step_2_v5_raw_rho_list...
        Rscript "$FIT_SCRIPTS/script_step2_part3_v5.R" "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_i" "$n_ij" >> "$outfile" 2>&1
        
        step_2_end_time=$(date +%s)
        step_2_runtime=$((step_2_end_time - step_2_start_time))
      
        echo "[STEP 2] Finished" | tee -a "$outfile"
        echo "Total elapsed time: ${step_2_runtime} seconds (~$((step_2_runtime/60)) minutes)." >> "$outfile"
        echo "" | tee -a "$outfile"
        echo "===========================================" >> "$outfile"
        echo "" | tee -a "$outfile"
      
        
        step_3_start_time=$(date +%s)
      
        echo "[PART 5] Downstream estimation for all y_c values ..." >> "$outfile"
        echo "" | tee -a "$outfile"
        
        echo "Number of queries: $n_queries" >> "$outfile"
        echo "" | tee -a "$outfile"
        
        

        # Loop over all subsets of y_c
  
        for j in $(seq 1 "$n_queries"); do
        

            echo "[START] y_c = $j" >> "$outfile"
  
            # gets weights and pads rho_i and rho_ii
            # temp_data/simu/step_2_rho_list...
            # temp_data/simu/part2_...
            Rscript "$FIT_SCRIPTS/script_step2_part4_v5.R" "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" "$y_c_bandwidth" >> "$outfile" 2>&1
            
    
            # estimation until GIC
            # temp_data/simu/part2b...
            Rscript "$THRESHOLD_SCRIPTS/script_fit_mice_data_part2b_before_GIC.R" "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" "$eigen_setting" "$y_c_bandwidth" >> "$outfile" 2>&1
            
            
            # ---------------------------------------------------------
            # GIC_LOCAL PROCEDURE
            # ---------------------------------------------------------
            
            echo "At GIC local" >> "$outfile"

            # PART 1: Precompute tau_c quantiles and get num_k
            # data stored in temp_data/simu/GIC_local_folder/file.name
            output=$(Rscript "$THRESHOLD_SCRIPTS/script_GIC_local_part1.R" \
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
            

                echo "[GIC] y_c = $j, suffix ID = $id_suffix" >> "$outfile"
            
                for ((k=1; k<=num_k; k++)); do
                
                  wait_for_slot
                  # Launch one R process per tau_c. 
                  # Inside this R script, it will loop through all tau_p (l=1...num_l)
                  # data stored in temp_data/simu/GIC_local_folder/file.name
                  Rscript "$THRESHOLD_SCRIPTS/script_GIC_local_part2and3_serial.R" \
                      "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" \
                      "$j" "$id_suffix" "$k" "${min_connect_pcts[@]}" >> "$outfile" 2>&1 &                              

                        
                done
                wait
                echo "[GIC DONE] Processing suffix ID: $id_suffix" >> "$outfile"

            done
            wait
            
            echo "All GIC local tasks complete. Combining." >> "$outfile"
            
            # PART 4: Finalize and combine results
            # temp_data/simu/GIC_local_folder/GIC_final_combined.RData
            Rscript "$THRESHOLD_SCRIPTS/script_GIC_local_part4.R" "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" "${min_connect_pcts[@]}" >> "$outfile" 2>&1

            # ---------------------------------------------------------
            

            
            echo "" | tee -a "$outfile"
            echo "===================================================" >> "$outfile"
            echo "" | tee -a "$outfile"
            echo "Merging GIC local info with part2b" >> "$outfile"
            
            # estimation after GIC
            # merge part2b with GIC_final results
            # temp_data/simu/part3...
            Rscript "$THRESHOLD_SCRIPTS/script_fit_mice_data_part2b_after_GIC.R" "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" "${min_connect_pcts[@]}" >> "$outfile" 2>&1
            
              
            
            
            echo "[END] y_c = $j" >> "$outfile"
                


        
            
        done
        wait 

        

        
        # ----------------
        # Part 2c - joint calculation of tau_c and tau_p across all timepoints                      
        #           
        #           AND/OR
        #
        #         - joint calculation of global tau_c and individual tau_p across all timepoints
        # ----------------
        
        
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
            output=$(Rscript "$THRESHOLD_SCRIPTS/script_GIC_global_part1.R" \
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
                        Rscript "$THRESHOLD_SCRIPTS/script_GIC_global_part2and3_serial.R" \
                            "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" \
                            "$id_suffix" "$k" >> "$outfile" 2>&1 &
                    fi
                    
                    # Run HYBRID if requested
                    wait_for_slot
                    
                    if [[ "$global_thresh_method" == "tau_c" || "$global_thresh_method" == "both" ]]; then
                        Rscript "$THRESHOLD_SCRIPTS/script_GIC_hybrid_part2and3_serial.R" \
                            "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" \
                            "$id_suffix" "$k" >> "$outfile" 2>&1 &
                    fi

                done

            done
            wait
            
            echo "All GIC global/hybrid tasks complete. Combining." >> "$outfile"
            
            # PART 4: Finalize and combine results
            
            if [[ "$global_thresh_method" == "joint" || "$global_thresh_method" == "both" ]]; then
                Rscript "$THRESHOLD_SCRIPTS/script_GIC_global_part4.R" "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_queries" >> "$outfile" 2>&1
            fi
            
            if [[ "$global_thresh_method" == "tau_c" || "$global_thresh_method" == "both" ]]; then
                Rscript "$THRESHOLD_SCRIPTS/script_GIC_hybrid_part4.R" "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_queries" >> "$outfile" 2>&1
            fi
            
            
            

        fi
        
        
        # ----------------
        # Part 2d - get step_12 and step_12b, ROC and edge set after all w_mats have been calculated                     
        # ----------------
        
        echo "At part 2d" >> "$outfile"
        
        
        # updating /part3 with more info
        Rscript "$THRESHOLD_SCRIPTS/script_fit_mice_data_part2d.R" "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" \
                                              "$n_queries" "${min_connect_pcts[@]}" >> "$outfile" 2>&1
              
        
        # ----------------
        # Part 3- when all part 2's are done, do part 3
        # ----------------
        
        
        echo "At part 3" >> "$outfile"
        
        # mice_results/adj_type/CPGM/... .RData
        Rscript "$SAVE_SCRIPTS/script_fit_mice_data_part3.R" "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" \
                                             "$n_queries" "$y_c_bandwidth" "$region" "$experiment_folder" "${min_connect_pcts[@]}" >> "$outfile" 2>&1
        
        echo "" | tee -a "$outfile"
        echo "[DONE] Estimating all y_cs" >> "$outfile"
        
        step_3_end_time=$(date +%s)
        step_3_runtime=$((step_3_end_time - step_3_start_time))
            
        echo "[STEP 3] Finished" | tee -a "$outfile"
        echo "Total elapsed time: ${step_3_runtime} seconds (~$((step_3_runtime/60)) minutes)." >> "$outfile"
        echo "" | tee -a "$outfile"
        echo "===========================================" >> "$outfile"
        echo "" | tee -a "$outfile"
      
        wait
        
        # ----------------
        # Last step - visualize results
        # ----------------
        
        echo "" | tee -a "$outfile"
        echo "===================================================" >> "$outfile"
        echo "" | tee -a "$outfile"
        echo "[STEP 3] Graphing results... " >> "$outfile"
        echo "" | tee -a "$outfile"
        
    
        Rscript "$UNPACK_SCRIPTS/script_unpack_mice_results.R" "$ID" "$time_scale" "$method" "$mov" "$vr" "$eigen_setting" "${results_folders[@]}" >> "$outfile" 2>&1
        
        echo "Deleting Files ..." >> "$outfile"
        
        MICE_FOLDER="temp_data/mice/${ID}_m${mov}vr${vr}_t${time_scale}"
        
        rm -rf "$MICE_FOLDER"
        
        echo "Cleanup complete for $ID, m${mov}vr${vr}, t=$time_scale." >> "$outfile"
        
        # ============================================================
        # END TIMER
        # ============================================================
        
        echo "" | tee -a "$outfile"
        echo "===================================================" >> "$outfile"
        
        
        total_end_time=$(date +%s)
        total_runtime=$((total_end_time - step_1_start_time))
        
        echo "Pipeline finished at: $(date)" >> "$outfile"
        echo "Total runtime: ${total_runtime} seconds (~$((total_runtime/60)) minutes)." >> "$outfile"
        
 
    done  # discrete loop
    
done # mouse loop

echo "Done. Results: ${results_folders[*]}"
echo "Next: make analysis-figures"
