#!/bin/bash
#$ -cwd
#$ -l h_rt=48:00:00              # walltime
#$ -l h_data=6G                  # memory per job - adjust as needed
#$ -pe shared 18                 # number of cores - match your ncores in R
# Email address to notify
#$ -M $USER@mail #don't change this line, finds your email in the system 
# Notify when
#$ -m bea

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

cd ..   # go one level up (from /scripts to /)

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:

module load apptainer
module load R
MEM_PER_SLOT="2G"
cluster="andrew" # hoffman or andrew
model_type="mice"  # simu or mice

eigen_setting="only_joint"  #only_joint, trig_and_joint, trig_simple
global_thresh_method="neither" # both, joint, tau_c, or neither

y_c_structure="week_only"
method="CPGM"
min_connect_pcts=(0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10)    # 1 percent = 0.01
min_connect_pcts=(0.13 0.25)    # 1 percent = 0.01

m=30
time_scale=10 
min_freq=0.5
min_events=$(echo "$time_scale * $min_freq" | bc | xargs printf "%.0f")


max_processes=10 #10000 
n_weeks=2  #50

if [ "$cluster" == "hoffman" ]; then

    module load R
    ID=$1
    mov=$2
    vr=$3
    y_c_bandwidth=$4
    region=$5
    
    mkdir -p ../../../project-biostat-chair/script_outputs
    mkdir -p ../../../project-biostat-chair/script_outputs/mice
    
else
    IDs=("Tau1" "Tau2" "Tau3" "WT1" "WT2" "WT3")
    movement=(0 1)
    VR=(1 1)
    
    IDs=("Tau1")
    movement=(1)
    VR=(1)
    
    # finished WT3 m0vr0
    
    # EHC, 0.001,  WT2_m0vr1
    # HIP, 0.001,  Tau2_m0vr1, WT1_m0vr0, WT2_m0vr1
    # EHC, 0.0003, WT2_m0vr1
    # HIP, 0.0003, Tau2_m0vr1, WT1_m0vr0, WT2_m0vr1
    y_c_bandwidth='0.1_default' # usually its 0.3, exp(-gamma * y_c_diff^2), use NULL for default
    
    # y_c_bandwidth='1_default'       --> 1 * default
    # y_c_bandwidth='0.1_default'     --> 0.1 * default
    # y_c_bandwidth='0.01_default'    --> 0.01 * default
    # y_c_bandwidth=0.001             --> 0.001
    # y_c_bandwidth=0.0001            --> 0.0001
    

    # I have done 0.001 and 0.0003, next try 0.003
    region="BOTH_150"   # "HIP", "EHC", or "HIP_EHC" "BOTH_100", "BOTH_150"
    
    max_jobs=72
    
    mkdir -p script_outputs
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
    
fi




print_bar() {
    local current=$1 total=$2 width=40
    local filled=$(( current * width / total ))
    local empty=$(( width - filled ))
    printf "\r["
    printf "%${filled}s" | tr ' ' '#'
    printf "%${empty}s" | tr ' ' '-'
    printf "] %3d / %3d" "$current" "$total"
}



# -------------------------------------------------
# preprocess
# -------------------------------------------------

if [ "$cluster" == "hoffman" ]; then

        
    outfile="../../../project-biostat-chair/script_outputs/mice/${ID}_m${mov}vr${vr}_t${time_scale}.log"   #Tau1_m0vr0_t10.log
    rm -f "$outfile"   # delete old log if it exists
    
    echo "===================================================" | tee -a "$outfile"
    echo "Starting full pipeline for mouse=$ID, region=$region, bandwidth=$y_c_bandwidth, movement=$mov, VR=$vr, time_scale=$time_scale" | tee -a "$outfile"
    echo "Logging to: $outfile" | tee -a "$outfile"
    echo "Start time: $(date)" | tee -a "$outfile"
    echo "---------------------------------------------------" | tee -a "$outfile"
    echo "Job ID: $JOB_ID" | tee -a "$outfile"
    echo "Host: $HOSTNAME" | tee -a "$outfile"
    echo "Cores: $NSLOTS" | tee -a "$outfile"
    total_mem=$(echo "$NSLOTS * ${MEM_PER_SLOT%G}" | bc)
    echo "Memory per slot: $MEM_PER_SLOT  =>  Total: ${total_mem}G" | tee -a "$outfile"
    echo "Walltime: 48:00:00" | tee -a "$outfile"
    echo "===================================================" | tee -a "$outfile"
    echo "" | tee -a "$outfile"
    
    step_1_start_time=$(date +%s)
    
    # -------------------
    # Step 1: Preprocess - divide the dataset by discrete covariates
    # -------------------
    
    echo "[STEP 1] Preprocess dataset..." | tee -a "$outfile"
    
    rm -rf "mice_data/$y_c_structure"   # delete old datasets
    
    # creates mice_data/week_only/Data.RData
    Rscript script_preprocess_mice_data.R "$ID" "$y_c_structure" "$time_scale" "$method" "$m" "$mov" "$vr" "$region" "$min_events" "$max_processes" "$n_weeks" "$cluster" >> "$outfile" 2>&1 
    
    
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
    output=$(Rscript script_fit_mice_data_part1.R \
              "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$cluster" \
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
    
    echo "[PART 2+3+4] Calculating Rho_i and Rho_ij in parallel and merging" >> "$outfile"
    echo "" | tee -a "$outfile"

    # it may have aborted
    Rscript script_step2_hoffman.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_i" "$n_ij" >> "$outfile" 2>&1   
    if [ $? -ne 0 ]; then
        echo "[ERROR] script_step2_hoffman.R failed. Aborting." | tee -a "$outfile"
        exit 1
    fi


    step_2_end_time=$(date +%s)
    step_2_runtime=$((step_2_end_time - step_2_start_time))
    
    echo "[STEP 2] Finished" | tee -a "$outfile"
    echo "" | tee -a "$outfile"
    echo "Step 2 elapsed time: ${step_2_runtime} seconds (~$((step_2_runtime/60)) minutes)." >> "$outfile"
    echo "" | tee -a "$outfile"
    echo "===================================================" >> "$outfile"
    echo "" | tee -a "$outfile"

    # on hoffman, run y_c queries in series, parallelizing GIC part 2+3 over k using ncores
    
    step_3_start_time=$(date +%s)
    echo "[PARTS 2-4 + GIC] Running full pipeline for all y_c queries (hoffman) ..." >> "$outfile"
    echo "" | tee -a "$outfile"

    Rscript script_step3_hoffman.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_queries" "$y_c_bandwidth" "$eigen_setting" >> "$outfile" 2>&1

    step_3_end_time=$(date +%s)
    step_3_runtime=$((step_3_end_time - step_3_start_time))
    
    
    echo "[STEP 3] Finished" | tee -a "$outfile"
    echo "" | tee -a "$outfile"
    echo "Step 3 elapsed time: ${step_3_runtime} seconds (~$((step_3_runtime/60)) minutes)." >> "$outfile"
    echo "" | tee -a "$outfile"
    echo "===================================================" >> "$outfile"
    echo "" | tee -a "$outfile"

    # ----------------
    # Part 2d - get step_12 and step_12b, ROC and edge set after all w_mats have been calculated                     
    # ----------------
    
    step_4_start_time=$(date +%s)
    echo "At part 2d" >> "$outfile"
    
    
    # updating /part3 with more info
    Rscript script_fit_mice_data_part2d.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_queries" "$cluster" >> "$outfile" 2>&1
          
    
    # ----------------
    # Part 3- when all part 2's are done, do part 3
    # ----------------
    
    
    echo "At part 3" >> "$outfile"
    
    # mice_results/adj_type/CPGM/... .RData
    Rscript script_fit_mice_data_part3.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_queries" "$y_c_bandwidth" "$region" "$cluster" "$min_connect_pct" >> "$outfile" 2>&1
    
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
    

    Rscript script_unpack_mice_results.R "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$eigen_setting"  >> "$outfile" 2>&1
    
    step_4_end_time=$(date +%s)
    step_4_runtime=$((step_4_end_time - step_4_start_time))
    
    echo "[STEP 4] Finished" | tee -a "$outfile"
    echo "Step 4 elapsed time: ${step_4_runtime} seconds (~$((step_4_runtime/60)) minutes)." >> "$outfile"
    
    # -------------------------
    # Deleting files
    # -------------------------
    
    echo "Deleting Files ..." >> "$outfile"

    MICE_FOLDER="temp_data/mice/${ID}_m${mov}vr${vr}_t${time_scale}"

    rm -rf "$MICE_FOLDER"
    
    echo "Cleanup complete for $ID, m${mov}vr${vr}, t=$time_scale." >> "$outfile"
    
    # ============================================================
    # END TIMER
    # ============================================================
    
    echo "" | tee -a "$outfile"
    echo "===================================================" >> "$outfile"
    
    
    final_end_time=$(date +%s)
    final_runtime=$((final_end_time - step_1_start_time))
    
    echo "Pipeline finished at: $(date)" >> "$outfile"
    echo "Total runtime: ${final_runtime} seconds (~$((final_runtime/60)) minutes)." >> "$outfile"

else

    # Loop over all mice
    for ID in "${IDs[@]}"; do

        # Loop over all combinations of movement and VR
        for i in "${!movement[@]}"; do
        
    # for i in "${!IDs[@]}"; do
    # 
    #     ID=${IDs[i]}
    
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
            Rscript script_preprocess_mice_data.R "$ID" "$y_c_structure" "$time_scale" "$method" "$m" "$mov" "$vr" "$region" "$min_events" "$max_processes" "$n_weeks" "$cluster" >> "$outfile" 2>&1 
    
            
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
            output=$(Rscript script_fit_mice_data_part1.R \
                      "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$cluster" \
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
            wait
            
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
                Rscript script_step2_part4_v5.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" "$y_c_bandwidth" >> "$outfile" 2>&1
                
        
                # estimation until GIC
                # temp_data/simu/part2b...
                Rscript script_fit_mice_data_part2b_before_GIC.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" "$eigen_setting" "$y_c_bandwidth" >> "$outfile" 2>&1
                
                
                # ---------------------------------------------------------
                # GIC_LOCAL PROCEDURE
                # ---------------------------------------------------------
                
                echo "At GIC local" >> "$outfile"
    
                # PART 1: Precompute tau_c quantiles and get num_k
                # data stored in temp_data/simu/GIC_local_folder/file.name
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
                
    
                    echo "[GIC] y_c = $j, suffix ID = $id_suffix" >> "$outfile"
                
                    for ((k=1; k<=num_k; k++)); do
                    
                      wait_for_slot
                      # Launch one R process per tau_c. 
                      # Inside this R script, it will loop through all tau_p (l=1...num_l)
                      # data stored in temp_data/simu/GIC_local_folder/file.name
                      Rscript script_GIC_local_part2and3_serial.R \
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
                Rscript script_GIC_local_part4.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" "${min_connect_pcts[@]}" >> "$outfile" 2>&1
    
                # ---------------------------------------------------------
                
    
                
                echo "" | tee -a "$outfile"
                echo "===================================================" >> "$outfile"
                echo "" | tee -a "$outfile"
                echo "Merging GIC local info with part2b" >> "$outfile"
                
                # estimation after GIC
                # merge part2b with GIC_final results
                # temp_data/simu/part3...
                Rscript script_fit_mice_data_part2b_after_GIC.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" "${min_connect_pcts[@]}" >> "$outfile" 2>&1
                
                  
                
                
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
            
            
            # updating /part3 with more info
            Rscript script_fit_mice_data_part2d.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" \
                                                  "$n_queries" "$cluster" "${min_connect_pcts[@]}" >> "$outfile" 2>&1
                  
            
            # ----------------
            # Part 3- when all part 2's are done, do part 3
            # ----------------
            
            
            echo "At part 3" >> "$outfile"
            
            # mice_results/adj_type/CPGM/... .RData
            Rscript script_fit_mice_data_part3.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" \
                                                 "$n_queries" "$y_c_bandwidth" "$region" "$cluster" "${min_connect_pcts[@]}" >> "$outfile" 2>&1
            
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
            
            
            total_end_time=$(date +%s)
            total_runtime=$((total_end_time - step_1_start_time))
            
            echo "Pipeline finished at: $(date)" >> "$outfile"
            echo "Total runtime: ${total_runtime} seconds (~$((total_runtime/60)) minutes)." >> "$outfile"
            
     
        done  # discrete loop
        
    done # mouse loop

fi

