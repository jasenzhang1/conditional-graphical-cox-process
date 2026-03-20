#!/bin/bash

# ---------------------------------------------------------------------------
# Generate AND fit finite basis data — single unified log per dataset
# 
# 11/19/2025 - v2: parallelize bivariate estimation
#            - v3: parallelize data generation
# 
#            - v5: calculate rho_i and rho_ij before normalizing them by different n and y_c
#            - v6: parallelize GIC
# ---------------------------------------------------------------------------





cd "$(dirname "$0")/.."  # go one level up (from /scripts to /)

# c1 = c2 = 2
# c3 = c4 = 0.8
# mu_0 = 5.5

adj_type_params=(
  #"indep_c2 0 1"
  #"single_c2 0 1 0.5"
  #"single_v2 0 1 0.7"
  #"single_j2 0 1 0.5 0.3 0.7"
  #"banded_c2 0 1 0.5"
  #"banded_trig2 0 1 0.3"
  #"sparse_v2 0 1 2 0.3 -1 2 0.01"
  #"block_banded_v2 0 1 0.4 0.8 2"
  #"block_banded_c0 0.5 0.5 2"
  
                                                      # c1 = 4, c2 = 3 always
                                                      # c3 < c1 / sqrt(3)
                                                      # c4 < c2 / sqrt(3)
  "hub_block_v2 0 1 4 4 3 0 1.2 0 0.9"                # for hub, we can allow [0, 1.2] and [0, 0.9]
  "hub_block_j2 0 1 4 0.5 4 3 0 1.2 0 0.9"            # jump                  [0, 1.2] and [0, 0.9]
  #"hub_block_c2 0 1 4 2 2 0.7 0.7"
  #"hub_block_c0 0.5 4 2 2 0.7 0.7"
  
  
  
  "complete_block_v2 0 1 4 4 3 0 1.2 0 0.9"             # for complete, c3 < c1/3, c4 < c2/3 [0, 1.2] and [0, 0.9]
  "complete_block_j2 0 1 4 0.5 4 3 0 1.2 0 0.9"         # for jump                           [0, 1.2] and [0, 0.9]
  #"complete_block_c0 0.5 4 2 2 0.7 0.7"
  #"complete_block_c2 0 1 4 2 2 0.7 0.7"
  
                                                         # c3 < c1 / 1.618
                                                         # c4 < c2 / 1.618
  "flexible_block_banded_v2 0 1 4 4 3 0 1.2 0 0.9"       # for flexible, a bit more chill [0, 1.2] and [0, 0.9]   
  "flexible_block_banded_j2 0 1 4 0.5 4 3 0 1.2 0 0.9"   # jump from                      [0, 1.2] and [0, 0.9]
  #"flexible_block_banded_c2 0 1 4 2 2 0.7 0.7"
  #"flexible_block_banded_c0 0.5 4 2 2 0.7 0.7"  
  
)

n_large=4000
ns=(500 1000 2000 4000)

n_large=100
ns=(100)

n_reps=1 # 50
y_c_bandwidth=""
#n_large=20000
#ns=(100 250 500 1000 2500 5000 10000 20000)

n_group=10
groups=$(( n_large / n_group ))

method="CPGM"
model_type="simu"  # simu or mice
max_jobs=80
min_events=10
max_events=5000

p=16 # 16
d=2
n_query=6
beta_0=4.7
beta_truth="F"
X_truth="F"
eigen_setting="trig_simple" #only_joint, trig_and_joint, trig_simple, mfpca
global_thresh_method="neither" #both, joint, tau_c, neither   both = do joint and tau_c



function wait_for_slot {
    while true; do
        running=$(pgrep -u "$USER" -x "R|Rscript" | wc -l)
        if (( running < max_jobs )); then
            # Small sleep to allow the process you're about to launch 
            # to actually show up in the process table for the next check.
            sleep 0.5 
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
mkdir -p script_outputs/simu



# ---------------------------------------------------------------------------
# Generate data + Fit model (single unified log file per adj_type)
# ---------------------------------------------------------------------------

for entry in "${adj_type_params[@]}"; do

    wait_for_slot
    (

        read -r -a fields <<< "$entry"
      
        adj_type="${fields[0]}"         # first field
        adj_params=("${fields[@]:1}")   # all fields after the first
        
        for rep_i in $(seq 1 "$n_reps"); do
    
          outfile="script_outputs/simu/${adj_type}_n_${n_large}_rep_${rep_i}.log"
          rm -f "$outfile"   # delete old log if it exists
        
          echo "===================================================" | tee -a "$outfile"
          echo "Starting full pipeline for adj_type=$adj_type, n=$n_large" | tee -a "$outfile"
          echo "Logging to: $outfile" | tee -a "$outfile"
          echo "Start time: $(date)" | tee -a "$outfile"
          echo "===================================================" | tee -a "$outfile"
          echo "" | tee -a "$outfile"
        
          start_time=$(date +%s)
        
          # -------------------
          # Step 1: Generate
          # -------------------
          echo "[STEP 1] Generating dataset..." | tee -a "$outfile"
          step1_start=$(date +%s)
        
          # temp_data/simu_data/part0...
          Rscript script_generate_finite_basis_data_part0.R "$p" "$d" "$n_large" "$rep_i" "$n_query" "$beta_0" "$adj_type" "${adj_params[@]}" >> "$outfile" 2>&1
          
          for group_idx in $(seq 1 "$groups"); do
              print_bar "$group_idx" "$groups"   # ← live bar on screen
              echo "[ $(date '+%F %T') ] Starting group $group_idx / $groups" >> "$outfile"
              wait_for_slot
              (
                  # temp_data/simu_data/parts1_and_2...
                  Rscript script_generate_finite_basis_data_parts_1_and_2.R "$n_large" "$rep_i" "$adj_type" "$group_idx" "$n_group" "$min_events" "$max_events" >> "$outfile" 2>&1
                  
                  echo "[ $(date '+%F %T') ] Finished group $group_idx" >> "$outfile"
              ) &
              
          done  
          wait
          
          echo "[STEP 1] Finished generating events" | tee -a "$outfile"
          echo "" | tee -a "$outfile"
          echo "===================================================" >> "$outfile"
          echo "" | tee -a "$outfile"
          
          
          # temp_data/simu_data/dataset...
          echo "[STEP 1] Starting merging events" | tee -a "$outfile"
          echo "" | tee -a "$outfile"
          
          Rscript script_generate_finite_basis_data_part3.R "$n_large" "$rep_i" "$adj_type" "$groups" >> "$outfile" 2>&1
          
          echo "[STEP 1] Finished merging events" | tee -a "$outfile"
          echo "" | tee -a "$outfile"
          echo "===================================================" >> "$outfile"
          echo "" | tee -a "$outfile"
          
          # ---------------------------------------------------
          # Step 1b: Get truths for each cont_ind query
          # ---------------------------------------------------
          
          echo "[STEP 1] Getting truths" | tee -a "$outfile"
          echo "" | tee -a "$outfile"
          
          for cont_ind in $(seq 1 "$n_query"); do
              print_bar "$cont_ind" "$n_query"   # ← live bar on screen
              echo "[ $(date '+%F %T') ] Starting query $cont_ind / $n_query" >> "$outfile"
              wait_for_slot
              (
                  # temp_data/simu_data/truths...
                  Rscript script_generate_finite_basis_data_part4.R "$n_large" "$rep_i" "$adj_type" "$cont_ind" "$beta_truth" "${adj_params[@]}" >> "$outfile" 2>&1
                  echo "[ $(date '+%F %T') ] Finished query $cont_ind" >> "$outfile"
              ) & 
          done  
          wait
          
          # simu_data/adj_type_n.RData
          # simu_data/adj_type_n_truths.RData
          Rscript script_generate_finite_basis_data_part5.R "$n_large" "$rep_i" "$adj_type" "$n_query" "$groups" >> "$outfile" 2>&1
          
          step1_end=$(date +%s)
          step1_elapsed=$(( step1_end - step1_start ))
          echo "[DONE] Generation complete. Elapsed: ${step1_elapsed}s" | tee -a "$outfile"
          echo "" | tee -a "$outfile"
          echo "===================================================" >> "$outfile"
          echo "" | tee -a "$outfile"
          
          # ============================================================================
          # Step 2: Fit
          # ============================================================================
        
          
          echo "[STEP 2] Fitting dataset ..." | tee -a "$outfile"
          echo "" | tee -a "$outfile"
          echo "===================================================" >> "$outfile"
          echo "" | tee -a "$outfile"
          
          # ----------------
          # Part 1 - get all outer products before re-weighting by y_c and n
          # ----------------
          
          echo "[PART 1] Collecting parameters ..." >> "$outfile"
          
        
          # temp_data/simu/part1...
          output=$(Rscript script_fit_mice_data_part1.R \
                    "$model_type" "$n_large" "$n_large" "$rep_i" "$adj_type" "$method" "$X_truth" \
                    2>&1 | tee -a "$outfile")
        
          
          echo "" | tee -a "$outfile"
          echo "===================================================" >> "$outfile"
          echo "" | tee -a "$outfile"
          
          # Extract n_i and n_ij from output
        
          n_i=$(echo "$output" | grep "n_processes" | awk -F= '{print $2}')
          n_i=$(echo "$n_i" | xargs)
          
          n_ij=$(echo "$output" | grep "n_bivariate_processes" | awk -F= '{print $2}')
          n_ij=$(echo "$n_ij" | xargs)
          
        
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
              Rscript script_step2_part1_v5.R "$model_type" "$n_large" "$n_large" "$rep_i" "$adj_type" "$method" "$X_truth" "$k" >> "$outfile" 2>&1 &
          done
        
        
          echo "[PART 3] Calculating Rho_ij ..." >> "$outfile"
          echo "" | tee -a "$outfile"
        
          for kl in $(seq 1 "$n_ij"); do
              wait_for_slot
              
              # temp_data/simu/step_2_rho_ij...
              Rscript script_step2_part2_v5.R "$model_type" "$n_large" "$n_large" "$rep_i" "$adj_type" "$method" "$X_truth" "$kl" >> "$outfile" 2>&1 &
          done
          wait  
          
          echo "[PART 4] Merging all Rho_i and Rho_ij ..." >> "$outfile"
          echo "" | tee -a "$outfile"
          
          
          # temp_data/simu/step_2_v5_raw_rho_list...
          Rscript script_step2_part3_v5.R "$model_type" "$n_large" "$n_large" "$rep_i" "$adj_type" "$method" "$n_i" "$n_ij" >> "$outfile" 2>&1
          
        
          echo "[PART 5] Downstream estimation for all sub_n and y_c values ..." >> "$outfile"
          echo "" | tee -a "$outfile"
          
          echo "All n: ${ns[*]}" >> "$outfile"
          echo "Number of queries: $n_query" >> "$outfile"
          echo "" | tee -a "$outfile"
          
          # Loop over all subsets of n
          for n in "${ns[@]}"; do
              for j in $(seq 1 "$n_query"); do
                  echo "[START] n = $n, y_c = $j" >> "$outfile"
        
                  wait_for_slot
                  (
                      # gets weights and pads rho_i and rho_ii
                      # temp_data/simu/step_2_rho_list...
                      Rscript script_step2_part4_v5.R "$model_type" "$n_large" "$n" "$rep_i" "$adj_type" "$method" "$j" "$y_c_bandwidth" >> "$outfile" 2>&1
                      
            
                      # estimation until GIC
                      # temp_data/simu/part2b...
                      Rscript script_fit_mice_data_part2b_before_GIC.R "$model_type" "$n_large" "$n" "$rep_i" "$adj_type" "$method" "$X_truth" "$j" "$eigen_setting" "$y_c_bandwidth" >> "$outfile" 2>&1
                      
                      
                      # ---------------------------------------------------------
                      # GIC_LOCAL PROCEDURE
                      # ---------------------------------------------------------
                      
                      echo "At GIC local" >> "$outfile"
          
                      # PART 1: Precompute tau_c quantiles and get num_k
                      output=$(Rscript script_GIC_local_part1.R \
                               "$model_type" "$n_large" "$n" "$rep_i" "$adj_type" "$method" "$X_truth" "$j" \
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
                              echo "[GIC] Processing suffix ID: $id_suffix" >> "$outfile"
                          
                              for ((k=1; k<=num_k; k++)); do
                              
                                wait_for_slot
                                # Launch one R process per tau_c. 
                                # Inside this R script, it will loop through all tau_p (l=1...num_l)
                                Rscript script_GIC_local_part2and3_serial.R \
                                    "$model_type" "$n_large" "$n" "$rep_i" "$adj_type" "$method" \
                                    "$j" "$id_suffix" "$k" >> "$outfile" 2>&1 &                              
    
                                      
                              done
                              wait
                              echo "[GIC DONE] Processing suffix ID: $id_suffix" >> "$outfile"
                          ) &
                      done
                      wait
                      
                      echo "All GIC local tasks complete. Combining." >> "$outfile"
                      
                      # PART 4: Finalize and combine results
                      Rscript script_GIC_local_part4.R "$model_type" "$n_large" "$n" "$rep_i" "$adj_type" "$method" "$j" >> "$outfile" 2>&1
          
                      # ---------------------------------------------------------
                      
                      # estimation after GIC
                      # temp_data/simu/part3...
                      
                      echo "" | tee -a "$outfile"
                      echo "===================================================" >> "$outfile"
                      echo "" | tee -a "$outfile"
                      echo "Merging GIC local info with part2b" >> "$outfile"
                      
                      Rscript script_fit_mice_data_part2b_after_GIC.R "$model_type" "$n_large" "$n" "$rep_i" "$adj_type" "$method" "$X_truth" "$j" >> "$outfile" 2>&1
                      
                      
                      echo "[END] n = $n, y_c = $j" >> "$outfile"
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
                  echo "Running global/hybrid script..." >> "$outfile"
                  
      
                  # PART 1: Precompute tau_c quantiles and get num_k
                  output=$(Rscript script_GIC_global_part1.R \
                           "$model_type" "$n_large" "$n" "$rep_i" "$adj_type" "$method" "$n_query" \
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
                                  "$model_type" "$n_large" "$n" "$rep_i" "$adj_type" "$method" \
                                  "$id_suffix" "$k" >> "$outfile" 2>&1 &
                          fi
                          
                          # Run HYBRID if requested
                          wait_for_slot
                          
                          if [[ "$global_thresh_method" == "tau_c" || "$global_thresh_method" == "both" ]]; then
                              Rscript script_GIC_hybrid_part2and3_serial.R \
                                  "$model_type" "$n_large" "$n" "$rep_i" "$adj_type" "$method" \
                                  "$id_suffix" "$k" >> "$outfile" 2>&1 &
                          fi
    
                      done
    
                  done
                  wait
                  
                  echo "All GIC global/hybrid tasks complete. Combining." >> "$outfile"
                  
                  # PART 4: Finalize and combine results
                  
                  if [[ "$global_thresh_method" == "joint" || "$global_thresh_method" == "both" ]]; then
                      Rscript script_GIC_global_part4.R "$model_type" "$n_large" "$n" "$rep_i" "$adj_type" "$method" "$n_query" >> "$outfile" 2>&1
                  fi
                  
                  if [[ "$global_thresh_method" == "tau_c" || "$global_thresh_method" == "both" ]]; then
                      Rscript script_GIC_hybrid_part4.R "$model_type" "$n_large" "$n" "$rep_i" "$adj_type" "$method" "$n_query" >> "$outfile" 2>&1
                  fi
                  
                  
                  
    
              fi
              
              
              # ----------------
              # Part 2d - get step_12 and step_12b, ROC and edge set after all w_mats have been calculated                     
              # ----------------
              
              echo "At part 2d" >> "$outfile"
              
              wait_for_slot
              
              Rscript script_fit_mice_data_part2d.R "$model_type" "$n_large" "$n" "$rep_i" "$adj_type" "$method" "$n_query" >> "$outfile" 2>&1
              
              # ----------------
              # Part 3- when all part 2's are done, do part 3
              # ----------------
              
              wait_for_slot
              
              echo "At part 3" >> "$outfile"
              
              # simu_results/adj_type/CPGM/... .RData
              Rscript script_fit_mice_data_part3.R "$model_type" "$n_large" "$n" "$rep_i" "$adj_type" "$method" "$n_query" >> "$outfile" 2>&1
              
              echo "" | tee -a "$outfile"
              echo "===================================================" >> "$outfile"
              echo "" | tee -a "$outfile"
              echo "[DONE] Estimating n=$n in rep $rep_i out of $n_reps" >> "$outfile"
              echo "" | tee -a "$outfile"
              echo "===================================================" >> "$outfile"
              echo "" | tee -a "$outfile"
              
          done # end of each n
          
          echo "Deleting Files ..." >> "$outfile"
          
          
          BASE_DIR="./temp_data"
          
          # Construct the folder names
          FOLDER1="temp_data/simu_data_${adj_type}_n_${n_large}_rep_${rep_i}"
          FOLDER2="temp_data/simu_${adj_type}_n_${n_large}_rep_${rep_i}"
          
          # Execute the deletion
          # Using -v (verbose) so you can see exactly what is being deleted
          rm -rf "$FOLDER1"
          rm -rf "$FOLDER2"
          
          echo "Cleanup complete for $adj_type, n=$n_large, rep=$rep_i." >> "$outfile"
          
          # ============================================================
          # END TIMER
          # ============================================================
      
          end_time=$(date +%s)
          runtime=$((end_time - start_time))
          
          echo "Pipeline finished at: $(date)" >> "$outfile"
          echo "Total runtime: ${runtime} seconds (~$((runtime/60)) minutes)." >> "$outfile"
              
          
    
          
          
        done  # end of each rep
        
        wait
            
        # ----------------
        # Last step - visualize results
        # ----------------
        
        viz_log="script_outputs/simu/${adj_type}_visualization.log"
        rm -f "$viz_log"   # delete old log if it exists
    
        echo "===================================================" | tee -a "$viz_log"
        echo "" | tee -a "$viz_log"
        echo "[PART 3] Visualization ..." | tee -a "$viz_log"
        echo "" | tee -a "$viz_log"
        echo "===================================================" | tee -a "$viz_log"
        echo "" | tee -a "$viz_log"
    
        
        
        wait_for_slot
    
    
        Rscript script_unpack_finite_basis_results.R "$n_large" "${ns[*]}" "$n_reps" "$method" "$X_truth" "$beta_truth" "$eigen_setting" "$adj_type" "${adj_params[@]}" >> "$viz_log" 2>&1
    
        echo "" | tee -a "$viz_log"
        echo "===================================================" | tee -a "$viz_log"
        echo "" | tee -a "$viz_log"
        echo "[PART 3] DONE" | tee -a "$viz_log"
        
    ) &
    

done

wait


