#!/bin/bash

# ---------------------------------------------------------------------------
# Generate AND fit finite basis data — single unified log per dataset
# 
# 11/19/2025 - v2: parallelize bivariate estimation
#            - v3: parallelize data generation
# 
#            - v5: calculate rho_i and rho_ij before normalizing them by different n and y_c
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
  
  #"hub_block_v2 0 1 4 4 3 0.1 1.5 0.1 1"
  #"hub_block_c2 0 1 4 2 2 0.7 0.7"
  #"hub_block_c0 0.5 4 2 2 0.7 0.7"
  #"hub_block_j2 0 1 4 0.5 2 2 0.7 0.7"
  
  
  "complete_block_v2 0 1 4 4 3 0.1 1.1 0.1 0.8"         # for complete, c3 < c1/3, c4 < c2/3 [0.1, 1.1] and [0.1, 0.8]
  #"complete_block_c0 0.5 4 2 2 0.7 0.7"
  #"complete_block_c2 0 1 4 2 2 0.7 0.7"
  #"complete_block_j2 0 1 4 0.5 2 2 0.7 0.7"
  
  "flexible_block_banded_v2 0 1 4 4 3 0.1 1.6 0.1 1.2"   # for flexible, a bit more chill [0.1, 1.6] and [0.1, 1.2]
  #"flexible_block_banded_c2 0 1 4 2 2 0.7 0.7"
  #"flexible_block_banded_c0 0.5 4 2 2 0.7 0.7"  
  #"flexible_block_banded_j2 0 1 4 0.5 2 2 0.7 0.7"
)


n_large=20000
ns=(100 250 500 1000 2500 5000 10000 20000)

n_group=10
groups=$(( n_large / n_group ))

method="CPGM"
model_type="simu"  # simu or mice
max_jobs=60
min_events=10
max_events=5000

p=12
d=2
n_query=4
beta_0=4.7
beta_truth="T"
X_truth="T"
eigen_setting="trig_and_joint" #only_joint, trig_and_joint
global_thresh_method="both" #both, joint, tau_c, neither   both = do joint and tau_c



function wait_for_slot {
    # 1. Use pgrep to find processes named exactly "R" owned by the current user
    # 2. This creates a global throttle across all subshells
    while true; do
        # Count R and Rscript processes
        running=$(pgrep -u "$USER" -x "R|Rscript" | wc -l)
        
        if (( running < max_jobs )); then
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

    read -r -a fields <<< "$entry"
  
    adj_type="${fields[0]}"         # first field
    adj_params=("${fields[@]:1}")   # all fields after the first

    wait_for_slot
    (

        outfile="script_outputs/simu/${adj_type}_n_${n_large}.log"
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
        Rscript script_generate_finite_basis_data_part0.R "$p" "$d" "$n_large" "$n_query" "$beta_0" "$adj_type" "${adj_params[@]}" >> "$outfile" 2>&1
        
        for group_idx in $(seq 1 "$groups"); do
            print_bar "$group_idx" "$groups"   # ← live bar on screen
            echo "[ $(date '+%F %T') ] Starting group $group_idx / $groups" >> "$outfile"
            wait_for_slot
            (
                # temp_data/simu_data/parts1_and_2...
                Rscript script_generate_finite_basis_data_parts_1_and_2.R "$n_large" "$adj_type" "$group_idx" "$n_group" "$min_events" "$max_events" >> "$outfile" 2>&1
                
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
        
        Rscript script_generate_finite_basis_data_part3.R "$n_large" "$adj_type" "$groups" >> "$outfile" 2>&1
        
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
                Rscript script_generate_finite_basis_data_part4.R "$n_large" "$adj_type" "$cont_ind" "$beta_truth" "${adj_params[@]}" >> "$outfile" 2>&1
                echo "[ $(date '+%F %T') ] Finished query $cont_ind" >> "$outfile"
            ) & 
        done  
        wait
        
        # simu_data/adj_type_n.RData
        # simu_data/adj_type_n_truths.RData
        Rscript script_generate_finite_basis_data_part5.R "$n_large" "$adj_type" "$n_query" "$groups" >> "$outfile" 2>&1
        
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
                  "$model_type" "$n_large" "$n_large" "$adj_type" "$method" "$X_truth" \
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
            Rscript script_step2_part1_v5.R "$model_type" "$n_large" "$n_large" "$adj_type" "$method" "$X_truth" "$k" >> "$outfile" 2>&1 &
        done
      
      
        echo "[PART 3] Calculating Rho_ij ..." >> "$outfile"
        echo "" | tee -a "$outfile"
      
        for kl in $(seq 1 "$n_ij"); do
            wait_for_slot
            
            # temp_data/simu/step_2_rho_ij...
            Rscript script_step2_part2_v5.R "$model_type" "$n_large" "$n_large" "$adj_type" "$method" "$X_truth" "$kl" >> "$outfile" 2>&1 &
        done
        wait  
        
        echo "[PART 4] Merging all Rho_i and Rho_ij ..." >> "$outfile"
        echo "" | tee -a "$outfile"
        
        
        # temp_data/simu/step_2_v5_raw_rho_list...
        Rscript script_step2_part3_v5.R "$model_type" "$n_large" "$n_large" "$adj_type" "$method" "$n_i" "$n_ij" >> "$outfile" 2>&1
        
      
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
                    Rscript script_step2_part4_v5.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$j" >> "$outfile" 2>&1
                    
          
                    # estimation from step3 onwards
                    # temp_data/simu/part3...
                    Rscript script_fit_mice_data_part2b.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$X_truth" "$j" "$eigen_setting" >> "$outfile" 2>&1
                    
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
            
            Rscript script_fit_mice_data_part2c.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$n_query" "$global_thresh_method" >> "$outfile" 2>&1
            
            # ----------------
            # Part 3- when all part 2's are done, do part 3
            # ----------------
            
            wait_for_slot
            
            # simu_results/adj_type/CPGM/... .RData
            Rscript script_fit_mice_data_part3.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$n_query" >> "$outfile" 2>&1
            
            echo "" | tee -a "$outfile"
            echo "[DONE] Estimating n=$n" >> "$outfile"
            
            
        done
        
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
    
        # temp_data/simu/step_2_v5_rho_i...
        Rscript script_unpack_finite_basis_results.R "$n_large" "${ns[*]}" "$method" "$X_truth" "$beta_truth" "$eigen_setting" "$adj_type" "${adj_params[@]}" >> "$outfile" 2>&1
        
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
    

    
done

wait

