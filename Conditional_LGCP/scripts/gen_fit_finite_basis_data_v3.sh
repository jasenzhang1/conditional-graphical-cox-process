#!/bin/bash

# ---------------------------------------------------------------------------
# Generate AND fit finite basis data — single unified log per dataset
# 
# 11/19/2025 - v2: parallelize bivariate estimation
#            - v3: parallelize data generation
# ---------------------------------------------------------------------------





cd "$(dirname "$0")/.."  # go one level up (from /scripts to /)

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
  "flexible_block_banded_c0 0.5 6 2 2.9 0.9"
  "flexible_block_banded_c2 0 1 6 2 2.9 0.9"
  "flexible_block_banded_v2 0 1 6 2 1.5 2.9 0.5 0.9"
)

# comment 1

n_large=1500
ns=(500 1000 1500)
n_group=10
groups=$(( n_large / n_group ))
method="CPGM"
model_type="simu"  # simu or mice
max_jobs=60
min_events=10
max_events=20000
n_query=2
beta_0=4.8
beta_truth="T"
X_truth="T"

function wait_for_slot {
    # Wait until the number of background jobs is strictly less than max_jobs
    while true; do
        running=$(jobs -rp | wc -l)
        if (( running < max_jobs )); then
            break
        fi
        sleep 0.5
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
  Rscript script_generate_finite_basis_data_part0.R "$n_large" "$n_query" "$beta_0" "$adj_type" "${adj_params[@]}" >> "$outfile" 2>&1
  
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
  
  # temp_data/simu_data/dataset...
  Rscript script_generate_finite_basis_data_part3.R "$n_large" "$adj_type" "$groups" >> "$outfile" 2>&1
  echo "[STEP 1] Finished merging events" | tee -a "$outfile"
  
  # ---------------------------------------------------
  # Step 1b: Get truths for each cont_ind query
  # ---------------------------------------------------
  
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
  
  # ----------------
  # Step 2: Fit
  # ----------------
  
  echo "[STEP 2] Fitting dataset..." | tee -a "$outfile"
  echo "" | tee -a "$outfile"
  echo "Estimating the following sample sizes: ${ns[@]}"
  echo "" | tee -a "$outfile"
  
  # Loop over all subsets of n
  for n in "${ns[@]}"; do
  

      
      echo "Part 1 of Estimating n=$n Starting" >> "$outfile"
      
  
      
      # ----------------
      # Part 1
      # ----------------

  
      wait_for_slot
      
      # temp_data/simu/part1...
      output=$(Rscript script_fit_mice_data_part1.R \
                "$model_type" "$n_large" "$n" "$adj_type" "$method" "$X_truth" \
                2>&1 | tee -a "$outfile")
      wait
      
      echo "" | tee -a "$outfile"
      echo "===================================================" >> "$outfile"
      
      
      # Extract n_queries from output
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
      
      # ----------------
      # Part 2 - parallelize each y_c_query
      # ----------------
      
      echo "Part 2 of Estimating n=$n Starting" >> "$outfile"
      for j in $(seq 1 "$n_queries"); do
  
          echo "Query $j out of $n_queries" >> "$outfile"
          
  
          
          # ----------------
          # Part 2a - within each fitting procedure, do rho_i and rho_ij estimation all together, and then collect
          # ----------------       
          
          wait_for_slot
          
          # temp_data/simu/part2...
          Rscript script_step2_part0.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$j" >> "$outfile" 2>&1
          
          for k in $(seq 1 "$n_i"); do
              wait_for_slot
              
              # temp_data/simu/step_2_rho_i...
              Rscript script_step2_part1.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$X_truth" "$j" "$k" >> "$outfile" 2>&1 &
          done
          
          echo "Query $j out of $n_queries [1/4] done with rho_i" >> "$outfile"
          
          for kl in $(seq 1 "$n_ij"); do
              wait_for_slot
              
              # temp_data/simu/step_2_rho_ij...
              Rscript script_step2_part2.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$X_truth" "$j" "$kl" >> "$outfile" 2>&1 &
          done
          wait
          echo "Query $j out of $n_queries [2/4] done with rho_ij" >> "$outfile"

          wait_for_slot
          
          # temp_data/simu/step_2_rho_list...
          Rscript script_step2_part3.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$X_truth" "$j" "$n_i" "$n_ij" >> "$outfile" 2>&1
          echo "Query $j out of $n_queries [3/4] done with merger" >> "$outfile"
          
          # ----------------
          # Part 2b - now continue for the rest of the estimation
          # ---------------- 
          
          wait_for_slot
          
          # temp_data/simu/part3...
          Rscript script_fit_mice_data_part2b.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$X_truth" "$j" >> "$outfile" 2>&1
          echo "Query $j out of $n_queries [4/4] finished" >> "$outfile"
      done
      
      wait
      
      echo "" | tee -a "$outfile"
      echo "===================================================" >> "$outfile"
      
      
      # ----------------
      # Part 3- when all part 2's are done, do part 3
      # ----------------
      
      echo "Part 3 of Estimating n=$n Starting" >> "$outfile"
      
      wait_for_slot
      
      # simu_results/adj_type/CPGM/... .RData
      Rscript script_fit_mice_data_part3.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$n_queries" >> "$outfile" 2>&1
      
      echo "[DONE] Estimating n=$n" >> "$outfile"
  
  done
  
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
done


