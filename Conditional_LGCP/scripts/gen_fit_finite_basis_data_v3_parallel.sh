#!/bin/bash

# ---------------------------------------------------------------------------
# Generate AND fit finite basis data — Parallelized Adj_Types
# ---------------------------------------------------------------------------

cd "$(dirname "$0")/.."  

adj_type_params=(
  "hub_block_v2 0 1 4 2 2 0.5 0.9 0.5 0.9"
  "hub_block_c2 0 1 4 2 2 0.7 0.7"
  "hub_block_c0 0.5 4 2 2 0.7 0.7"
  "hub_block_j2 0 1 4 0.5 2 2 0.7 0.7"
  "complete_block_c0 0.5 4 2 2 0.7 0.7"
  "complete_block_c2 0 1 4 2 2 0.7 0.7"
  "complete_block_v2 0 1 4 2 2 0.5 0.9 0.5 0.9"
  "complete_block_j2 0 1 4 0.5 2 2 0.7 0.7"
  "flexible_block_banded_v2 0 1 2 2 0.5 0.9 0.5 0.9"
  "flexible_block_banded_c2 0 1 2 2 0.7 0.7"
  "flexible_block_banded_c0 0.5 2 2 0.7 0.7"  
  "flexible_block_banded_j2 0 1 0.5 2 2 0.7 0.7"
)

n_large=1000
ns=(1000)
n_group=10
groups=$(( n_large / n_group ))
method="CPGM"
model_type="simu"
max_jobs=60  # Total global background processes allowed
min_events=10
max_events=20000
n_query=8
beta_0=5.4
beta_truth="T"
X_truth="T"
same_basis="T"
constant_d=2

function wait_for_slot {
    while true; do
        running=$(jobs -rp | wc -l)
        if (( running < max_jobs )); then
            break
        fi
        sleep 1
    done
}

mkdir -p script_outputs/simu

# ---------------------------------------------------------------------------
# MAIN LOOP: Parallelized across adj_types
# ---------------------------------------------------------------------------

for entry in "${adj_type_params[@]}"; do
  
  # Wait for a slot before starting a NEW adjacency type pipeline
  wait_for_slot

  (
    read -r -a fields <<< "$entry"
    adj_type="${fields[0]}"
    adj_params=("${fields[@]:1}")
    
    outfile="script_outputs/simu/${adj_type}_n_${n_large}.log"
    rm -f "$outfile"

    echo "[$adj_type] Starting pipeline..."
    start_time=$(date +%s)

    {
      echo "==================================================="
      echo "Starting full pipeline for adj_type=$adj_type, n=$n_large"
      echo "Start time: $(date)"
      echo "==================================================="

      # --- STEP 1: GENERATE ---
      echo "[STEP 1] Generating dataset..."
      Rscript script_generate_finite_basis_data_part0.R "$n_large" "$n_query" "$beta_0" "$adj_type" "${adj_params[@]}"
      
      for group_idx in $(seq 1 "$groups"); do
          wait_for_slot
          (
              Rscript script_generate_finite_basis_data_parts_1_and_2.R "$n_large" "$adj_type" "$group_idx" "$n_group" "$min_events" "$max_events"
          ) &
      done  
      wait
      
      Rscript script_generate_finite_basis_data_part3.R "$n_large" "$adj_type" "$groups"
      
      # --- STEP 1b: TRUTHS ---
      for cont_ind in $(seq 1 "$n_query"); do
          wait_for_slot
          (
              Rscript script_generate_finite_basis_data_part4.R "$n_large" "$adj_type" "$cont_ind" "$beta_truth" "${adj_params[@]}"
          ) & 
      done  
      wait
      
      Rscript script_generate_finite_basis_data_part5.R "$n_large" "$adj_type" "$n_query" "$groups"

      # --- STEP 2: FIT ---
      for n in "${ns[@]}"; do
          output=$(Rscript script_fit_mice_data_part1.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$X_truth" 2>&1)
          
          n_queries=$(echo "$output" | grep "n_queries" | awk -F= '{print $2}' | xargs)
          n_i=$(echo "$output" | grep "n_processes" | awk -F= '{print $2}' | xargs)
          n_ij=$(echo "$output" | grep "n_bivariate_processes" | awk -F= '{print $2}' | xargs)

          for j in $(seq 1 "$n_queries"); do
              Rscript script_step2_part0.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$j"
              
              for k in $(seq 1 "$n_i"); do
                  wait_for_slot
                  Rscript script_step2_part1.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$X_truth" "$j" "$k" &
              done
              
              for kl in $(seq 1 "$n_ij"); do
                  wait_for_slot
                  Rscript script_step2_part2.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$X_truth" "$j" "$kl" &
              done
              wait
              
              Rscript script_step2_part3.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$X_truth" "$j" "$n_i" "$n_ij"
              Rscript script_fit_mice_data_part2b.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$X_truth" "$j" "$same_basis" "$constant_d"
          done
          
          Rscript script_fit_mice_data_part3.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$n_queries"
      done

      end_time=$(date +%s)
      runtime=$((end_time - start_time))
      echo "Total runtime for $adj_type: ${runtime} seconds."
    } >> "$outfile" 2>&1

    echo "[$adj_type] Finished."
  ) &

done

wait
echo "All datasets processed."