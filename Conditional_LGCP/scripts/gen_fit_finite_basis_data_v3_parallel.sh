#!/bin/bash

# ---------------------------------------------------------------------------
# Generate AND fit finite basis data — Global Semaphore Parallelization
# ---------------------------------------------------------------------------

cd "$(dirname "$0")/.."  # go one level up

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
max_jobs=60  # GLOBAL LIMIT
min_events=10
max_events=20000
n_query=8
beta_0=5.4
beta_truth="T"
X_truth="T"
same_basis="T"
constant_d=2

# ---------------------------------------------------------------------------
# GLOBAL SEMAPHORE SETUP
# ---------------------------------------------------------------------------
# This creates a "bucket" of tokens. Every process must take one to run.
# This works across subshells because they all read from the same file descriptor.

res_fifo="/tmp/fifo.$$"
mkfifo "$res_fifo"
exec 100<>"$res_fifo"
rm -f "$res_fifo"

# Fill the semaphore with 'tokens' (newlines)
for ((i=0; i<max_jobs; i++)); do
    echo >&100
done

# Function to run a command using a semaphore token
# Usage: run_globally <command...>
run_globally() {
    read -u 100 # Claim a token (blocks if none available)
    (
        "$@"
        echo >&100 # Release token back to pool
    ) &
}

# ---------------------------------------------------------------------------
# PREPARE OUTPUTS
# ---------------------------------------------------------------------------
mkdir -p script_outputs/simu

# ---------------------------------------------------------------------------
# MAIN LOOP: Parallelized across adj_types
# ---------------------------------------------------------------------------

for entry in "${adj_type_params[@]}"; do
  
  # Start the adjacency type logic in the background immediately
  # It will internally respect the global max_jobs
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
      step1_start=$(date +%s)

      # Part 0 is small/setup, run directly
      Rscript script_generate_finite_basis_data_part0.R "$n_large" "$n_query" "$beta_0" "$adj_type" "${adj_params[@]}"
      
      # Parallel Group Generation
      for group_idx in $(seq 1 "$groups"); do
          run_globally Rscript script_generate_finite_basis_data_parts_1_and_2.R "$n_large" "$adj_type" "$group_idx" "$n_group" "$min_events" "$max_events"
      done  
      wait # Wait for all group generation for THIS adj_type to finish
      
      Rscript script_generate_finite_basis_data_part3.R "$n_large" "$adj_type" "$groups"
      
      # Parallel Truth Calculation
      for cont_ind in $(seq 1 "$n_query"); do
          run_globally Rscript script_generate_finite_basis_data_part4.R "$n_large" "$adj_type" "$cont_ind" "$beta_truth" "${adj_params[@]}"
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
              
              # Parallel rho_i estimation
              for k in $(seq 1 "$n_i"); do
                  run_globally Rscript script_step2_part1.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$X_truth" "$j" "$k"
              done
              
              # Parallel rho_ij estimation
              for kl in $(seq 1 "$n_ij"); do
                  run_globally Rscript script_step2_part2.R "$model_type" "$n_large" "$n" "$adj_type" "$method" "$X_truth" "$j" "$kl"
              done
              wait # Wait for all rho estimates for this query/adj_type
              
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
exec 100>&- # Close the semaphore file descriptor
echo "All datasets processed successfully."