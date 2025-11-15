#!/bin/bash

# ---------------------------------------------------------------------------
# Generate AND fit finite basis data — single unified log per dataset
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
  "block_banded_v2 0 1 0.4 0.8 2"
)

n_large=1024
ns=(128 256 512 1024)
method="CPGM"
max_jobs=30

mkdir -p script_outputs

# ---------------------------------------------------------------------------
# Generate data + Fit model (single unified log file per adj_type)
# ---------------------------------------------------------------------------

for entry in "${adj_type_params[@]}"; do
  adj_type=$(echo "$entry" | awk '{print $1}')
  outfile="script_outputs/${adj_type}_n_${n_large}.log"
  rm -f "$outfile"   # delete old log if it exists

  echo "===================================================" | tee -a "$outfile"
  echo "Starting full pipeline for adj_type=$adj_type, n=$n_large" | tee -a "$outfile"
  echo "Logging to: $outfile" | tee -a "$outfile"
  echo "Start time: $(date)" | tee -a "$outfile"
  echo "===================================================" | tee -a "$outfile"
  echo "" | tee -a "$outfile"

  total_start=$(date +%s)

  # -------------------
  # Step 1: Generate
  # -------------------
  echo "[STEP 1] Generating dataset..." | tee -a "$outfile"
  step1_start=$(date +%s)

  Rscript script_generate_finite_basis_data.R $n_large $entry >> "$outfile" 2>&1

  step1_end=$(date +%s)
  step1_elapsed=$(( step1_end - step1_start ))
  echo "[DONE] Generation complete. Elapsed: ${step1_elapsed}s" | tee -a "$outfile"
  echo "" | tee -a "$outfile"

  # -------------------
  # Step 2: Fit for multiple n
  # -------------------
  
  echo "[STEP 2] Start" | tee -a "$outfile"
  step2_very_beginning=$(date +%s)
  
  for n in "${ns[@]}"; do
      # Wait if max_jobs reached
      while (( $(jobs -r | wc -l) >= max_jobs )); do
          sleep 1
      done
  
      # Coarse fit
      (
          step_start=$(date +%s)
          {
              echo "[STEP 2-coarse] Fitting for n=$n, method=$method"
              Rscript script_fit_generated_finite_basis_data.R \
                  "$n_large" "$n" "$adj_type" "$method"
              step_end=$(date +%s)
              echo "[DONE] Coarse fit for n=$n complete. Elapsed: $((step_end - step_start))s"
          } >> "$outfile" 2>&1
      ) &
  
      # Fine fit
      (
          step_start=$(date +%s)
          {
              echo "[STEP 2-fine] Fitting for n=$n, method=$method"
              Rscript script_fit_generated_finite_basis_data_finer.R \
                  "$n_large" "$n" "$adj_type" "$method"
              step_end=$(date +%s)
              echo "[DONE] Fine fit for n=$n complete. Elapsed: $((step_end - step_start))s"
          } >> "$outfile" 2>&1
      ) &
  done
  
  wait
  
  step2_very_end=$(date +%s)
  echo "[DONE] STEP 2 End" | tee -a "$outfile"
  step2_total_elapsed=$(( step2_very_end - step2_very_beginning ))
  echo "[DONE] Fitting for all n. Elapsed: ${step2_total_elapsed}s" | tee -a "$outfile"



  total_end=$(date +%s)
  total_elapsed=$(( total_end - total_start ))
  echo "" | tee -a "$outfile"
  echo "---------------------------------------------------" | tee -a "$outfile"
  echo "Pipeline completed for $adj_type." | tee -a "$outfile"
  echo "Total time: ${total_elapsed}s" | tee -a "$outfile"
  echo "---------------------------------------------------" | tee -a "$outfile"
  echo "" | tee -a "$outfile"
done

echo ""
echo "===================================================" 
echo "All datasets finished generating and fitting."
echo "===================================================" 
