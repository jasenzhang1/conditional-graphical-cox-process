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
  "block_banded_v2 0 1 0 0.5"
)

n=300
methods=("CPGM")
n_large=300
max_jobs=30

mkdir -p script_outputs

# ---------------------------------------------------------------------------
# Generate data + Fit model (single unified log file per adj_type)
# ---------------------------------------------------------------------------

for entry in "${adj_type_params[@]}"; do
  adj_type=$(echo "$entry" | awk '{print $1}')
  outfile="script_outputs/${adj_type}_n_${n}.log"

  echo "===================================================" | tee -a "$outfile"
  echo "Starting full pipeline for adj_type=$adj_type, n=$n" | tee -a "$outfile"
  echo "Logging to: $outfile" | tee -a "$outfile"
  echo "===================================================" | tee -a "$outfile"
  echo "" | tee -a "$outfile"

  # -------------------
  # Step 1: Generate
  # -------------------
  echo "[STEP 1] Generating dataset..." | tee -a "$outfile"
  Rscript script_generatefinite_basis_data.R $n $entry >> "$outfile" 2>&1
  echo "[DONE] Generation complete." | tee -a "$outfile"
  echo "" | tee -a "$outfile"

  # -------------------
  # Step 2: Fit
  # -------------------
  for method in "${methods[@]}"; do
    echo "[STEP 2] Fitting method=$method, n_large=$n_large, adj_type=$adj_type" | tee -a "$outfile"

    # Wait until fewer than max_jobs are running
    while (( $(jobs -r | wc -l) >= max_jobs )); do
      sleep 1
    done

    Rscript script_fit_generated_finite_basis_data.R "$n_large" "$n" "$adj_type" "$method" >> "$outfile" 2>&1 &
  done
done

wait

echo ""
echo "===================================================" 
echo "All datasets finished generating and fitting."
echo "===================================================" 
