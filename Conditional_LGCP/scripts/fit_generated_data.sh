#!/bin/bash

# Loop over dataset indices

cd "$(dirname "$0")/.." || exit 1   # go one level up (from /scripts to /) and exit if fails

methods=("OG" "CPGM")

adj_types=(
  #"banded_c2"
  #"banded_trig2"
  #"sparse_v2"
  #"indep_2"
  #"single_c2"
  #"single_v2"
  "single_j2"
)

n_large=2000

ns=(500 1000 1500 2000)
max_jobs=12  

mkdir -p script_outputs

for n in "${ns[@]}"; do

  for adj_type in "${adj_types[@]}"; do
    mkdir -p script_outputs/${adj_type}
    
    for method in "${methods[@]}"; do
      mkdir -p script_outputs/${adj_type}/${method}
      
      outfile="script_outputs/${adj_type}/${method}/fit_log_n${n}.txt"
      echo "Fitting dataset: n_large=$n_large, n=$n, adj_type=$adj_type, method=$method"
      echo "Logging to: $outfile"
      
      # Wait until fewer than max_jobs are running
      while (( $(jobs -r | wc -l) >= max_jobs )); do
        sleep 1
      done    
      
      Rscript script_fit_generated_data.R "$n_large" "$n" "$adj_type" "$method" > "$outfile" 2>&1 &
    done
  done
done

wait

echo "All datasets finished fitting."

