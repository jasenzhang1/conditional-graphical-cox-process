#!/bin/bash

# Loop over dataset indices

cd "$(dirname "$0")/.." || exit 1   # go one level up (from /scripts to /) and exit if fails

methods=("OG")

adj_types=(
  "single_c2"
  "single_v2"
)


max_jobs=12  
i=1
j=2

mkdir -p "script_outputs"



for adj_type in "${adj_types[@]}"; do
  mkdir -p "script_outputs/${adj_type}"
  
  for method in "${methods[@]}"; do
    mkdir -p "script_outputs/${adj_type}/${method}"
    
    outfile="script_outputs/${adj_type}/${method}/unpack_log.txt"
    echo "Unpacking dataset: adj_type=$adj_type, method=$method"
    echo "Logging to: $outfile"
    
    # Wait until fewer than max_jobs are running
    while (( $(jobs -r | wc -l) >= max_jobs )); do
      sleep 1
    done    
    
    Rscript script_unpack_results.R "simu_results" "$adj_type" "$method" "$i" "$j" > "$outfile" 2>&1 &
  done
done


wait

echo "All datasets finished unpacking."

