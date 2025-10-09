#!/bin/bash

# Loop over dataset indices

cd "$(dirname "$0")/.."   # go one level up (from /scripts to /)

adj_type_params=(
  "banded_c2 0 1 0.3"
  "banded_trig2 0 1 0.3"
  "sparse_v2 0 1 2 0.3 -1 2 0.01"
)
n=2000

mkdir -p script_outputs

for entry in "${adj_type_params[@]}"; do
  adj_type=$(echo "$entry" | awk '{print $1}')
  outfile="script_outputs/generate_data_${adj_type}_n_${n}.txt"
  echo "Generating dataset: n=$n, adj_type=$adj_type"
  echo "Logging to: $outfile"
  Rscript script_generate_data.R $n $entry > "$outfile" 2>&1 &
done

wait

echo "All datasets finished."

