#!/bin/bash

# Loop over dataset indices

cd "$(dirname "$0")/.."   # go one level up (from /scripts to /)

IDs=(346 351 366 361 362 368)
time_scales=(1 2 5 10)
max_jobs=12

mkdir -p script_outputs

outfile="script_outputs/generate_mice_data.txt"
echo "Logging to: $outfile"
# start fresh each run
> "$outfile"

for ID in "${IDs[@]}"; do
  for time_scale in "${time_scales[@]}"; do
  
    # Wait until fewer than max_jobs are running
    while (( $(jobs -r | wc -l) >= max_jobs )); do
      sleep 1
    done  
  
    echo "Generating dataset: mouse=$ID, timescale=$time_scale"
    
    # Run Rscript in background and measure time
    (
      start_time=$(date +%s)
      Rscript script_mice_data.R $ID $time_scale >> "$outfile" 2>&1
      end_time=$(date +%s)
      elapsed=$((end_time - start_time))
      minutes=$((elapsed / 60))
      seconds=$((elapsed % 60))
      echo "Mouse $ID, timescale $time_scale finished in ${minutes}m ${seconds}s" >> "$outfile"
    ) &  
    
  done
done

wait

echo "Finished generating all mice datasets."

