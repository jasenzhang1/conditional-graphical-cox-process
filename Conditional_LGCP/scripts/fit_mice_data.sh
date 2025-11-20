#!/bin/bash

# ============================================================
# START TIMER
# ============================================================
sh_outfile="script_outputs/mice/fit_mice.txt"
rm -f "$sh_outfile"

start_time=$(date +%s)
echo "Pipeline started at: $(date)" >> "$sh_outfile"

cd "$(dirname "$0")/.." || exit 1   # go one level up (from /scripts to /) and exit if fails

method="CPGM"

max_jobs=30

# function wait_for_slot {
#     while (( $(jobs -rp | wc -l) >= max_jobs )); do
#         sleep 1
#     done
# }

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

ID="Tau3"
y_c_structure="week_only"
time_scale=10
m=20
movement=(0 0 1 1)
VR=(0 1 0 1)
max_processes=500

mkdir -p script_outputs
mkdir -p script_outputs/mice

# -------------------------------------------------
# preprocess
# -------------------------------------------------


# Loop over all combinations of movement and VR
for i in "${!movement[@]}"; do
  mov=${movement[i]}
  vr=${VR[i]}
  
  outfile="script_outputs/mice/preprocess_m${mov}vr${vr}.txt"
  echo "Preprocessing for movement=$mov, VR=$vr" >> "$sh_outfile"
  echo "Logging to: $outfile" >> "$sh_outfile"
  
  # Wait if max_jobs are running
  while (( $(jobs -r | wc -l) >= max_jobs )); do
    sleep 1
  done
  
  # Launch in background
  Rscript script_preprocess_mice_data.R "$ID" "$y_c_structure" "$time_scale" "$method" "$m" "$mov" "$vr" "$max_processes" > "$outfile" 2>&1 &
done

# Wait for all background jobs to finish
wait



echo "All preprocessing jobs finished." >> "$sh_outfile"

end_time=$(date +%s)
runtime=$((end_time - start_time))

echo "Total elapsed time: ${runtime} seconds (~$((runtime/60)) minutes)." >> "$sh_outfile"
echo "===========================================" >> "$sh_outfile"

# ----------
# fit
# ----------

# Loop over all combinations of movement and VR
for i in "${!movement[@]}"; do

    mov=${movement[i]}
    vr=${VR[i]}
    
    outfile="script_outputs/mice/fit_m${mov}vr${vr}.txt"
    echo "Part 1 of Strata $i Starting" >> "$sh_outfile"
    

    
    
    # ----------------
    # Part 1
    # ----------------
    
    wait_for_slot
    output=$(Rscript script_fit_mice_data_part1.R \
              "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" \
              2>&1 | tee "$outfile")
    wait
    
    echo "=========================================" >> "$outfile"
    
    # Extract n_queries from output
    n_queries=$(echo "$output" | grep "n_queries" | awk -F= '{print $2}')
    n_queries=$(echo "$n_queries" | xargs)

    n_i=$(echo "$output" | grep "n_i" | awk -F= '{print $2}')
    n_i=$(echo "$n_i" | xargs)
    
    n_ij=$(echo "$output" | grep "n_ij" | awk -F= '{print $2}')
    n_ij=$(echo "$n_ij" | xargs)
    
    # ----------------
    # Part 2 - parallelize each y_c_query
    # ----------------
    
    echo "Part 2 of Strata $i Starting" >> "$sh_outfile"
    for j in $(seq 1 "$n_queries"); do

        echo "Query $j out of $n_queries" >> "$sh_outfile"
        

        
        # ----------------
        # Part 2a - within each fitting procedure, do rho_i and rho_ij estimation all together, and then collect
        # ----------------       
        
        wait_for_slot
        Rscript script_step2_part0.R "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" >> "$outfile" 2>&1
        
        for k in $(seq 1 "$n_i"); do
            wait_for_slot
            Rscript script_step2_part1.R "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" "$k" >> "$outfile" 2>&1 &
        done
        
        echo "Query $j out of $n_queries done with rho_i" >> "$sh_outfile"
        
        for kl in $(seq 1 "$n_ij"); do
            wait_for_slot
            Rscript script_step2_part2.R "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" "$kl" >> "$outfile" 2>&1 &
        done
        
        echo "Query $j out of $n_queries done with rho_ij" >> "$sh_outfile"
        
        wait
        
        wait_for_slot
        Rscript script_step2_part3.R "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_i" "$n_ij" >> "$outfile" 2>&1 &
        
        # ----------------
        # Part 2b - now continue for the rest of the estimation
        # ---------------- 
        
        wait_for_slot
        Rscript script_fit_mice_data_part2b.R "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" >> "$outfile" 2>&1 &
        

    done
    
    wait
    
    echo "=========================================" >> "$outfile"
    
    # ----------------
    # Part 3- when all part 2's are done, do part 3
    # ----------------
    
    echo "Part 3 of Strata $i Starting" >> "$sh_outfile"
    
    wait_for_slot
    Rscript script_fit_mice_data_part3.R "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_queries" >> "$outfile" 2>&1 &

done

wait

# ============================================================
# END TIMER
# ============================================================

echo "===========================================" >> "$sh_outfile"

end_time=$(date +%s)
runtime=$((end_time - start_time))

echo "Pipeline finished at: $(date)" >> "$sh_outfile"
echo "Total runtime: ${runtime} seconds (~$((runtime/60)) minutes)." >> "$sh_outfile"
