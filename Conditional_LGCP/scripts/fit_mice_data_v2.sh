#!/bin/bash

# ------------------------------------------------------------------------------
#
# GOAL: Fit data replicates for mice data
# 
# v2:   (1/12/2026)
#       Incorporate parallelized estimation code 
#
# 
# inputs:
#
# 
# - IDs          (vector)   mice ID that we want to make data for
# - time_scales  (vector)   how long is each data replicate? 2, 5, or 10 seconds
# - m            (integer)  granularity of time_grid
# - max_processes  (integer)  do we truncate the number of neurons?
# ------------------------------------------------------------------------------

cd "$(dirname "$0")/.."   # go one level up (from /scripts to /)


IDs=("Tau3")
y_c_structure="week_only"
method="CPGM"
model_type="mice"  # simu or mice
max_jobs=30

time_scale=10 
m=20
movement=(0 0 1 1)
VR=(0 1 0 1)
max_processes=500

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

mkdir -p script_outputs
mkdir -p script_outputs/mice

# -------------------------------------------------
# preprocess
# -------------------------------------------------


# Loop over all combinations of movement and VR
for ID in "${IDs[@]}"; do

  outfile="script_outputs/mice/${ID}_t${time_scale}.log"   #Tau1_t10.log
  rm -f "$outfile"   # delete old log if it exists
  
  
  echo "===================================================" | tee -a "$outfile"
  echo "Starting full pipeline for mouse=$ID, time_scale=$time_scale" | tee -a "$outfile"
  echo "Logging to: $outfile" | tee -a "$outfile"
  echo "Start time: $(date)" | tee -a "$outfile"
  echo "===================================================" | tee -a "$outfile"
  echo "" | tee -a "$outfile"
  
  start_time=$(date +%s)

  # -------------------
  # Step 1: Preprocess
  # -------------------
  echo "[STEP 1] Preprocess dataset..." | tee -a "$outfile"

  for i in "${!movement[@]}"; do
    mov=${movement[i]}
    vr=${VR[i]}
    
    # Wait if max_jobs are running
    while (( $(jobs -r | wc -l) >= max_jobs )); do
      sleep 1
    done
    
    # Launch in background
    Rscript script_preprocess_mice_data.R "$ID" "$y_c_structure" "$time_scale" "$method" "$m" "$mov" "$vr" "$max_processes" > "$outfile" 2>&1 &
  done
  
  # Wait for all background jobs to finish
  wait

  
  end_time=$(date +%s)
  runtime=$((end_time - start_time))
  
  echo "[STEP 1] All preprocessing jobs finished." | tee -a "$outfile"
  echo "Total elapsed time: ${runtime} seconds (~$((runtime/60)) minutes)." >> "$outfile"
  echo "===========================================" >> "$outfile"
  echo "" | tee -a "$outfile"
  
  # -----------------------
  # Step 2: Fit
  # -----------------------
  
  echo "[STEP 2] Fitting dataset..." | tee -a "$outfile"
  echo "" | tee -a "$outfile"
  
  # Loop over all combinations of movement and VR
  for i in "${!movement[@]}"; do
  
      mov=${movement[i]}
      vr=${VR[i]}
      

      echo "Starting movement $mov, vr $vr" >> "$outfile"
      
      
      # ----------------
      # Part 1
      # ----------------
      
      wait_for_slot
      output=$(Rscript script_fit_mice_data_part1.R \
                "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" \
                2>&1 | tee "$outfile")
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
      
      echo "Part 2 of Strata $i (mov=$mov, vr=$vr) Starting" >> "$outfile"
      for j in $(seq 1 "$n_queries"); do
  
          echo "Query $j out of $n_queries" >> "$outfile"
          
  
          
          # ----------------
          # Part 2a - within each fitting procedure, do rho_i and rho_ij estimation all together, and then collect
          # ----------------       
          
          wait_for_slot
          Rscript script_step2_part0.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" >> "$outfile" 2>&1
          
          for k in $(seq 1 "$n_i"); do
              wait_for_slot
              Rscript script_step2_part1.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" "$k" >> "$outfile" 2>&1 &
          done
          
          echo "Query $j out of $n_queries done with rho_i" >> "$outfile"
          
          for kl in $(seq 1 "$n_ij"); do
              wait_for_slot
              Rscript script_step2_part2.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" "$kl" >> "$outfile" 2>&1 &
          done
          
          echo "Query $j out of $n_queries done with rho_ij" >> "$outfile"
          
          wait
          
          wait_for_slot
          Rscript script_step2_part3.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" "$n_i" "$n_ij" >> "$outfile" 2>&1 &
          
          # ----------------
          # Part 2b - now continue for the rest of the estimation
          # ---------------- 
          
          wait_for_slot
          Rscript script_fit_mice_data_part2b.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$j" >> "$outfile" 2>&1 &
          
  
      done
      
      wait
      
      echo "=========================================" >> "$outfile"
      
      # ----------------
      # Part 3- when all part 2's are done, do part 3
      # ----------------
      
      echo "Part 3 of Strata $i Starting" >> "$outfile"
      
      wait_for_slot
      Rscript script_fit_mice_data_part3.R "$model_type" "$ID" "$y_c_structure" "$time_scale" "$method" "$mov" "$vr" "$n_queries" >> "$outfile" 2>&1 &
  
  done  # discrete strata loop
done  # mouse loop

wait

# ============================================================
# END TIMER
# ============================================================

echo "===========================================" >> "$outfile"

end_time=$(date +%s)
runtime=$((end_time - start_time))

echo "Pipeline finished at: $(date)" >> "$outfile"
echo "Total runtime: ${runtime} seconds (~$((runtime/60)) minutes)." >> "$outfile"
