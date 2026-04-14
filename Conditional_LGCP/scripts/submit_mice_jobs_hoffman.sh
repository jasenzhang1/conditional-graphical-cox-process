# for ID in "Tau1" "Tau2" "Tau3" "WT1" "WT2" "WT3"; do
#   for mov in "0" "1"; do
#     for vr in "0" "1"; do
#       qsub -N "${ID}_m${mov}_vr${vr}" \
#            -v ID="$ID",MOV="$mov",VR="$vr",BW="0.001",REGION="HIP" \
#            fit_mice_data_v6.sh
#     done
#   done
# done

for ID in "Tau1" "Tau2" "Tau3" "WT1" "WT2" "WT3"; do
  for mov in "0" "1"; do
    for vr in "0" "1"; do
      qsub -N "${ID}_m${mov}vr${vr}" \
           fit_mice_data_v6.sh "$ID" "$mov" "$vr" "0.001" "HIP"
    done
  done
done

# to run:
# qsub submit_mice_jobs_hoffman.sh