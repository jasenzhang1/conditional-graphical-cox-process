#!/bin/bash

Rscript Alzheimers_run_05_t1.R > cox_output_05_t1.txt 2>&1 &
Rscript Alzheimers_run_05_t2.R > cox_output_05_t2.txt 2>&1 &
Rscript Alzheimers_run_05_t5.R > cox_output_05_t5.txt 2>&1 &
Rscript Alzheimers_run_05_t10.R > cox_output_05_t10.txt 2>&1 &

wait  # Wait for all background jobs to finish

