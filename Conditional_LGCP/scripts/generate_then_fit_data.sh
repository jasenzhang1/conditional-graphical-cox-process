#!/bin/bash

# Run first script
./generate_data.sh

# Run second script after first finishes
./fit_generated_data.sh
