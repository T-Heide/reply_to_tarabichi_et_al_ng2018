#!/bin/bash

# Get arguments:
sc_mutation_rate=$1
sc_birth_rate=$2
seed=$3
start_time=$4

# Program path:
progPath="../gillespie_simulation"

# Output options:
output_path="results/simulations/mmr_$sc_mutation_rate/mbr_$sc_birth_rate/"
out_file_prefix="$output_path/simulation-mmr_$sc_mutation_rate"
out_file_prefix+="-mbr_$sc_birth_rate-seed_$seed-cst_$start_time-"

# Create output dir: 
mkdir -p "$output_path"

# Run simulation:
"$progPath" -B "$sc_birth_rate" -M "$sc_mutation_rate" -s "$seed" \
            -t "$start_time" -o "$out_file_prefix"

