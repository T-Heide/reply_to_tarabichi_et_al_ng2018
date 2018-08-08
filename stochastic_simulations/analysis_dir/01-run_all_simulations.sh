#!/bin/bash

# Number of simulations to run in parallel:
NUMPROC=4

# Run simulation for all possible combinations of defined parameters:
parallel --jobs "$NUMPROC" -N1 ./run_simulation.sh \
  :::: selected_sc_mutation_rates.txt \
  :::: selected_sc_birth_rates.txt \
  ::: {1000..1199} \
  ::: 256

