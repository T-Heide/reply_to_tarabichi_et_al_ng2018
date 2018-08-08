1) Compilation of simulation program:
    - Requires the boost library, GNU-make, GCC. 
    - Compiled by executing 'make' in this directory. 
    - After the compilation execute the './gillespie_simulation -h'
      to see a help page explaining the program arguments.


2) Produce simulation results:
     - Bash and R scripts to reproduce the results are in './analysis_dir'.
     - The Bash script to run the simulations requires GNU-parallel.
     - The R scripts require the following R packages:
       neutralitytestr, dplyr, cowplot, reshape2.
