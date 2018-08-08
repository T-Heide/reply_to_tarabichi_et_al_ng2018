/*
    Simulation of a ABCG (Asymetric Branching Cell Growth)
    Copyright (C) 2018 Timon Heide (timon.heide@icr.ac.uk)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifdef _DEBUG_
#define D(x) x
#else
#define D(x)
#endif

#include <getopt.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <time.h>
#include "Phylogeny.hpp"
#include "Universe.hpp"
#include "Cell.hpp"
#include "CellType.hpp"
#include "extern_global_variables.hpp"


// Params of old function simulateTumorGrowth2D
#define WT_MU 16        // background mutation rate [int]
#define MU_MU 16        // mutant mutation rate [int]
#define WT_BR 1.0       // background birth rate [double]
#define MU_BR 1.0       // mutant birth rate   [double]
#define WT_DR 0.2       // background death rate [double]
#define MU_DR 0.2       // mutant death rate   [double]
#define CLST 256        // clone start time [int]
#define SIET 1048576    // simulation end time [int]
#define SEED 1          // random seed [time_t]
#define MINVAF 0.01
#define DEPTH 100

#define OUTPUT_PREFIX "./"
#define SEQUENCING_OUPUT_FILE "simulated_sequencing.tsv"
#define COUNT_OUPUT_FILE "cell_number.tsv"


boost::random::mt19937_64 rng;         // produces randomness out of thin air

int main(int argc, char* argv[]) {

  /****************************************************************************/
  /********************** Argument parsing ************************************/
  /****************************************************************************/

  // Cell parameters:
  double wildtype_mutation_rate = WT_MU;
  double mutant_mutation_rate = MU_MU;
  double wildtype_birthrate = WT_BR;
  double mutant_birthrate = MU_BR;
  double wildtype_deathrate = WT_DR;
  double mutant_deathrate = MU_DR;
  // Simulation parameters:
  int clone_start_time = CLST;
  int simulation_end_time = SIET;
  time_t seed = SEED;
  // Sequencing parameters:
  double minVaf = MINVAF;
  int depth = DEPTH;
  // Output options:
  std::string output_prefix = OUTPUT_PREFIX;

 // The possible arguments for getopts
  static struct option long_options[] = {
    // Cell parameters:
    {"wildtype_mutation_rate", optional_argument, 0, 'm'},
    {"mutant_mutation_rate", optional_argument, 0, 'M'},
    {"wildtype_birthrate",  optional_argument, 0, 'b'},
    {"mutant_birthrate",  optional_argument, 0, 'B'},
    {"wildtype_deathrate",  optional_argument, 0, 'd'},
    {"mutant_deathrate",  optional_argument, 0, 'D'},
    // Simulation parameters:
    {"clone_start_time",    optional_argument, 0, 't'},
    {"simulation_end_time",    optional_argument, 0, 'T'},
    {"seed",    optional_argument, 0,  's'},
    // Sequencing parameters:
    {"minVaf",    optional_argument, 0,  'f'},
    {"depth",    optional_argument, 0,  'x'},
    // Output options:
    {"output_prefix",    optional_argument, 0,  'o'},
    //
    {"help",    optional_argument, 0,  'h'},
    {NULL, 0, NULL, 0}
  };

  // Parse commandline arguments:
  while (true) {
    int c = getopt_long(argc, argv, "m:M:b:B:d:D:t:T:s:f:x:o:h",
                        long_options, NULL);

    if (c == -1)
      break;

    switch (c) {
      // Cell parameters:
      case 'm':
        wildtype_mutation_rate = atof(optarg);
          if (wildtype_mutation_rate < 0) {
            std::cout << "Error: " << std::endl;
            std::cout << "  Mutation rate (m) should be >= 0." << std::endl;
            exit(EXIT_FAILURE);
          }
        break;
      case 'M':
        mutant_mutation_rate = atof(optarg);
        if (mutant_mutation_rate < 0) {
          std::cout << "Error: " << std::endl;
          std::cout << "  Mutation rate (M) should be >= 0." << std::endl;
          exit(EXIT_FAILURE);
        }
        break;
      case 'b':
        wildtype_birthrate = atof(optarg);
        if (wildtype_birthrate < 0.0) {
          std::cout << "Error: " << std::endl;
          std::cout << "  Birth rate (b) should be >= 0." << std::endl;
          exit(EXIT_FAILURE);
        }
        break;
      case 'B':
        mutant_birthrate = atof(optarg);
        if (mutant_birthrate < 0.0) {
          std::cout << "Error: " << std::endl;
          std::cout << "  Birth rate (B) should be >= 0." << std::endl;
          exit(EXIT_FAILURE);
        }
        break;
      case 'd':
        wildtype_deathrate = atof(optarg);
        if (wildtype_deathrate < 0.0 || wildtype_deathrate > 1.0) {
          std::cout << "Error: " << std::endl;
          std::cout << "  Death rate (d) should be 0 <= d <= 1." << std::endl;
          exit(EXIT_FAILURE);
        }
        break;
      case 'D':
        mutant_deathrate = atof(optarg);
        if (mutant_deathrate < 0.0 || mutant_deathrate > 1.0) {
          std::cout << "Error: " << std::endl;
          std::cout << "  Death rate (D) should be 0 <= D <= 1." << std::endl;
          exit(EXIT_FAILURE);
        }
      break;
      // Simulation parameters:
      case 't':
        clone_start_time = atof(optarg);
        if (clone_start_time < 0) {
          std::cout << "Error: " << std::endl;
          std::cout << "  Clone start time (t) should be >= 0." << std::endl;
          exit(EXIT_FAILURE);
        }
        break;
      case 'T':
        simulation_end_time = atof(optarg);
        if (simulation_end_time < 0) {
          std::cout << "Error: " << std::endl;
          std::cout << "  Simulation end time (T) should be >= 0."
            << std::endl;
          exit(EXIT_FAILURE);
        }
        if (simulation_end_time <= clone_start_time) {
          std::cout << "Warning: " << std::endl;
          std::cout << "  Simulation end time (T) is >= clone start time (t)."
            << std::endl;
        }
        break;
      case 's':
        seed = atoi(optarg);
        break;
      // Sequencing parameters:
      case 'f':
        minVaf = atof(optarg);
        if (minVaf < 0.0 || minVaf > 1.0) {
          std::cout << "Error: " << std::endl;
          std::cout << "  VAF cutoff (f) should be 0 <= D <= 1." << std::endl;
          exit(EXIT_FAILURE);
        }
        break;
      case 'x':
        depth = atof(optarg);
        if (depth < 0) {
          std::cout << "Error: " << std::endl;
          std::cout << "  Sequencing depth (x) should be >= 0." << std::endl;
          exit(EXIT_FAILURE);
        }
      break;
      // Output options:
      case 'o':
        output_prefix = optarg;
        break;
      case ':':
        break;
      case 'h':
      case '?':
        std::cout << "\n";
        std::cout << "Simulation of a ABCG (Asymetric Branching Cell Growth)\n";
        std::cout << "Copyright (C) 2018, Timon Heide (timon.heide@icr.ac.uk)\n";
        std::cout << "\n";
        std::cout << "This program is free software; you can redistribute it and/or modify\n";
        std::cout << "it under the terms of the GNU General Public License as published by\n";
        std::cout << "the Free Software Foundation; either version 3 of the License, or\n";
        std::cout << "(at your option) any later version.";
        std::cout << "\n";
        std::cout << "This program is distributed in the hope that it will be useful,\n";
        std::cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
        std::cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n";
        std::cout << "GNU General Public License for more details.\n";
        std::cout << "\n";
        std::cout << "You should have received a copy of the GNU General Public License\n";
        std::cout << "along with this program. If not, see http://www.gnu.org/licenses/.\n";
        std::cout << "\n";
        std::cout << "\n";
        std::cout <<  "Usage: " << argv[0] << " [options]\n";
        std::cout <<  "  Options: \n";
        std::cout <<  "\n";
        // Cell parameters:
        std::cout <<  "    -m MU, --wildtype_mutation_rate=MU\n";
        std::cout <<  "        Mutation rate of wildtype cells (blue).\n";
        std::cout <<  "         Number of new mutations\n";
        std::cout <<  "        will drawn randomly from a Poisson \n";
        std::cout <<  "        distribution with lambda = MU during\n";
        std::cout <<  "        each division.\n";
        std::cout <<  "          default: " << WT_MU << std::endl;
        std::cout <<  "\n";
        std::cout <<  "    -M MU, --mutant_mutation_rate=MU\n";
        std::cout <<  "        Mutation ratemutated cells (red).\n";
        std::cout <<  "        See above option -m for details.\n";
        std::cout <<  "          default: " << MU_MU << std::endl;
        std::cout <<  "\n";
        std::cout <<  "    -b F, --wildtype_birthrate=F\n";
        std::cout <<  "        Birth rate of wildtype cells (blue).\n";
        std::cout <<  "          default: " << WT_BR << std::endl;
        std::cout <<  "\n";
        std::cout <<  "    -B F, --mutant_birthrate=F\n";
        std::cout <<  "        Birth rate of mutated cells (red).\n";
        std::cout <<  "          default: " << MU_BR << std::endl;
        std::cout <<  "\n";
        std::cout <<  "    -d F, --wildtype_deathrate=F\n";
        std::cout <<  "        Death rate of wildtype cells (blue).\n";
        std::cout <<  "          default: " << WT_DR << std::endl;
        std::cout <<  "\n";
        std::cout <<  "    -D F, --mutant_deathrate=F\n";
        std::cout <<  "        Death rate of mutated cells (red).\n";
        std::cout <<  "          default: " << MU_DR << std::endl;
        std::cout <<  "\n";
        // Simulation parameters:
        std::cout <<  "    -t F, --clone_start_time=F\n";
        std::cout <<  "        Time mutant clone is introduced.\n";
        std::cout <<  "          default: " << CLST << std::endl;
        std::cout <<  "\n";
        std::cout <<  "    -T F, --simulation_end_time=F\n";
        std::cout <<  "        Time simulations terminates.\n";
        std::cout <<  "          default: " << SIET << std::endl;
        std::cout <<  "\n";
        std::cout <<  "    -s SEED, --seed=SEED\n";
        std::cout <<  "        Random seed.\n";
        std::cout <<  "          default: " << SEED << std::endl;
        std::cout <<  "\n";
        // Sequencing parameters:
        std::cout <<  "    -f F, --min_vaf=F\n";
        std::cout <<  "        Sequencing detection limit.\n";
        std::cout <<  "          default: " << MINVAF << std::endl;
        std::cout <<  "\n";
        std::cout <<  "    -x F, --depth=F\n";
        std::cout <<  "        Sequencing coverage.\n";
        std::cout <<  "          default: " << DEPTH << std::endl;
        std::cout <<  "\n";
        // Output options:
        std::cout <<  "    -o DIR, --output_prefix=DIR\n";
        std::cout <<  "        Output prefix.\n";
        std::cout <<  "          default: " << OUTPUT_PREFIX << std::endl;
        std::cout <<  "\n";
        //
        std::cout <<  "    -h, --help\n";
        std::cout <<  "        Print this help\n";
        std::cout << std::endl;
        exit(EXIT_FAILURE);
	      break;
      default:
        std::cerr << " returned character code 0" << c << std::endl;
    }
  }

  // Tell user about the help page:
  if (argc <= 1) {
     std::cout << "\n";
     std::cout << "Simulation of a ABCG (Asymetric Branching Cell Growth)\n";
     std::cout << "Copyright (C) 2018, Timon Heide (timon.heide@icr.ac.uk)\n";
     std::cout << "\n";
     std::cout << "Programm will run with default parameters.\n";
     std::cout << "Please type '" << argv[0] << " -h' for help.\n"; 
     std::cout << std::endl;
  }

  // Print arguments:
  std::cout << std::endl;
  std::cout << "########## Options #############" << "\n";
  std::cout << "  Wildtype mutation rate: " << wildtype_mutation_rate << "\n";
  std::cout << "  Mutant mutation rate: " << mutant_mutation_rate << "\n";
  std::cout << "  Wildtype birth rate: " << wildtype_birthrate << "\n";
  std::cout << "  Mutant birth rate: " << mutant_birthrate << "\n";
  std::cout << "  Wildtype death rate: " << wildtype_deathrate << "\n";
  std::cout << "  Mutant death rate: " << mutant_deathrate << "\n";
  std::cout << "\n";
  std::cout << "  Clone start time: " << clone_start_time << "\n";
  std::cout << "  Simulation end time: " << simulation_end_time << "\n";
  std::cout << "  Random seed: " << seed << "\n";
  std::cout << "\n";
  std::cout << "  VAF cutoff: " << minVaf << "\n";
  std::cout << "  Sequencing depth: " << depth << "\n";
  std::cout << "\n";
  std::cout << "  Output prefix: " << output_prefix << "\n";
  std::cout << "################################" << std::endl;
  std::cout << std::endl;

  /****************************************************************************/
  /************************ Initialization ************************************/
  /****************************************************************************/

  // randomize seed
  rng.seed(seed);

  // Create the universe, a reusable pointer to handled cells:
  Universe* pUniverse = new Universe;
  Cell* pCell;

  // Other variables used for simulations:
  int generation = 0;
  double dt = 0.0;
  int action = -1;

  // Create two cell types (red and blue):
  CellType* pTypeBlue = new CellType(wildtype_birthrate,
                                     wildtype_deathrate,
                                     wildtype_mutation_rate);

  CellType* pTypeRed = new CellType(mutant_birthrate,
                                    mutant_deathrate,
                                    mutant_mutation_rate);

  // Create and insert a new blue cell:
  pCell = new Cell(pTypeBlue);
  pUniverse->InsertCell(pCell);


  // Open output files:
  std::string sequencing_output_file = output_prefix;
  sequencing_output_file += SEQUENCING_OUPUT_FILE;
  std::ofstream sequencing_output_stream(sequencing_output_file);
    if (!sequencing_output_stream.is_open()) {
      std::cout << "Error." << std::endl;
      std::cout << "  Unable to open sequencing output file:" << std::endl;
      std::cout << "    " << sequencing_output_file << std::endl;
      exit(EXIT_FAILURE);
    }

  std::string count_output_file = output_prefix;
  count_output_file += COUNT_OUPUT_FILE;
  std::ofstream count_output_stream(count_output_file);
    if (!count_output_stream.is_open()) {
      std::cout << "Error." << std::endl;
      std::cout << "  Unable to open cell count output file:" << std::endl;
      std::cout << "    " << count_output_file<< std::endl;
      exit(EXIT_FAILURE);
    }


  /****************************************************************************/
  /*************************** Simulation *************************************/
  /****************************************************************************/


  // Run till insertion of new cell type after clone_start_time should occure:
  while (++generation < clone_start_time) {
    pUniverse->NextReaction(&dt, &action)->RandomMember()->DoAction(&action);
    pUniverse->IncrementTimeBy(dt);
  } // stop running after reaching clone_start_time


  // Make conversion of a random blue cell:
  pTypeBlue->RandomMember()->Type(pTypeRed);


  // Run till simulation_end_time is reached
  while(++generation < simulation_end_time) {
    pUniverse->NextReaction(&dt, &action)->RandomMember()-> DoAction(&action);
    pUniverse->IncrementTimeBy(dt);
  } // stop running after reaching simulation_end_time


  /****************************************************************************/
  /***************************** Sampling *************************************/
  /****************************************************************************/


  // Get total cell number:
  int total_cells = pTypeBlue->NumMembers() + pTypeRed->NumMembers();

  // Print result:
  std::cout << std::endl;
  std::cout << "########## Result ##############" << std::endl;
  std::cout << "Total number of cells: " << total_cells << std::endl;
  std::cout << "Subclone I: " << pTypeBlue->NumMembers() << std::endl;
  std::cout << "Subclone II: " << pTypeRed->NumMembers() << std::endl;
  std::cout << "################################" << std::endl;
  std::cout << std::endl;

  // Put clone counts to ouput file and close the stream:
  count_output_stream << "total\tclone1\tclone2" << std::endl;
  count_output_stream << total_cells << "\t" << pTypeBlue->NumMembers() << "\t"
    << pTypeRed->NumMembers() << std::endl;
  count_output_stream.close();

  // Do sampling
  pUniverse->Sample(sequencing_output_stream, minVaf, depth);

  // Print result:
  std::cout << std::endl;
  std::cout << "Cell counts output: " << std::endl;
  std::cout << "  ->  " << count_output_file << std::endl;
  std::cout << std::endl;
  std::cout << "Sequencing result output: " << std::endl;
  std::cout << "  -> " << sequencing_output_file << std::endl;
  std::cout << std::endl;

  // Exit:
  return EXIT_SUCCESS;
}
