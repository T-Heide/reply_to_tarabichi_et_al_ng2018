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

#include "extern_global_variables.hpp"
#include "CellType.hpp"
#include "Cell.hpp"
#include "Universe.hpp"
#include "Phylogeny.hpp"
#include <stdlib.h>
#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <fstream>
#include <cmath>
#include <boost/random/uniform_real.hpp>


// Universe ////////////////////////////////////////////////////////////////////


// Constructor:
Universe::Universe()
  : mTime(0.0) {}


// Destructor:
Universe::~Universe(){

  // Remove all CellTypes:
  for(std::vector<CellType*>::size_type i = 0; i < mpTypes.size(); i++) {
    delete mpTypes[i++];
  }

  // Remove the complete phylogeny:
  for(std::vector<PhylogenyRoot*>::size_type i = 0; i < mpPhylogenies.size();) {
    delete mpPhylogenies[i++];
  }
}

// Getter functions:
double Universe::Time() { return mTime; }

CellType* Universe::NextReaction(double* r_delta_time, int* action){
  // Samples and returns the next reaction that occures.

  // Exit if the universe contains no types:
  if (mpTypes.size() == 0) {
    std::cerr << "In Universe::NextReaction: No type to choose." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Debug messages:
  D(std::cout << std::endl;)
  D(std::cout << "########## Sampling ###########" << std::endl;)
  D(std::cout << "  Next reaction type:" << std::endl;)
  D(std::cout << std::endl;)


  // Keep track of minum index and dt:
  double delta_time_minimum;
  int index_minimum;

  boost::uniform_real<double> runi_dist(0.0, 1.0);

  // Calculate dt for all types in the universe:
  for(std::vector<CellType*>::size_type i = 0; i < mpTypes.size(); i++) {
    // Debug messages:
    D(std::cout << "    Candidate " << i <<  ":" << std::endl;)
    D(std::cout << "      No: " << i << std::endl;)
    D(std::cout << "      Type: " << mpTypes[i]->Id() << std::endl;)
    D(std::cout << "      BR: " << mpTypes[i]->Birthrate() << std::endl;)
    D(std::cout << "      N: " << mpTypes[i]->NumMembers() << std::endl;)


    double runi = runi_dist(rng);
    double c_birthrate = mpTypes[i]->Birthrate();
    double c_number = mpTypes[i]->NumMembers();
    double delta_time_current = (log(1) - log(runi)) / (c_number * c_birthrate);

    // Debug messages:
    D(std::cout << "      X: " << runi << std::endl;)
    D(std::cout << "      dt: " <<  delta_time_current <<  std::endl;)
    D(std::cout << std::endl;)

    // Keep track of minum index and dt:
    if (i == 0 || delta_time_minimum > delta_time_current) {
      index_minimum = i;
      delta_time_minimum = delta_time_current;
    }
  }

  // Debug messages:
  D(std::cout << "  Selected:" << std::endl;)
  D(std::cout << "    Type No: #" << index_minimum << std::endl;)
  D(std::cout << "    Type: " << mpTypes[index_minimum]->Id() << std::endl;)
  D(std::cout << "    dt: " << delta_time_minimum << std::endl;)
  D(std::cout << "###############################" << std::endl;)

  // Return results:
  *r_delta_time = delta_time_minimum;
  *action = 1; // action divide

  return(mpTypes[index_minimum]);
}



bool Universe::Sample(std::vector <int> &clone , std::vector <double>& vaf,
                      double minVAF, int depth) {

  // Determine total number of cells:
  int ncells = 0;
  for (std::vector<CellType*>::size_type i = 0; i < mpTypes.size(); i++){
    ncells += mpTypes[i]->NumMembers();
  }

  // Traverse the Phylogenies:
  int ncells2 = 0;
  for (std::vector<CellType*>::size_type i = 0; i < mpPhylogenies.size(); i++){
    ncells2 += mpPhylogenies[i]->Root()->SampleNode(clone, vaf, minVAF,
                                                    depth, ncells);
  }

  if (ncells != ncells2) {
    std::cerr << "sum_cells != total_cells: " << ncells << " vs " <<
      ncells2 << std::endl;
  }

  return true;
}

bool Universe::Sample(std::ofstream &ostream, double minVAF, int depth) {
  ostream << "VAF\tALT\tDP\tCLONE" << std::endl;

  // Determine total number of cells:
  int ncells = 0;
  for (std::vector<CellType*>::size_type i = 0; i < mpTypes.size(); i++){
    ncells += mpTypes[i]->NumMembers();
  }

  // Traverse the Phylogenies:
  int ncells2 = 0;
  for (std::vector<CellType*>::size_type i = 0; i < mpPhylogenies.size(); i++){
    ncells2 += mpPhylogenies[i]->Root()->SampleNode(ostream, minVAF,
                                                    depth, ncells);
  }

  if (ncells != ncells2) {
    std::cerr << "sum_cells != total_cells: " << ncells << " vs " <<
      ncells2 << std::endl;
  }

  return true;
}

// Setter functions:
void Universe::IncrementTimeBy(double Delta){ mTime += Delta; }


bool Universe::InsertCell(Cell* pCell) {
  return InsertCell(pCell, true);
}

bool Universe::InsertCell(Cell* pCell, bool is_new_lineage) {

  // Debug messages:
  D(std::cout << std::endl;)
  D(std::cout << "########## Insertion ##########" << std::endl;)
  D(std::cout << "   ID: " << pCell->Id() << std::endl;)
  D(std::cout << "   Location:" << std::endl;)
  D(std::cout << "       Universe: " << this << std::endl;)
  D(std::cout << "###############################" << std::endl;)

  // Insert type into universe
  if (is_new_lineage) {
    this->RegisterType(pCell->Type());
  }

  // Set universe variable in cell:
  pCell->AssociatedUniverse(this);

  // Update universe:
  number_of_cells++;

  // Register new phylogeny:
  if (pCell->AssociatedNode() == 0) {
    PhylogenyRoot* pNewPhylo = new PhylogenyRoot(pCell);
    mpPhylogenies.push_back(pNewPhylo);
  }

  return true;
}

void Universe::RegisterType(CellType* p_new_type){
  std::vector<CellType *>::size_type i = 0;

  while(i < mpTypes.size()) {
    D(std::cout << "Comparing new cell type " <<  p_new_type;)
    D(std::cout << " with " << mpTypes[i] << std::endl;)
    if (mpTypes[i] == p_new_type)
      break;
    i++;
  }

  if (i == mpTypes.size()) { // reached end.
    D(std::cout << "Registering new cell type " <<  p_new_type << std::endl;)
    mpTypes.push_back(p_new_type);
  }
}
