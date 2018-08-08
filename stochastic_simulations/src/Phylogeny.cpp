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
#include "Phylogeny.hpp"
#include "Cell.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>

// PhylogenyRoot ///////////////////////////////////////////////////////////////

// Statics:
unsigned int PhylogenyRoot::msNextId = 0;

// Constructor:
PhylogenyRoot::PhylogenyRoot(Cell* pCell) : mId(msNextId++) {
  mpRoot = new PhylogenyNode(pCell);
}

PhylogenyNode* PhylogenyRoot::Root(){ return mpRoot; }

//Destructor
PhylogenyRoot::~PhylogenyRoot() {
  delete mpRoot;
}

// PhylogenyNode ///////////////////////////////////////////////////////////////

// Statics:
unsigned long PhylogenyNode::msNextId = 0;

// Constructors
PhylogenyNode::PhylogenyNode(Cell* pCell) //
  : mId(msNextId++),
    mpCell(pCell),
    mpUp(0),
    mpLeft(0),
    mpRight(0),
    mGeneration(1),
    mNumMutsGeneration(0)
  {
    if (pCell->AssociatedNode() != 0) {
      pCell->AssociatedNode()->AssociatedCell(0);
    }
    pCell->AssociatedNode(this);

    mTypeId = pCell->Type() == 0 ? 0 : pCell->Type()->Id();
  }

PhylogenyNode::PhylogenyNode(Cell* pCell, PhylogenyNode* pUp)
  : mId(msNextId++),
  mpCell(pCell),
  mpUp(pUp),
  mpLeft(0),
  mpRight(0),
  mNumMutsGeneration(0)
  {
    mGeneration = pUp->Generation() + 1;
    if (pCell->AssociatedNode() != 0) {
      pCell->AssociatedNode()->AssociatedCell(0);
    }
    pCell->AssociatedNode(this);

    mTypeId = pCell->Type() == 0 ? 0 : pCell->Type()->Id();
  }

//Destructor
PhylogenyNode::~PhylogenyNode() {

  // Delete left node:
  if (mpLeft != 0) {
    delete mpLeft;
    mpLeft = 0;
  }

  // Delete right node:
  if (mpRight != 0) {
    delete mpRight;
    mpRight = 0;
  }

  // Delete associated cell:
  if (mpCell != 0) {
    delete mpCell;
    mpCell = 0;
  }

  // Clean pointers in the up node:
  if (mpUp != 0) {
    if (mpUp->LeftNode() == this) {
      mpUp->LeftNode(0);
    } else if (mpUp->RightNode() == this) {
      mpUp->RightNode(0);
    }
  }
}

int PhylogenyNode::TypeId() {return mTypeId;}


// Getters:
unsigned long PhylogenyNode::Id() {return mId;}
unsigned int PhylogenyNode::Generation(){return mGeneration;}
unsigned int PhylogenyNode::NumMutations(){return mNumMutsGeneration;}
Cell* PhylogenyNode::AssociatedCell() {return mpCell;}
PhylogenyNode* PhylogenyNode::UpNode() {return mpUp;}
PhylogenyNode* PhylogenyNode::LeftNode() {return mpLeft;}
PhylogenyNode* PhylogenyNode::RightNode() {return mpRight;}

std::vector <PhylogenyNode*> PhylogenyNode::NodeAncestry(){
  std::vector <PhylogenyNode*> result;
  PhylogenyNode* pCurrent = this;
  do {
    result.push_back(pCurrent);
  } while ((pCurrent = pCurrent->UpNode()) != 0);
  return result;
}


void simulateNGS(double vaf, int dp, double& sim_vaf, int& sim_alt, int& sim_dp)
{
  boost::random::poisson_distribution<int> dist_depth(dp);
  sim_dp = dist_depth(rng);


  boost::random::binomial_distribution<int> binom_alt(sim_dp, vaf);
  sim_alt = binom_alt(rng);
  sim_vaf = sim_alt * 1.0 / sim_dp;
  return;
}

int PhylogenyNode::SampleNode(std::vector <int>&clone,
                              std::vector <double>& vaf,
                              double minVAF, int depth, int total_cells)
{
  int n_cells = 0;

  if (mpLeft != 0) {
    n_cells += mpLeft->SampleNode(clone, vaf, minVAF, depth, total_cells);
  }

  if (mpRight != 0) {
    n_cells += mpRight->SampleNode(clone, vaf, minVAF, depth, total_cells);
  }

  if (mpCell != 0) {
    n_cells += 1;
  }

  if (n_cells != 0) {
    double exp_vaf = 0.5 * n_cells / total_cells;
    double sim_vaf = 0.0;
    int sim_alt = 0;
    int sim_dp = 0;

    for (unsigned int i = 0; i < mNumMutsGeneration; i++) {
      simulateNGS(exp_vaf, depth, sim_vaf, sim_alt, sim_dp);
      if (sim_vaf >= minVAF) {
        vaf.push_back(sim_vaf);
        clone.push_back(this->TypeId());
      }
    }
  }

  return n_cells;
}

int PhylogenyNode::SampleNode(std::ofstream& ostream,
                              double minVAF, int depth, int total_cells)
{
  int n_cells = 0;

  if (mpLeft != 0) {
    n_cells += mpLeft->SampleNode(ostream, minVAF, depth, total_cells);
  }

  if (mpRight != 0) {
    n_cells += mpRight->SampleNode(ostream, minVAF, depth, total_cells);
  }

  if (mpCell != 0) {
    n_cells += 1;
  }

  if (n_cells != 0) {
    double exp_vaf = 0.5 * n_cells / total_cells;
    double sim_vaf = 0.0;
    int sim_alt = 0.0;
    int sim_dp = 0;

    for (unsigned int i = 0; i < mNumMutsGeneration; i++) {
      simulateNGS(exp_vaf, depth, sim_vaf, sim_alt, sim_dp);
      if (sim_vaf >= minVAF) {
        std::stringstream ss;
        ss << std::scientific << std::setprecision(3) << sim_vaf << "\t";
        ss << sim_alt << "\t" << sim_dp << "\t" << this->TypeId();
        ostream << ss.str()  << std::endl;
      }
    }
  }

  return n_cells;
}


// Setters:
void PhylogenyNode::AddNewMutations(int num_new_mutations) {
  mNumMutsGeneration += num_new_mutations;
}

void PhylogenyNode::AssociatedCell(Cell* pCell) {mpCell = pCell;}

void PhylogenyNode::UpNode(PhylogenyNode* pNode) {mpUp = pNode;}
void PhylogenyNode::LeftNode(PhylogenyNode* pNode) {mpLeft = pNode;}
void PhylogenyNode::RightNode(PhylogenyNode* pNode) {mpRight = pNode;}
void PhylogenyNode::TypeId(int newTypeId) {mTypeId = newTypeId;}



// Output:
void PhylogenyNode::Print() {
  unsigned int up = (mpUp == 0) ? 0 : mpUp->Id();
  unsigned int left = (mpLeft == 0) ? 0 : mpLeft->Id();
  unsigned int right = (mpRight == 0) ? 0 : mpRight->Id();

  std::cout << "###### Phylogenetic Node ######" << std::endl;
  std::cout << "   ID: " << mId << std::endl;
  std::cout << "   Up: " << up << std::endl;
  std::cout << "   Left: " << left << std::endl;
  std::cout << "   Right: " << right << std::endl;
  std::cout << "   Generation: " << mGeneration << std::endl;
  std::cout << "   Mutations: " << mNumMutsGeneration << std::endl;
  if (mpCell != 0) { std::cout << "   Cell: " << mpCell << std::endl; }
  std::cout << "###############################" << std::endl;
}

void PhylogenyNode::PrintAncestry() {

  if (mpLeft == 0 && mpRight == 0) { // If this is a leaf
    std::string cell_id = (mpCell == 0) ? "NA" : std::to_string(mpCell->Id());
    std::cout << std::endl;
    std::cout << "###### Ancestry of cell #######" << std::endl;
    std::cout << "  Cell ID: " << cell_id << std::endl;
    std::cout << std::endl;
  }

  this->Print();

  if (mpUp != 0) { // If next is not root
    std::cout << "              | " << std::endl;
    std::cout << "              | " << std::endl;
    mpUp->PrintAncestry();
  } else {
    std::cout << std::endl;
    std::cout << "###############################" << std::endl;
  }
}
