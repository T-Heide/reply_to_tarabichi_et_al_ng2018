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
#include "Cell.hpp"
#include "CellType.hpp"
#include "Phylogeny.hpp"
#include "Universe.hpp"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <iostream>
#include <array>
#include <vector>


// Statics:
unsigned long Cell::msNextId = 1;


// Constructors:
Cell::Cell () // Default definition without given location and types.
  : mId(msNextId++),
    mpUniverse(0),
    mpType(0),
    mpNode(0)
  {}

Cell::Cell (CellType* pType) // Mostly default but a given cell type.
  : mId(msNextId++),
    mpUniverse(0),
    mpType(pType),
    mpNode(0)
  {
    mpType->RegisterMember(this);
  }


// Destructor:
Cell::~Cell(){ // Proper deletion of a cell.
  D(std::cout << "Cell " << mId << "(" << this << ") dies!" << std::cout;) // Debug message
  mpType->DeregisterMember(this);      // Let the associated type forget.
  mpNode->AssociatedCell(0);           // Unlink from the associated node.
  // Removal of cell complete.
}

int Cell::TypeIndex() {return mTypeIndex;}


// Getter functions:
unsigned long Cell::Id() {return(mId);}
CellType* Cell::Type() {return mpType;}

bool Cell::AsProgenitorDies() { // Random realization to die with prob alpha.
  if (mpType->NumMembers() == 1) {
    return false;
  } else {
    double alpha = mpType->Alpha();
    boost::random::bernoulli_distribution<> bern_alpha(alpha);
    bool res = bern_alpha(rng);
    return res;
  }
}

PhylogenyNode* Cell::AssociatedNode() {return mpNode;} // Assoc. node
Universe* Cell::AssociatedUniverse() {return mpUniverse;} // Assoc. nodeå

// Setter functions:
void Cell::Type(CellType* pNewType) {
  mpUniverse->RegisterType(pNewType);
  mpType->DeregisterMember(this);
  pNewType->RegisterMember(this);
  mpType = pNewType;

  if (mpNode != 0) { // If the cell has a associated node, member of a universe,
    mpNode->TypeId(pNewType->Id());
  }

  // Update the TypeId of all ancestral nodes to 0 (undertermined/mix type):
  PhylogenyNode *pCurrentNode = mpNode;
  while ((pCurrentNode = pCurrentNode->UpNode()) != 0) {
    pCurrentNode->TypeId(0);
  }
}

void Cell::AssociatedNode(PhylogenyNode* new_node) { mpNode = new_node; }
void Cell::AssociatedUniverse(Universe* new_universe) {
  mpUniverse = new_universe;
}

void Cell::TypeIndex(int newIndex) {mTypeIndex = newIndex;}


// Other functions:
void Cell::MutateCell() {
  Cell::MutateCell(mpType->Mu());
}

void Cell::MutateCell(double mu){
  // Sample number of new muts from poisson:
  boost::random::poisson_distribution<int> dist_mutations(mu);
  int new_muts = dist_mutations(rng);

  // Append muts to node:
  mpNode->AddNewMutations(new_muts);
}


void Cell::Print() {
  std::cout << std::endl;
  std::cout << "############ Cell ##############" << std::endl;
  std::cout << "   ID: " << mId << std::endl;
  std::cout << "   Location:" << std::endl;
  std::cout << "       Universe: " << mpUniverse << std::endl;
  std::cout << "   Type: " << mpType->Id() << std::endl;
  std::cout << "       Birth rate: " << mpType->Birthrate() << std::endl;
  std::cout << "###############################" << std::endl;
}

void Cell::DoAction(int *action) {
  switch(*action) {
    case 1: // Action 'Divide'
      this->Divide();
      break;
    default:
      std::cerr << "Unknow action '" << action << "'" << std::endl;
      exit(EXIT_FAILURE);
      break;
  }
}

void Cell::Divide() {

  // Debug messages:
  D(std::cout << std::endl;)
  D(std::cout << "########### Division ###########" << std::endl;)
  D(std::cout << "  Cell:" << this->Id() << std::endl;)

  // Cells that were not introduced into a universe can't divide!
  if (mpUniverse == 0) {
    std::cerr << "Cells not introduced into universe can't divide!" << std::endl;
    return;
  }

  if(this->AsProgenitorDies()){
    D(std::cout << "  Cell died during division" << std::endl;)
    D(std::cout << "###############################" << std::endl;)
    delete this;
  } else {

    // Store current associated node, then branch to the left and mutate:
    PhylogenyNode* old_node = this->mpNode;
    PhylogenyNode* new_left_node = new PhylogenyNode(this, old_node);
    old_node->LeftNode(new_left_node);
    this->MutateCell();

    // Create daughter, insert on branch to the right of old node and mutate:
    Cell* pDaughter = new Cell(mpType);
    PhylogenyNode* new_right_node = new PhylogenyNode(pDaughter, old_node);
    old_node->RightNode(new_right_node);
    mpUniverse->InsertCell(pDaughter, false);
    pDaughter->MutateCell();

    D(std::cout << "  Cell divided." << std::endl;)
    D(std::cout << "###############################" << std::endl;)
  }
}
