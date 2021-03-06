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
#include <iostream>
#include <boost/random/uniform_int_distribution.hpp>

// CellType ///////////////////////////////////////////////////////////////////

// Statics:
unsigned int CellType::msNextId = 1;

// Constructors:
CellType::CellType(double birthrate, double alpha, double mu)
    : mId(msNextId++),
      mBirthRate(birthrate),
      mAlpha(alpha),
      mMu(mu),
      mNumMembers(0),
      mpMembers(0)
    {}


// Destructor:
CellType::~CellType() {
  while(mpMembers.size() != 0) {
    Cell* pCell = mpMembers.back();
    delete pCell;
  }
}


// Getter functions:
Cell* CellType::RandomMember() {

  //:ToDo: Replace with  proper sampling of index
  unsigned long max_i = mNumMembers - 1;
  boost::random::uniform_int_distribution< unsigned long> unif_member(0, max_i);
  unsigned long rm = unif_member(rng);

  // Debug messages:
  D(std::cout << std::endl;)
  D(std::cout << "########## Sampling ###########" << std::endl;)
  D(std::cout << "  Next reaction member:" << std::endl;)
  D(std::cout << "    Member No: " << mpMembers[rm]->Id() << std::endl;)
  D(std::cout << "    Member: " << mpMembers[rm] << std::endl;)
  D(std::cout << "    Type: " << this->Id() << std::endl;)
  D(std::cout << "###############################" << std::endl;)

  return mpMembers[rm];
}

// Setter functions:
void CellType::RegisterMember(Cell* pCell){
  pCell->TypeIndex(mpMembers.size());
  mpMembers.push_back(pCell); // Append pointer to new cell to member vector
  mNumMembers++;              // Increase the population count
}

void CellType::DeregisterMember(Cell* pCell){
  D(std::cout << "Deregister member " << pCell << std::endl;)
  std::vector<Cell*>::size_type i = pCell->TypeIndex();  // iterator for members vector.
  std::vector<Cell*>::size_type j = mpMembers.size() - 1;  // iterator for members vector.

  pCell->TypeIndex(-1);

  if (i != j) {
    mpMembers[j]->TypeIndex(i);
    mpMembers[i] = mpMembers[j];
  }
  mpMembers.pop_back();
  mNumMembers--; // Decrease the population count
}

// Print functions:
void CellType::Print(){
  std::cout << std::endl;
  std::cout << "########## Cell type ##########" << std::endl;
  std::cout << "  ID: " << mId << std::endl;
  std::cout << "  Birth rate: " << mBirthRate << std::endl;
  std::cout << "  Alpha: " << mAlpha << std::endl;
  std::cout << "  Members: " << mNumMembers << std::endl;
  std::cout << "###############################" << std::endl;
}

void CellType::PrintAllMembers(){
  std::cout << std::endl;
  std::cout << "#### Cell type members ########" << std::endl;
  std::cout << "  ID:" << mId << std::endl;
  std::cout << std::endl;
  for(std::vector<Cell*>::size_type i = 0; i != mpMembers.size(); i++) {
    std::cout << "  Member " << i << ": " << std::endl;
    std::cout << "    " << mpMembers[i]->Id() << std::endl;
  }
  std::cout << "###############################" << std::endl;
}
