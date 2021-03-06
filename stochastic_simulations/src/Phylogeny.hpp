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


#ifndef PHYLOGENY_H
#define PHYLOGENY_H


// Forward declerations: ///////////////////////////////////////////////////////
class PhylogenyNode;
class Cell;
class Universe;

// Includes: ///////////////////////////////////////////////////////////////////
#include <vector>
#include <fstream>

// Phylo ///////////////////////////////////////////////////////////////////////

class PhylogenyRoot {
    const unsigned int mId;
    static unsigned int msNextId; // Increments from one.
    PhylogenyNode* mpRoot;

  public:
    //Constructor
    PhylogenyRoot(Cell*); // Create root node with cell.
    PhylogenyNode* Root();
    //Destructor
    ~PhylogenyRoot(); // Create root node with cell.
};


class PhylogenyNode {
  const unsigned long mId;
  static unsigned long msNextId; // Increments from one.

  Cell* mpCell;
  int mTypeId;
  PhylogenyNode* mpUp;
  PhylogenyNode* mpLeft;
  PhylogenyNode* mpRight;

  unsigned int mGeneration;
  unsigned int mNumMutsGeneration;

  public:
    // Constructors:
    PhylogenyNode(Cell*);
    PhylogenyNode(Cell*, PhylogenyNode*);
    //Destructor
    ~PhylogenyNode();
    // Getters:
    unsigned long Id();
    unsigned int Generation();
    unsigned int NumMutations();
    Cell* AssociatedCell();
    PhylogenyNode* UpNode();
    PhylogenyNode* LeftNode();
    PhylogenyNode* RightNode();
    std::vector <PhylogenyNode*> NodeAncestry();
    int SampleNode(std::vector <int>&, std::vector <double>&, double, int, int);
    int SampleNode(std::ofstream&, double, int, int);
    int TypeId();

    // Setters:
    void AddNewMutations(int);
    void AssociatedCell(Cell*);
    void UpNode(PhylogenyNode*);
    void LeftNode(PhylogenyNode*);
    void RightNode(PhylogenyNode*);
    void TypeId(int);

    // Output:
    void Print();
    void PrintAncestry();

};

#endif // PHYLOGENY_H
