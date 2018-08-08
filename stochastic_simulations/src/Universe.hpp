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


#ifndef UNIVERSE_H
#define UNIVERSE_H

// Forward declerations: ///////////////////////////////////////////////////////
//class Phylogeny_Node;
class Cell;
class CellType;
class PhylogenyRoot;
class Shape;

// Includes: ///////////////////////////////////////////////////////////////////
#include <vector>

// Universe ////////////////////////////////////////////////////////////////////

class Universe {
    double mTime;
    unsigned int number_of_cells;
    std::vector<CellType*> mpTypes;
    std::vector<PhylogenyRoot*> mpPhylogenies;

  public:
    // Constructor:
    Universe();

    // Destructor:
    ~Universe();

    // Getter functions:
    double Time();
    class CellType* NextReaction(double*, int*);
    bool Sample(std::vector <int>&, std::vector <double>&, double, int);
    bool Sample(std::ofstream&, double, int);

    // Setter functions:
    void IncrementTimeBy(double);
    bool InsertCell(Cell*);
    bool InsertCell(Cell*, bool);
    void RegisterType(CellType *);

    // Output Functions:
    //void Print();
    //void TypesToCsvFile(std::string, int tum_id);
    //void PutNodesToStream(PhylogenyNode*, std::ofstream&);
    //void PhylogeniesToFile(std::string, int tum_id);

};

#endif // UNIVERSE_H
