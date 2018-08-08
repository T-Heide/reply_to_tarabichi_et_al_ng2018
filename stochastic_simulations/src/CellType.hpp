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


#ifndef CELLTYPE_H
#define CELLTYPE_H

// Forward declerations: ///////////////////////////////////////////////////////
class Phylogeny_Node;
class Cell;
class Universe;

// Includes: ///////////////////////////////////////////////////////////////////
#include <vector>

// CellType ///////////////////////////////////////////////////////////////////

class CellType {
    // Identifiers:
    const unsigned int mId;
    static unsigned int msNextId; // Increments from one.

    // Properties:
    double mBirthRate;
    double mAlpha;
    double mMu;

    // Member specific:
    unsigned long mNumMembers;
    std::vector<Cell*> mpMembers;

  public:
    // Constructors:
    CellType(double, double, double);

    // Destructor:
    ~CellType();

    // Getter functions:
    unsigned int Id(){return mId;}
    double Birthrate(){return mBirthRate;}
    double Alpha(){return mAlpha;};
    double Mu(){return mMu;};
    unsigned long NumMembers() {return mNumMembers;}
    Cell* RandomMember();

    // Setter functions:
    void RegisterMember(Cell*);
    void DeregisterMember(Cell*);
    void Alpha(double new_alpha){mAlpha=new_alpha;};

    // Print functions:
    void Print();
    void PrintAllMembers();
};

#endif // CELLTYPE_H
