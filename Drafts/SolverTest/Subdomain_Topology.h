//#####################################################################
// Copyright (c) 2018, Eftychios Sifakis, Qisi Wang
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Subdomain_Topology_h__
#define __Subdomain_Topology_h__

#include "Macroblock_Utilities.h"

#include <vector>
#include <array>
#include <unordered_map>

namespace Macroblock_Utilities{
//#####################################################################
// Class Subdomain_Topology
//#####################################################################
struct Subdomain_Topology
{
    static constexpr int d=3;
    using Coord_Type = std::array<int,d>;
    using Hashtable_Type = std::unordered_map<Coord_Type,int,ArrayHasher<d>>;
    using Coord_Array_Type = std::vector<Coord_Type>;

    Coord_Type subdomain_size;
    int Nvariables;
    int Nsubdomain;
    Hashtable_Type index;
    Coord_Array_Type coordinates; // a reverse mapping for filling in the data, only necesary for subdomain
                                  // TODO: Figure out if this is really necessary

    CSC_Topology A;               // TODO: Fill out if really necessary
    CSC_Topology L;

//#####################################################################
    void Initialize_Index_3x3x3();
    void Initialize_Coordinates();
    void Initialize_Matrix_Topology();
//#####################################################################
};
}
#endif
