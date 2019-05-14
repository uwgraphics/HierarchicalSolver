//#####################################################################
// Copyright (c) 2018, Eftychios Sifakis, Qisi Wang
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Aggregate_Topology_h__
#define __Aggregate_Topology_h__

#include "Macroblock_Utilities.h"

#include <vector>
#include <array>
#include <unordered_map>

namespace Macroblock_Utilities{
    // forward declaration
    struct Subdomain_Topology;
//#####################################################################
// Class Aggregate_Topology
//#####################################################################
    struct Aggregate_Topology
    {
        static constexpr int d=3;
        using Coord_Type = std::array<int,d>;
        using Hashtable_Type = std::unordered_map<Coord_Type,int,ArrayHasher<d>>;
        using Coord_Array_Type = std::vector<Coord_Type>;

        int level;
        Coord_Type aggregate_size;
        const int Nvariables;
        const int Ninterface;
        Hashtable_Type index;
        Coord_Array_Type coordinates; // a reverse mapping for filling in the data, only necesary for subdomain
        // TODO: Figure out if this is really necessary
        std::vector<std::array<int,2>> map;

//#####################################################################
        void Initialize_Coordinates();
        /* static Aggregate_Topology* Aggregate(const Subdomain_Topology& subdomain); */
        /* static Aggregate_Topology* Aggregate(const Aggregate_Topology& child); */

        static int Aggregate_Axis(const int level);
        static std::array<int,d> Aggregation_In_Each_Dimension(const int level);


        static Aggregate_Topology* Aggregate(const Subdomain_Topology& subdomain, const int depth);
        static Aggregate_Topology* Aggregate_Morton(const Subdomain_Topology& subdomain, const int depth);

        static Aggregate_Topology* Aggregate(const Aggregate_Topology& child, const int depth);
        static Aggregate_Topology* Aggregate_Morton(const Aggregate_Topology& child, const int depth);

        // the mapping depend on aggregation axis, thus is put in this class
        static void Get_sbd_Mapping_Morton(std::array<int, d>& min_corner, std::array<int, d>& Sstep, const std::array<int, d>& sbd_size, const int rank);

        static std::array<std::array<int,3>,2> Aggregate_Mapping(int level, int rank, const std::array<int, 3>& sbd_size);
        static std::array<std::array<int,3>,2> Level_Mapping(int level, int rank, const std::array<int, 3>& sbd_size);

        static std::array<std::array<int,3>,2> Level_Range(int level, int rank, int depth);
        static void Get_sbd_Mapping(std::array<int, d>& min_corner, std::array<int, d>& Sstep, const std::array<int, d>& sbd_size, const int rank); // not needed
//#####################################################################
    };
}
#endif
