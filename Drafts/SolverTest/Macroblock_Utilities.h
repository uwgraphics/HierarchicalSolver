//#####################################################################
// Copyright (c) 2018, Eftychios Sifakis, Qisi Wang
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Macroblock_Utilities_h__
#define __Macroblock_Utilities_h__

#include <array>
#include <vector>
#include <unordered_map>

namespace Macroblock_Utilities{
//#####################################################################
// Class ArrayHasher
//#####################################################################

struct CSC_Topology {
    int n;
    std::vector<int> offsets;
    std::vector<int> row;
};

//#####################################################################
// Class ArrayHasher
//#####################################################################

template<int d>
struct ArrayHasher {
    std::size_t operator()(const std::array<int,d>& a) const {
        std::size_t h = 0;
        for (auto e : a) h ^= std::hash<int>{}(e)  + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};



//#####################################################################
void Fill_Dense_CSC_Topology(CSC_Topology& topology, const int n);
//#####################################################################
}
#endif
