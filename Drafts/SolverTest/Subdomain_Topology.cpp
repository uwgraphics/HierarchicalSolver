//#####################################################################
// Copyright (c) 2018, Eftychios Sifakis, Qisi Wang
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include "Subdomain_Topology.h"
#include <iostream>
#include <vector>
#include <set>

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>

namespace Macroblock_Utilities{
//#####################################################################
// Function Subdomain::Initialize_Index
//#####################################################################
namespace{
int Subdomain_Indices_5x5x5[5][5][5] = {
    // -1,  -1,  -1,  -1,  -1,
    //   -1,  -1,  -1,  -1,  -1,
    //     -1,  -1,  -1,  -1,  -1,
    //       -1,  -1,  -1,  -1,  -1,
    //         -1,  -1,  -1,  -1,  -1,
    // -1,  -1,  -1,  -1,  -1,
    //   -1,   0,  18,   9,  -1,
    //     -1,   6,  24,  15,  -1,
    //       -1,   3,  21,  12,  -1,
    //         -1,  -1,  -1,  -1,  -1,
    // -1,  -1,  -1,  -1,  -1,
    //   -1,  2,  20,  11,  -1,
    //     -1,  8,  26,  17,  -1,
    //       -1,  5,  23,  14,  -1,
    //         -1,  -1,  -1,  -1,  -1,
    // -1,  -1,  -1,  -1,  -1,
    //   -1,   1,  19,  10,  -1,
    //     -1,  7,  25,  16,  -1,
    //       -1,  4,  22,  13,  -1,
    //         -1,  -1,  -1,  -1,  -1,
    // -1,  -1,  -1,  -1,  -1,
    //   -1,  -1,  -1,  -1,  -1,
    //     -1,  -1,  -1,  -1,  -1,
    //       -1,  -1,  -1,  -1,  -1,
    //         -1,  -1,  -1,  -1,  -1

    100,  105,  110,  115,  120,
      101,  106,  111,  116,  121,
        102,  107,  112,  117,  122,
          103,  108,  113,  118,  123,
            104,  109,  114,  119,  124,
     80,   84,   88,   92,   96,
       68,    0,   18,    9,   52,
         71,    6,   24,   15,   55,
           65,    3,   21,   12,   49,
             77,   40,   43,   37,   61,
     81,   85,   89,   93,   97,
       69,    2,   20,   11,   53,
         72,    8,   26,   17,   56,
           66,    5,   23,   14,   50,
             78,   41,   44,   38,   62,
     82,   86,   90,   94,   98,
       67,    1,   19,   10,   51,
         70,    7,   25,   16,   54,
           64,    4,   22,   13,   48,
             76,   39,   42,   36,   60,
     83,   87,   91,   95,   99,
       74,   31,   34,   28,   58,
         75,   32,   35,   29,   59,
           73,   30,   33,   27,   57,
             79,   46,   47,   45,   63
};
}
void Subdomain_Topology::
Initialize_Index_3x3x3()
{
    subdomain_size = Coord_Type{3,3,3};
    Nvariables=125;
    Nsubdomain=27;

    //int count = 27;
    for (int i=0; i<=subdomain_size[0]+1; i++)
    for (int j=0; j<=subdomain_size[1]+1; j++)
    for (int k=0; k<=subdomain_size[2]+1; k++)
        //       if (Subdomain_Indices_5x5x5[i][j][k] != -1)
            index.insert({{i,j,k}, Subdomain_Indices_5x5x5[i][j][k]});
        // else
        //     index.insert({{i,j,k}, count++});
}

//#####################################################################
// Function Subdomain::Initialize_Coordinates
//#####################################################################
    void Subdomain_Topology::
    Initialize_Coordinates()
    {
        coordinates.resize(index.size());
        for (auto it : index)
            coordinates[it.second] = it.first;
    }

//#####################################################################
// Function Subdomain::Initialize_Matrix_Topology
//#####################################################################
    void Subdomain_Topology::
    Initialize_Matrix_Topology()
{
    using Matrix_Topology_Type = std::vector<std::set<int>>;
    const int n_sbd=subdomain_size[0]*subdomain_size[1]*subdomain_size[2];
    Matrix_Topology_Type A_tmp{n_sbd};
    const int n_block=(subdomain_size[0]+2)*(subdomain_size[1]+2)*(subdomain_size[2]+2);


    for(int node1_i=1; node1_i<=subdomain_size[0]; node1_i++)
    for(int node1_j=1; node1_j<=subdomain_size[1]; node1_j++)
    for(int node1_k=1; node1_k<=subdomain_size[2]; node1_k++)
        for(int node2_i=node1_i-1; node2_i<=node1_i+1; node2_i++)
        for(int node2_j=node1_j-1; node2_j<=node1_j+1; node2_j++)
        for(int node2_k=node1_k-1; node2_k<=node1_k+1; node2_k++) {
            int J=index[{node1_i,node1_j,node1_k}];
            int I=index[{node2_i,node2_j,node2_k}];
            //if(J>=I) A_tmp[I].insert(J); // use this to create just the 27x27 matrix
            if(I>=J) A_tmp[J].insert(I);
        }

    // TODO: fill topology for A

    // symbolic cholesky
    for (int j=0; j<A_tmp.size()-1; j++)
    {
        auto it = A_tmp[j].begin();
        ++it;
        if (it == A_tmp[j].end()) continue;     // //assert(it != A_tmp[j].end());
        const int parent = *it;
        if (parent >= n_sbd) continue;
        for (++it;it!=A_tmp[j].end(); it++)
            A_tmp[parent].insert(*it);
    }

    // fill in topology for L
    //int n_entry =
    L.n = n_block;
    L.offsets.resize(n_block+1);
    L.offsets[0] = 0;
    for (int j=0; j<n_sbd; j++)
        L.offsets[j+1] = L.offsets[j] + A_tmp[j].size() - 1;
    for (int j=n_sbd; j<n_block; j++)
        L.offsets[j+1] = L.offsets[j] + n_block-j - 1;
    //std::cout<<(L.offsets[L.n])<<std::endl;

    L.row.resize(L.offsets[L.n]);
    int index = 0;
    for (int j=0; j<n_sbd; j++)
    {
        auto it = A_tmp[j].begin();
        for (++it; it!=A_tmp[j].end();it++)
            L.row[index++] = *it;
    }
    //std::cout<<(L.row.size())<<std::endl;

    for (int j=n_sbd; j<n_block; j++)
        for (int i=j+1; i<n_block; i++)
            L.row[index++] = i;

    //std::cout<<index<<std::endl;
    //std::cout<<"i = [ ";for(auto it : L.row) std::cout << it+1<<" "; std::cout<<"];"<<std::endl;


}



//#####################################################################
}
