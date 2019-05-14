//#####################################################################
// Copyright (c) 2018, Eftychios Sifakis, Qisi Wang
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Hierarchical_Cholesky_h__
#define __Hierarchical_Cholesky_h__

#include <vector>
#include <array>

#include "./Macroblock_Utilities/BlockedCSCSymmetricMatrix3.h"
#include "./Macroblock_Utilities.h"

#include "./Subdomain_Topology.h"
#include "./Subdomain_Data.h"
#include "./Aggregate_Topology.h"
#include "./Aggregate_Data.h"

namespace{
    template<class T>
        struct Blocked_Array_3D {
            using Matrix3_ptr        = T (*)[9];
            using Matrix3_type       = T (&)[9];
            using Matrix3_const_type = const T (&)[9];

            static constexpr int d = 3;
            size_t Nx, Ny, Nz;
            T* data;

             inline Matrix3_type operator()(size_t i, size_t j, size_t k, int di, int dj, int dk) {
                return reinterpret_cast<Matrix3_ptr>(data)[(i*Ny*Nz+j*Nz+k)*d*d*d+di*d*d+dj*d+dk];
            }

             inline const Matrix3_const_type operator()(size_t i, size_t j, size_t k, int di, int dj, int dk)  const{
                return reinterpret_cast<const Matrix3_ptr>(data)[(i*Ny*Nz+j*Nz+k)*d*d*d+di*d*d+dj*d+dk];
            }
        };
}

namespace Macroblock_Utilities{
//#####################################################################
// Class Hierarchical_Cholesky
//#####################################################################
template<class T>
class Hierarchical_Cholesky
{
public:
    static constexpr int d = 3;
    using Coord_Type = std::array<int,d>;

    Coord_Type domain_size;
    T* data;

    int depth;
    Subdomain_Topology sbd_topology;
    std::vector<Aggregate_Topology*> agg_topologies;

    std::vector<Subdomain_Data<T>> sbd_schur;         // L matrix for subdomain
    std::vector<std::vector<T*>>   agg_schur;       // shcur complements for interface
    CSC_Topology sbd_csc_topology;                       // topology for schur of subdomain
    //private:
    Subdomain_Data<T> sbd_workspace;
    std::vector<Aggregate_Data<T>> agg_workspace;
//#####################################################################
public:
Hierarchical_Cholesky(int Nx, int Ny, int Nz)
    : domain_size{Nx,Ny,Nz}
    {}
    void Initialize();
    void Fill_Subdomain_CSC_Data(BlockedCSCSymmetricMatrix3<T>& CSC, const int rank);
    void Allocate_Space_For_Schur_Complements();
    void Compute_Schur_Complements();

    //private:
    void Compute_Schur_Complements(int level, int rank);
    void Accumulate_to_Next_Level(int level, const int child_index);
    void Accumulate_to_Next_Level_Naive(int level, const int child_index);
    void Extract_Schur(const int level, const int rank);
    inline void Cholesky_on_Aggregate(const int level);
    // void Prepare_Leaf_Node_Data(const int level);
//#####################################################################
};
}
#endif
