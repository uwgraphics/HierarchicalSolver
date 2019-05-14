#ifndef __COMMON_UTIL_H__
#define __COMMON_UTIL_H__

#endif// __COMMON_UTIL_H__

#include "../Macroblock_Utilities/BlockedCSCSymmetricMatrix3_Reference.h"
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>

namespace PhysBAM{
template <class T, int d, size_t Nx, size_t Ny, size_t Nz>
        void K_Matrix_To_Reference(const T (&K_matrix)[Nx][Ny][Nz][3][3][3][d][d], const ARRAY<int, VECTOR<int, d> >& linear_indices, BlockedCSCSymmetricMatrix3_Reference<T>& K_Blocked_Reference)
    {
        LOG::SCOPE scope("Converting K Reference");
        //####################### Populating Reference for Cholesky #################
        using T_INDEX = VECTOR<int, d>;
        using T_RANGE = RANGE<T_INDEX>;
        using CONST_LOCAL_MATRIX = const T (&)[d][d];

        const T_INDEX size(Nx,Ny,Nz);

        for (RANGE_ITERATOR<d> iterator(T_INDEX(1), size); iterator.Valid(); iterator.Next()) {
            const T_INDEX& index = iterator.Index();
            const int i = linear_indices(index);
            CONST_LOCAL_MATRIX Kd_local = reinterpret_cast<CONST_LOCAL_MATRIX>(*K_matrix[index(1)-1][index(2)-1][index(3)-1][1][1][1]);
            K_Blocked_Reference.D(i-1) = SYMMETRIC_MATRIX<T, d>(Kd_local[0][0], Kd_local[1][0], Kd_local[2][0],
                                                                Kd_local[1][1], Kd_local[2][1], Kd_local[2][2]);

            T_RANGE neighbor_range(index-1, index+1);
            neighbor_range.min_corner = T_INDEX::Componentwise_Max(neighbor_range.min_corner, T_INDEX(1));
            neighbor_range.max_corner = T_INDEX::Componentwise_Min(neighbor_range.max_corner, size);
            for (RANGE_ITERATOR<d> neighbor_iterator(neighbor_range); neighbor_iterator.Valid(); neighbor_iterator.Next()) {
                const T_INDEX& neighbor_index = neighbor_iterator.Index();
                const int j = linear_indices(neighbor_index);
                if (i>j) {
                    CONST_LOCAL_MATRIX Kl_local = reinterpret_cast<CONST_LOCAL_MATRIX>(*K_matrix[index(1)-1][index(2)-1][index(3)-1][neighbor_index(1)-index(1)+1][neighbor_index(2)-index(2)+1][neighbor_index(3)-index(3)+1]);
                    K_Blocked_Reference.L_Or_Insert(i-1,j-1) = MATRIX<T, d>(Kl_local[0][0], Kl_local[1][0], Kl_local[2][0],
                                                                            Kl_local[0][1], Kl_local[1][1], Kl_local[2][1],
                                                                            Kl_local[0][2], Kl_local[1][2], Kl_local[2][2]);}}}
    }

    void Reordering_To_BSC_Topology(const ARRAY<int, VECTOR<int, 3> >& linear_indices, const ARRAY<VECTOR<int, 3> >& index_coord, Macroblock_Utilities::CSC_Topology& topology){
        LOG::SCOPE scope("generating topology");
        using namespace Macroblock_Utilities;
        constexpr int d = 3;
        using T_INDEX = VECTOR<int,d>;
        using T_RANGE = RANGE<T_INDEX>;

        const T_INDEX& size = linear_indices.Size();
        const T_RANGE coordinate_range(T_INDEX(1), size);
        topology.offsets.resize(topology.n+1);
        topology.offsets[0] = 0;
        for (int j=0; j<topology.n; j++) {
            const T_INDEX& index = index_coord(j+1);
            T_RANGE neighbor_range(index-1, index+1);
            neighbor_range.min_corner = T_INDEX::Componentwise_Max(neighbor_range.min_corner, coordinate_range.min_corner);
            neighbor_range.max_corner = T_INDEX::Componentwise_Min(neighbor_range.max_corner, coordinate_range.max_corner);
            for (RANGE_ITERATOR<d> neighbor_iterator(neighbor_range); neighbor_iterator.Valid(); neighbor_iterator.Next()) {
                const T_INDEX& neighbor_index = neighbor_iterator.Index();
                const int i = linear_indices(neighbor_index)-1;

                if (i>=j)
                    topology.row.push_back(i);
                std::sort(topology.row.begin()+topology.offsets[j],topology.row.end());
                topology.offsets[j+1] = topology.row.size();
            }
        }
    }
}
