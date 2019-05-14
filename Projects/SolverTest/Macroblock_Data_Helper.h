#ifndef __DATA_HELPER_H__
#define __DATA_HELPER_H__
#include <algorithm>
#include "./macro_block_solver.h"
#include "./Hierarchical_Cholesky.h"
#include "./Macroblock_Utilities/Matrix4.h"
#include "./SymmetricMatrix.h"

namespace Macroblock_Utilities{
    /*
     * containing functions for moving data from global mat to macroblocks
     * assuming that the z-axis(dim2) and x-axis(dim0) are flipped in macroblocks
     * compare to global mat space.
     * all the axes parameters represent aggregation order in global system
     */
    template <class T, size_t Nx, size_t Ny, size_t Nz>
        void Copy_Out_DOF(T (&x_variable)[Nx][Ny][Nz][3], const MACROBLOCK::DOF_INTERIOR_DATA* x, const int rank, const MACROBLOCK::MATRIX_STRUCTURE* matrix_structure, const std::vector<int>& axes);
    template <class T, size_t Nx, size_t Ny, size_t Nz>
        void Copy_In_DOF(MACROBLOCK::DOF_INTERIOR_DATA* x, const T (&y_variable)[Nx][Ny][Nz][3], const int rank, const MACROBLOCK::MATRIX_STRUCTURE* matrix_structure, const std::vector<int>& axes);
    void Clear_DOF(MACROBLOCK::DOF_INTERIOR_DATA* x);
    template <class T, size_t Nx, size_t Ny, size_t Nz>
        void Copy_In_Transfer(MACROBLOCK::MATRIX_DATA* matrix_data, const T(&K_matrix)[Nx][Ny][Nz][3][3][3][3][3], const int rank, const MACROBLOCK::MATRIX_STRUCTURE* matrix_structure, const std::vector<int>& axes);
    template <class T>
        void Copy_In_Diags(MACROBLOCK::MATRIX_DATA* matrix_data, Macroblock_Utilities::Hierarchical_Cholesky<T>& hierarchical_cholesky, const int rank, const MACROBLOCK::MATRIX_STRUCTURE* matrix_structure, const std::vector<int>& axes);

    int Convert_Hierarchical_Sbd_and_Macroblock_Sbd(int sbd_nubmer);
    void Cube_Coord_Mapping(std::array<int,3>& step, std::array<int,3>& min_corner, const int rank, const int cube_number, const std::vector<int> axes);
    void Macroblock_Rank_Coord_Mapping(std::array<int,3>& step, std::array<int,3>& min_corner, std::array<int,3>& box_size, const int rank);
    int Primary_Cube(const std::array<int,3>& canonical_coord, const int cube_n);
    std::vector<int> Adjacent_Subdomain_Ranks(const int level, const int rank, const std::vector<int>& axes);
    void Inplace_Conjugate_L(Matrix4<float,float>& diag);
}

namespace Macroblock_Utilities {




    template <class T, size_t Nx, size_t Ny, size_t Nz>
        void Copy_In_Transfer(MACROBLOCK::MATRIX_DATA* matrix_data, const T(&K_matrix)[Nx][Ny][Nz][3][3][3][3][3], const int rank, const MACROBLOCK::MATRIX_STRUCTURE* matrix_structure, const std::vector<int>& axes) {
        using namespace Macroblock_Utilities;
        constexpr int d = 3;
        // generate revert mapping of dof_index
        std::vector<std::array<int,d>> reverse_index(125);
        for (int i=0; i<5; i++)
            for (int j=0; j<5; j++)
                for (int k=0; k<5; k++)
                    // exchange x and z so the coordinate is in hierarchical solve order
                    reverse_index[matrix_structure->dof_index[i][j][k]] = {k,j,i};

        // Transfer mats
        for (int l=1; l<=4; l++) {
            float* ptr = matrix_data->TransferPointer(l);
            std::fill(ptr, ptr+matrix_data->TransferSize(l)/sizeof(float),0);
        }

#if 1
        {
            for (int sbd_rank=0; sbd_rank<1<<4; sbd_rank++)
                for (int level=1; level<=4; level++) {
                    std::array<int,d> step, min_corner;
                    const int cube_n = Convert_Hierarchical_Sbd_and_Macroblock_Sbd(sbd_rank);
                    Cube_Coord_Mapping(step, min_corner,rank,cube_n,axes);
                    const int row_offset = matrix_structure->TransferColumns(level);
                    BlockedCSRMatrix3<float[16]> transfer{matrix_structure->TransferRows(level), matrix_structure->TransferColumns(level), matrix_data->TransferPointer(level), matrix_structure->TransferOffsets(level), matrix_structure->TransferColumn(level)};
                    // Transfer matrices in CSR format
                    for (int i=0; i<transfer.rows; i++) {
                        for (int jj=transfer.offsets[i]; jj<transfer.offsets[i+1]; jj++) {
                            int j = transfer.column[jj];
                            const auto& node = reverse_index[i+row_offset];
                            const auto& neighbor_node = reverse_index[j];

                            std::array<int,d> indicator = node;
                            for (int v=0; v<d; v++)
                                if (!(indicator[v]&3))
                                    indicator[v] = neighbor_node[v];
                            const std::array<int,d> global_coord {node[0]*step[0]+min_corner[0], node[1]*step[1]+min_corner[1], node[2]*step[2]+min_corner[2]};
                            const std::array<int,d> global_neighbor_coord {neighbor_node[0]*step[0]+min_corner[0], neighbor_node[1]*step[1]+min_corner[1], neighbor_node[2]*step[2]+min_corner[2]};
                            if (cube_n==Primary_Cube({indicator[2],indicator[1],indicator[0]},cube_n)) {
                                if (global_coord[0]>0&&global_coord[0]<Nx+1&&global_coord[1]>0&&global_coord[1]<Ny+1&&global_coord[2]>0&&global_coord[2]<Nz+1&&
                                    global_neighbor_coord[0]>0&&global_neighbor_coord[0]<Nx+1&&global_neighbor_coord[1]>0&&global_neighbor_coord[1]<Ny+1&&global_neighbor_coord[2]>0&&global_neighbor_coord[2]<Nz+1) {

                                    for (int v=0; v<9; v++)
                                        transfer.E(jj)[v][cube_n] = K_matrix[global_coord[0]-1][global_coord[1]-1][global_coord[2]-1][global_neighbor_coord[0]-global_coord[0]+1][global_neighbor_coord[1]-global_coord[1]+1][global_neighbor_coord[2]-global_coord[2]+1][v%d][v/d];
                                }
                            }
                        }
                    }
                }
        }
#endif
    }

    template <class T>
        void Copy_In_Diags(MACROBLOCK::MATRIX_DATA* matrix_data, Macroblock_Utilities::Hierarchical_Cholesky<T>& hierarchical_cholesky, const int mb_rank, const MACROBLOCK::MATRIX_STRUCTURE* matrix_structure, const std::vector<int>& axes) {
        using namespace Macroblock_Utilities;
        constexpr int d = 3;
        constexpr int mb_level = 4;
        // fill matrices
        // L0
        CSC_Topology& topology = hierarchical_cholesky.sbd_csc_topology;
        for (int cube_n=0; cube_n<16; cube_n++) {
            const int h_rank = Convert_Hierarchical_Sbd_and_Macroblock_Sbd(cube_n)+(mb_rank<<mb_level);
            //            std::cout<<"h_rank = "<<h_rank<<std::endl;
            BlockedCSCSymmetricMatrix3<float> L{topology.n, hierarchical_cholesky.sbd_schur[h_rank].L_diagonal, hierarchical_cholesky.sbd_schur[h_rank].L_lower, &topology.offsets[0], &topology.row[0]};
            for (int j=0; j<topology.n; j++)
                for (int v=0; v<6; v++)
                    matrix_data->L0_diagonal_entries[j][v][cube_n] = L.D(j)[v];
            for (int ii=0; ii<topology.offsets[topology.n]; ii++)
                for (int v=0; v<9; v++)
                    matrix_data->L0_offdiagonal_entries[ii][v][cube_n] = L.L(ii)[v];
        }

        float* ptr;
        float identity4x4[] = {1,0,0,0,
                               0,1,0,0,
                               0,0,1,0,
                               0,0,0,1};


        // sigma_diagonal_entries
        std::array<int,d> sizes[] ={{3,3,1},
                                    {3,1,7},
                                    {1,7,7},
                                    {1,7,7}};

#if 1
        {
            const Matrix4<float, float> identity4x4 = {1,0,0,0,
                                                       0,1,0,0,
                                                       0,0,1,0,
                                                       0,0,0,1};

            const int rank_mapping[]{2,6,0,4,3,7,1,5};
            for (int level=1; level<=4; level++) {
                for (int rank=0; rank<(1<<(4-level)); rank++) {
                    const int N = hierarchical_cholesky.agg_topologies[level-1]->Ninterface*3;
                    // copy in data in 4x4 blocked form
                    const int N_full_col = N/4;
                    const int N_remain_col = N%4;
                    int N_total_col  = matrix_structure->SigmaN(level);
                    const SymmetricMatrix<float> scalar{N,hierarchical_cholesky.agg_schur[level-1][rank+(mb_rank<<(mb_level-level))]};
                    //std::cout<<"level = "<<level<<" rank = "<<rank+(mb_rank<<(mb_level-level))<<std::endl;
                    int mapped_rank = rank_mapping[Convert_Hierarchical_Sbd_and_Macroblock_Sbd(rank)>>level]>>(level-1);
                    BlockedMatrixNXN<float> blocked{N_total_col, matrix_data->SigmaDiagonal(level,mapped_rank), matrix_data->SigmaOffDiagonal(level,mapped_rank)};

                    // diagonal entries
                    {
                        int index = 0;
                        for (int bj=0; bj<N_full_col; bj++) {
                            Matrix4<float,float> diag{0};
                            for (int j=0; j<4; j++)
                                for (int i=j; i<4; i++)
                                    diag.data[4*j+i] = scalar(bj*4+i,bj*4+j);
                            diag.Invert();
                            for (int bi=bj+1; bi<N_full_col; bi++) {
                                Matrix4<float,float> off_diag{0};
                                for (int j=0; j<4; j++)
                                    for (int i=0; i<4; i++)
                                        off_diag.data[4*j+i] = scalar(bi*4+i,bj*4+j);
                                off_diag*=diag;
                                off_diag.Store(blocked.L(index++));
                            }
                            {
                                Matrix4<float,float> off_diag{0};
                                int bi = N_full_col;
                                for (int j=0; j<4; j++)
                                    for (int i=0; i<N_remain_col; i++)
                                        off_diag.data[4*j+i] = scalar(bi*4+i,bj*4+j);
                                for (int j=0; j<4; j++)
                                    for (int i=N_remain_col; i<4; i++)
                                        off_diag.data[4*j+i] = 0.f;
                                off_diag*=diag;
                                off_diag.Store(blocked.L(index++));

                            }
                            for (int bi=N_full_col+1; bi<blocked.n; bi++, index++)
                                for (int v=0; v<16; v++)
                                    blocked.L(index)[v] = 0.f;
                            Inplace_Conjugate_L(diag);
                            diag.Store(blocked.D(bj));
                        }
                        {
                            // 4x4 matrices that are partially filled
                            int bj=N_full_col;
                            Matrix4<float,float> diag{0};
                            for (int j=0; j<N_remain_col; j++)
                                for (int i=j; i<N_remain_col; i++) {
                                    diag.data[4*j+i] = scalar(bj*4+i,bj*4+j);
                                }
                            for (int j=N_remain_col; j<4; j++) {
                                diag.data[5*j] = 1.f;
                            }
                            diag.Invert();
                            Inplace_Conjugate_L(diag);
                            diag.Store(blocked.D(bj));

                            for (int bi=N_full_col+1; bi<blocked.n; bi++, index++)
                                for (int v=0; v<16; v++)
                                    blocked.L(index)[v] = 0.f;
                        }
                        // matrices that are completely out of range
                        for (int bj=N_full_col+1; bj<blocked.n; bj++) {
                            identity4x4.Store(blocked.D(bj));
                            for (int bi=bj+1; bi<blocked.n; bi++, index++)
                                for (int v=0; v<16; v++)
                                    blocked.L(index)[v] = 0.f;
                        }
                    }
                }
            }
        }
#endif
    }

    template <class T, size_t Nx, size_t Ny, size_t Nz>
        void Copy_In_DOF(MACROBLOCK::DOF_INTERIOR_DATA* x, const T (&y_variable)[Nx][Ny][Nz][3], const int rank, const MACROBLOCK::MATRIX_STRUCTURE* matrix_structure, const std::vector<int>& axes) {
        //std::cout<<"x+test_rank = "<<x<<std::endl;

        // axes - agg order of macroblocks should be yxzyxzyx...
        constexpr int d = 3;
        constexpr int mb_level = 4;
        constexpr std::array<int,mb_level> mb_axes{2,1,0,0}; // the agg order (in canonical coord) within the macroblock

        std::vector<int> all_axes(mb_level+axes.size());
        for (int l=0; l<all_axes.size(); l++)
            if (l<mb_level)
                all_axes[l] = mb_axes[l];
            else
                all_axes[l] = axes[l-mb_level];
        const std::vector<int> adj_rank = Adjacent_Subdomain_Ranks(mb_level, 0, all_axes);

#if 1
        // right now only works for 7x7x15 domain size
        //fill data by blocks
        for (int cube_n=0; cube_n<16; cube_n++) {
            //    std::cout<<"cube_n = "<<cube_n<<std::endl;
            std::array<int,d> step, min_corner;
            Cube_Coord_Mapping(step, min_corner, rank, cube_n, axes);
            //std::cout<<"min_corner = ["<<min_corner[0]<<" "<<min_corner[1]<<" "<<min_corner[2]<<"]"<<std::endl;
            for (int i=0; i<5; i++)
                for (int j=0; j<5; j++)
                    for (int k=0; k<5; k++) {
                        const std::array<int,d> coord{min_corner[0]+i*step[0],min_corner[1]+j*step[1],min_corner[2]+k*step[2]}; // in global space
                        //std::cout<<"coord = ["<<coord[0]<<" "<<coord[1]<<" "<<coord[2]<<"]"<<std::endl;
                        if (coord[0]>0&&coord[0]<=Nx&&coord[1]>0&&coord[1]<=Ny&&coord[2]>0&&coord[2]<=Nz) {

                            // need to determine level here
                            int dof_index = matrix_structure->dof_index[k][j][i];
                            int level=5;
                            for (int l=0; l<=5; l++)
                                if (dof_index<matrix_structure->TransferColumns(l)) {
                                    level = l-1;
                                    break;
                                }
                            dof_index-=matrix_structure->TransferColumns(level);
                            auto array = reinterpret_cast<float (*)[3][16]>(x->Level(level));

                            if (level<mb_level) {
                                if (cube_n==Primary_Cube({k,j,i},cube_n))
                                    for (int v=0; v<d; v++)
                                        array[dof_index][v][cube_n] = y_variable[coord[0]-1][coord[1]-1][coord[2]-1][v];
                            } else if (level == mb_level){
                                if (std::find(adj_rank.begin(), adj_rank.end(), Convert_Hierarchical_Sbd_and_Macroblock_Sbd(cube_n))!=adj_rank.end())
                                    if (cube_n==Primary_Cube({k,j,i},cube_n))
                                        for (int v=0; v<d; v++)
                                            array[dof_index][v][cube_n] = y_variable[coord[0]-1][coord[1]-1][coord[2]-1][v];
                            }
                        }
                    }
        }
#endif
    }

    template <class T, size_t Nx, size_t Ny, size_t Nz>
        void Copy_Out_DOF(T (&x_variable)[Nx][Ny][Nz][3], const MACROBLOCK::DOF_INTERIOR_DATA* x, const int rank, const MACROBLOCK::MATRIX_STRUCTURE* matrix_structure, const std::vector<int>& axes) {
        // axes - agg order of macroblocks should be yxzyxzyx...
        // right now only works for 7x7x15 domain size
// copy back
        constexpr int d = 3;
        constexpr int mb_level = 4;
        constexpr std::array<int,mb_level> mb_axes{2,1,0,0}; // the agg order (in canonical coord) within the macroblock

        std::vector<int> all_axes(mb_level+axes.size());
        for (int l=0; l<all_axes.size(); l++)
            if (l<mb_level)
                all_axes[l] = mb_axes[l];
            else
                all_axes[l] = axes[l-mb_level];
        const std::vector<int> adj_rank = Adjacent_Subdomain_Ranks(mb_level, 0, all_axes);

        for (int cube_n=0; cube_n<16; cube_n++) {
            std::array<int,d> step, min_corner;
            Cube_Coord_Mapping(step, min_corner, rank, cube_n, axes);
            for (int i=0; i<5; i++)
                for (int j=0; j<5; j++)
                    for (int k=0; k<5; k++) {
                        const std::array<int,d> coord{min_corner[0]+i*step[0],min_corner[1]+j*step[1],min_corner[2]+k*step[2]}; // in global space
                        if (coord[0]>0&&coord[0]<=Nx&&coord[1]>0&&coord[1]<=Ny&&coord[2]>0&&coord[2]<=Nz) {
                            int dof_index = matrix_structure->dof_index[k][j][i];
                            int level=5;
                            for (int l=0; l<=5; l++)
                                if (dof_index<matrix_structure->TransferColumns(l)) {
                                    level = l-1;
                                    break;
                                }
                            dof_index-=matrix_structure->TransferColumns(level);
                            auto array = reinterpret_cast<float (*)[3][16]>(x->Level(level));
                            if (level<mb_level) {
                                if (cube_n==Primary_Cube({k,j,i},cube_n))
                                    for (int v=0; v<d; v++)
                                        x_variable[coord[0]-1][coord[1]-1][coord[2]-1][v] = array[dof_index][v][cube_n];
                            } else if (level == mb_level){
                                if (std::find(adj_rank.begin(), adj_rank.end(), Convert_Hierarchical_Sbd_and_Macroblock_Sbd(cube_n))!=adj_rank.end())
                                    if (cube_n==Primary_Cube({k,j,i},cube_n))
                                        for (int v=0; v<d; v++)
                                            x_variable[coord[0]-1][coord[1]-1][coord[2]-1][v] = array[dof_index][v][cube_n];
                            }

                        }
                    }
        }
    }

}

#endif // __DATA_HELPER_H__
