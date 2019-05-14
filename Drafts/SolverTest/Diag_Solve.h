#ifndef __DIAGONAL_SOLVE_H__
#define __DIAGONAL_SOLVE_H__

#include "./Macroblock_Data_Helper.h"
#include "./Hierarchical_Cholesky.h"
#include "./macro_block_solver.h"
#include <mkl.h>

namespace Hierarchical_Solve {
    // forward declarations
    void Print_Out_Result(MACROBLOCK::DOF_INTERIOR_DATA* x);

    // macroblock solve
    template <class T, size_t Nx, size_t Ny, size_t Nz>
        void Diag_Solve(
            T (&x_variable)[Nx][Ny][Nz][3],
            const T (&y_variable)[Nx][Ny][Nz][3],
            const T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3],
            Macroblock_Utilities::Hierarchical_Cholesky<T> &hierarchical_cholesky,
            MACROBLOCK::MATRIX_DATA* matrix_data,
            MACROBLOCK::MATRIX_STRUCTURE* matrix_structure,
            MACROBLOCK::DOF_INTERIOR_DATA* x,
            MACROBLOCK::DOF_SCRATCH_DATA* scratch,
            int test_rank,
            std::vector<int>& mb_axes
            )
    {
        using namespace Macroblock_Utilities;
        constexpr int d = 3;
//        Print_Out_Result(x+test_rank);
        //Copy_In_DOF(x, y_variable, test_rank, matrix_structure, mb_axes);
        //std::cout<<"x = "<<x<<std::endl;
        /* Copy_In_Diags(matrix_data, hierarchical_cholesky, test_rank, matrix_structure, mb_axes); */
        /* Copy_In_Transfer(matrix_data, K_matrix, test_rank, matrix_structure, mb_axes); */
        MACROBLOCK::Solve_Level<4,true>( *matrix_structure,*matrix_data,*x, *x, *scratch);
        //Copy_Out_DOF(x_variable, x, test_rank, matrix_structure, mb_axes);
        //can print out solve result for test
        //Print_Out_Result(x);
    }

    template <class T, size_t Nx, size_t Ny, size_t Nz>
        void Diag_Solve(
            T (&x_variable)[Nx][Ny][Nz][3],
            const T (&y_variable)[Nx][Ny][Nz][3],
            const T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3],
            Macroblock_Utilities::Hierarchical_Cholesky<T> &hierarchical_cholesky,
            int test_rank,
            int test_level
            ) {
        //LOG::SCOPE scope("Schur Solve level "+std::to_string(test_level));
        using namespace Macroblock_Utilities;
        constexpr int d = 3;
        auto mapping = Aggregate_Topology::Aggregate_Mapping(test_level, test_rank, {3,3,3});
        std::array<int, d> step{mapping[0][0]>mapping[1][0]?-1:1,mapping[0][1]>mapping[1][1]?-1:1,mapping[0][2]>mapping[1][2]?-1:1};
        //std::cout<<"step = ["<<step[0]<<" "<<step[1]<<" "<<step[2]<<"]"<<std::endl;

        const int Matrix_K = hierarchical_cholesky.agg_topologies[test_level-1]->Ninterface;
        //std::cout<<"Matrix_K = "<<Matrix_K<<std::endl;
        float* RHS = 0;
        int status=posix_memalign((void**)(&RHS),64,Matrix_K*d*sizeof(T));
        if(status) { LOG::cout<<"posix_memalign failed with error value = "<<status<<std::endl; exit(1); }
        // grab RHS data
        auto& coordinates = hierarchical_cholesky.agg_topologies[test_level-1]->coordinates;
        for (int i=0; i<Matrix_K; i++) {
            const std::array<int,d> coord{coordinates[i][0]*step[0]+mapping[0][0],coordinates[i][1]*step[1]+mapping[0][1],coordinates[i][2]*step[2]+mapping[0][2]};
            //std::cout<<"coord = ["<<coord[0]<<" "<<coord[1]<<" "<<coord[2]<<"]"<<std::endl;
            for (int v=0; v<d; v++)
                RHS[i*3+v] = y_variable[coord[0]-1][coord[1]-1][coord[2]-1][v];
        }
        cblas_stpsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, Matrix_K*d, hierarchical_cholesky.agg_schur[test_level-1][test_rank], RHS, 1);
        cblas_stpsv(CblasColMajor, CblasLower, CblasTrans, CblasNonUnit, Matrix_K*d, hierarchical_cholesky.agg_schur[test_level-1][test_rank], RHS, 1);
        //copy back to x_variable
        for (int i=0; i<Matrix_K; i++) {
            const std::array<int,d> coord{coordinates[i][0]*step[0]+mapping[0][0],coordinates[i][1]*step[1]+mapping[0][1],coordinates[i][2]*step[2]+mapping[0][2]};
            for (int v=0; v<d; v++)
                x_variable[coord[0]-1][coord[1]-1][coord[2]-1][v] = RHS[i*3+v];
        }
        free(RHS);
    }



}
#endif //__DIAGONAL_SOLVE_H__
