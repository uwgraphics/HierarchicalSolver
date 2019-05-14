#ifndef __HIERARCHICAL_LEVEL_SOLVE_H__
#define __HIERARCHICAL_LEVEL_SOLVE_H__

#include "./Macroblock_Data_Helper.h"
#include "./Hierarchical_Cholesky.h"
#include "./Diag_Solve.h"
namespace {
    template <class T, size_t Nx, size_t Ny, size_t Nz>
    void Visualize_Array(const T (&x_variable)[Nx][Ny][Nz][3]) {
        // for debug purpose, used to visualize the data for 0th dim
    constexpr int d=3;
    for(int i=0;i<Nx;i++)
        for(int j=0;j<Ny;j++) {
            for (int v=0; v<j; v++)
                std::cout<<" ";
            for(int k=0;k<Nz;k++) {
                int v = 0;
                    std::cout<<std::setw(5)<<std::setprecision(3)<<x_variable[i][j][k][v];
            }
            std::cout<<std::endl;
        }
}
}
namespace Hierarchical_Solve {
    template <class T, size_t Nx, size_t Ny, size_t Nz>
        void Solve(
            T (&x_variable)[Nx][Ny][Nz][3],
            const T (&y_variable)[Nx][Ny][Nz][3],
            const T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3],
            const int (&level_map)[Nx][Ny][Nz],
            Macroblock_Utilities::Hierarchical_Cholesky<T> &hierarchical_cholesky,
            MACROBLOCK::MATRIX_DATA* matrix_data,
            MACROBLOCK::MATRIX_STRUCTURE* matrix_structure,
            MACROBLOCK::DOF_INTERIOR_DATA* x,
            MACROBLOCK::DOF_SCRATCH_DATA* scratch,
            std::vector<int>& mb_axes,
            int level
            ) {

        constexpr int d = 3;
        constexpr int mb_level = 4;
        using variable_type = T (&)[Nx][Ny][Nz][d];
        const int& depth = hierarchical_cholesky.depth;
        if (level>mb_level) {
            T *temp_raw = new T [Nx*Ny*Nz*d]();
            variable_type temp = reinterpret_cast<variable_type>(*temp_raw);
// 1 temp_i = inv(K_ii)*y_i
            Solve(temp,y_variable,K_matrix,level_map,hierarchical_cholesky,matrix_data, matrix_structure, x, scratch, mb_axes,level-1); // should this be temp or y?

// 2 temp_r=K_ri*temp_i; -temp_r+=y_r
            Transfer(temp, K_matrix, level_map, level);
#pragma omp parallel for
            for (int i=0; i<Nx; i++)
            for (int j=0; j<Ny; j++)
            for (int k=0; k<Nz; k++)
                if (level_map[i][j][k] == level)
                    for (int v=0; v<d; v++)
                        temp[i][j][k][v] = y_variable[i][j][k][v]-temp[i][j][k][v];

// 3 temp_r = inv(Sig)*temp_r
            for (int test_rank=0; test_rank<1<<(depth-level); test_rank++) {
                Diag_Solve(x_variable, temp, K_matrix, hierarchical_cholesky, test_rank, level);// need to change temp to x_variable later
            }
#pragma omp parallel for
            for (int i=0; i<Nx; i++)
            for (int j=0; j<Ny; j++)
            for (int k=0; k<Nz; k++)
                if (level_map[i][j][k] == level)
                    for (int v=0; v<d; v++)
                        temp[i][j][k][v] = x_variable[i][j][k][v];

            // 4 temp_i=K_ir*temp_r; -temp_i+=y_i
            TransferBack(temp, K_matrix, level_map, level);
#pragma omp parallel for
            for (int i=0; i<Nx; i++)
            for (int j=0; j<Ny; j++)
            for (int k=0; k<Nz; k++)
                if (level_map[i][j][k] < level)
                    for (int v=0; v<d; v++)
                        temp[i][j][k][v] = y_variable[i][j][k][v]-temp[i][j][k][v];

            // 5 x_i = inv(K_ii)*temp_i; x = temp
            Solve(x_variable,temp,K_matrix,level_map,hierarchical_cholesky,matrix_data, matrix_structure, x, scratch, mb_axes,level-1); // should this be temp or y?
            delete[] temp;
        } else if (level==mb_level) {
            using namespace Macroblock_Utilities;
            // solve all subdomains
            {
                //LOG::SCOPE scope("MB Solve");
                #pragma omp parallel for
                for (int test_rank=0; test_rank<1<<(depth-mb_level); test_rank++) {
                    Clear_DOF(x+test_rank);
                    Copy_In_DOF(x+test_rank, y_variable, test_rank, matrix_structure, mb_axes);
                    Diag_Solve(x_variable, y_variable, K_matrix, hierarchical_cholesky, matrix_data+test_rank, matrix_structure, x+test_rank, scratch+test_rank, test_rank, mb_axes);
                    Copy_Out_DOF(x_variable, x+test_rank, test_rank, matrix_structure, mb_axes);
                }
            }
        } else {
            std::cerr<<"level = "<<level<<" need to be greater than "<<mb_level<<std::endl;
            exit(1);
        }

    }
 #if 1
    template <class T, size_t Nx, size_t Ny, size_t Nz>
        void TransferBack(T (&temp)[Nx][Ny][Nz][3], const T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3], const int (&level_map)[Nx][Ny][Nz], int level){
        //LOG::SCOPE scope("Transferback level "+std::to_string(level));
        constexpr int d = 3;
#pragma omp parallel for
        for (int i=0; i<Nx; i++)
        for (int j=0; j<Ny; j++)
        for (int k=0; k<Nz; k++)
            if (level_map[i][j][k] < level)
                for (int v=0; v<d; v++) {
                    temp[i][j][k][v] = (T) 0;
                    for (int di=-1; di<=1; di++)
                    for (int dj=-1; dj<=1; dj++)
                    for (int dk=-1; dk<=1; dk++)
                        if ((i+di>=0) && (i+di<Nx) &&(j+dj>=0) && (j+dj<Ny) &&(k+dk>=0) && (k+dk<Nz) && (level_map[i+di][j+dj][k+dk]==level))
                            for (int w=0; w<d; w++)
                                temp[i][j][k][v] += K_matrix[i][j][k][di+1][dj+1][dk+1][v][w] * temp[i+di][j+dj][k+dk][w];
                        }

        }
#else
    // original implementation
    template <class T, size_t Nx, size_t Ny, size_t Nz>
        void TransferBack(T (&temp)[Nx][Ny][Nz][3], const T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3], const int (&level_map)[Nx][Ny][Nz], int level){
        //LOG::SCOPE scope("Transfer back level "+std::to_string(level));
        constexpr int d = 3;
#pragma omp parallel for
            for (int i=0; i<Nx; i++)
            for (int j=0; j<Ny; j++)
            for (int k=0; k<Nz; k++)
                if (level_map[i][j][k] < level)
                    for (int v=0; v<d; v++)
                        temp[i][j][k][v] = (T) 0;

            for (int i=0; i<Nx; i++)
            for (int j=0; j<Ny; j++)
            for (int k=0; k<Nz; k++)
                if (level_map[i][j][k] == level)
                    for (int v=0; v<d; v++)
                        for (int di=-1; di<=1; di++)
                        for (int dj=-1; dj<=1; dj++)
                        for (int dk=-1; dk<=1; dk++)
                            if ((i+di>=0) && (i+di<Nx) &&(j+dj>=0) && (j+dj<Ny) &&(k+dk>=0) && (k+dk<Nz) && (level_map[i+di][j+dj][k+dk]<level))
                                for (int w=0; w<d; w++)
                                    temp[i+di][j+dj][k+dk][v] += K_matrix[i][j][k][di+1][dj+1][dk+1][w][v] * temp[i][j][k][w];
        }
    #endif

    template <class T, size_t Nx, size_t Ny, size_t Nz>
        void Transfer(T (&temp)[Nx][Ny][Nz][3] ,const T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3],  const int (&level_map)[Nx][Ny][Nz], int level){
        //LOG::SCOPE scope("Transfer level "+std::to_string(level));
        constexpr int d = 3;
#pragma omp parallel for
        for (int i=0; i<Nx; i++)
        for (int j=0; j<Ny; j++)
        for (int k=0; k<Nz; k++)
            if (level_map[i][j][k] == level)
                for (int v=0; v<d; v++) {
                    temp[i][j][k][v] = (T) 0;
                    for (int di=-1; di<=1; di++)
                    for (int dj=-1; dj<=1; dj++)
                    for (int dk=-1; dk<=1; dk++)
                        if ((i+di>=0) && (i+di<Nx) &&(j+dj>=0) && (j+dj<Ny) &&(k+dk>=0) && (k+dk<Nz) && (level_map[i+di][j+dj][k+dk]<level))
                            for (int w=0; w<d; w++)
                                temp[i][j][k][v] += K_matrix[i][j][k][di+1][dj+1][dk+1][v][w] * temp[i+di][j+dj][k+dk][w];
                }
        }
}
#endif//__HIERARCHICAL_LEVEL_SOLVE_H__
