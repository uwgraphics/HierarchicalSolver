#define TEST_RESULT 1

#include "../Hierarchical_Solve.h"
#include "../Diag_Solve.h"
#include "../Hierarchical_Cholesky.h"
#include "../macro_block_solver.h"
#include "../Macroblock_Data_Helper.h"

#include "../Common_Utilities.h"

#include <algorithm> // fill
#include <iomanip> // setw
#include <iostream> // cout

#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>

template <class T, size_t Nx, size_t Ny, size_t Nz>
void Compare_To_Reference(const T (&x_ref)[Nx][Ny][Nz][3], const T (&x_variable)[Nx][Ny][Nz][3]);
template <class T, size_t Nx, size_t Ny, size_t Nz>
void IdentityForInterfaces(T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3], const std::array<std::array<int,3>,2>& range);
template <class T, size_t Nx, size_t Ny, size_t Nz>
void ClearTransferForNodes(T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3], const std::array<std::array<int,3>,2>& range);
template <class T, size_t Nx, size_t Ny, size_t Nz>
void ClearTransferForInterfaces(T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3], const std::array<std::array<int,3>,2>& range);
template <class T, size_t Nx, size_t Ny, size_t Nz>
void IdentityForNode(T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3], const std::array<int,3>& coord);
template <class T, size_t Nx, size_t Ny, size_t Nz>
void ClearTransferForNode(T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3], const std::array<int,3>& coord);
void Print_Indices(const std::unordered_map<std::array<int, 3>, int, Macroblock_Utilities::ArrayHasher<3>>& indices, const std::array<int, 3>& size);
template <class T, size_t Nx, size_t Ny, size_t Nz>
void Visualize_Error(const std::array<std::array<int,3>,2>& range, const T (&x_ref)[Nx][Ny][Nz][3], const T (&x_variable)[Nx][Ny][Nz][3]);
template <class T, size_t Nx, size_t Ny, size_t Nz>
void Visualize_Error(const T (&x_ref)[Nx][Ny][Nz][3], const T (&x_variable)[Nx][Ny][Nz][3]);

template <class T, size_t Nx, size_t Ny, size_t Nz>
T L_Inf_Distance(const T (&x_ref)[Nx][Ny][Nz][3], const T (&x_variable)[Nx][Ny][Nz][3]);

int main(int argc, char** argv ){
    using namespace Macroblock_Utilities;
    using namespace Hierarchical_Solve;
    using namespace PhysBAM;

    bool rK,rX;

        LOG::Initialize_Logging();
        PARSE_ARGS parse_args;
        parse_args.Add_Option_Argument("-rK", "Use random SPD stiffness mat");
        parse_args.Add_Option_Argument("-rX", "Use random solution vector");
        parse_args.Parse(argc,argv);
        rK = parse_args.Get_Option_Value("-rK");
        rX = parse_args.Get_Option_Value("-rX");


    constexpr size_t Nx=63;
    constexpr size_t Ny=63;
    constexpr size_t Nz=63;

    constexpr size_t d=3;
    constexpr int mb_level = 4;
    constexpr std::array<int,d> sbd_size{3,3,3};
    std::cout<<"size ["<<Nx<<", "<<Ny<<", "<<Nz<<"]"<<std::endl;;
    using T = float;
    using matrix_type = T (&)[Nx][Ny][Nz][3][3][3][d][d];

    int depth = 0;
    for (int b=Nx>>2; b; b>>=1, depth++);
    for (int b=Ny>>2; b; b>>=1, depth++);
    for (int b=Nz>>2; b; b>>=1, depth++);
    std::cout<<"depth = "<<depth<<std::endl;

    std::vector<int> axes(depth), mb_axes(depth-mb_level);
    for (int level=1; level<=depth; level++)
        axes[level-1] = Aggregate_Topology::Aggregate_Axis(level);
    for (int l=0; l<mb_axes.size(); l++)
        mb_axes[l] = axes[l+mb_level];

    // create stiffness matrix
    T* K_matrix_raw = new T [Nx*Ny*Nz*3*3*3*d*d];
    std::cout<<"K_mat size = "<<Nx*Ny*Nz*3*3*3*d*d<<" float"<<std::endl;
    matrix_type K_matrix = reinterpret_cast<matrix_type>(*K_matrix_raw);
    Populate_Stiffness_Matrix(K_matrix,rK);

    using level_map_type = int (&)[Nx][Ny][Nz];
    int* level_map_raw = new int [Nx*Ny*Nz];
    level_map_type level_map = reinterpret_cast<level_map_type>(*level_map_raw);
    std::fill(&level_map[0][0][0],&level_map[0][0][0]+Nx*Ny*Nz,0);
    for (int l=mb_level+1; l<=depth; l++)
        for (int r=0; r<1<<(depth-l); r++) {
            auto level_mapping = Aggregate_Topology::Level_Mapping(l, r, {3,3,3});
            const std::array<std::array<int,d>,2> level_range {
                min(level_mapping[0][0],level_mapping[1][0]),min(level_mapping[0][1],level_mapping[1][1]),min(level_mapping[0][2],level_mapping[1][2]),
                    max(level_mapping[0][0],level_mapping[1][0]),max(level_mapping[0][1],level_mapping[1][1]),max(level_mapping[0][2],level_mapping[1][2])};
            for (int i=level_range[0][0]; i<=level_range[1][0]; i++)
                for (int j=level_range[0][1]; j<=level_range[1][1]; j++)
                    for (int k=level_range[0][2]; k<=level_range[1][2]; k++) {
                        level_map[i-1][j-1][k-1] = l;
                    }
        }


    Hierarchical_Cholesky<T> hierarchical_cholesky{Nx, Ny, Nz};
    hierarchical_cholesky.data = K_matrix_raw;
    hierarchical_cholesky.Initialize();
    hierarchical_cholesky.Allocate_Space_For_Schur_Complements();

    // cholesky with hierarchical ordering right now...
    {
        LOG::SCOPE scope("cholesky");
        hierarchical_cholesky.Compute_Schur_Complements();
    }

#if 1
    using variable_type = T (&)[Nx][Ny][Nz][d];
    T *x_variable_raw = new T [Nx*Ny*Nz*d];
    T *y_variable_raw = new T [Nx*Ny*Nz*d];
    T *x_ref_raw = new T [Nx*Ny*Nz*d];

    variable_type x_variable = reinterpret_cast<variable_type>(*x_variable_raw);
    variable_type y_variable = reinterpret_cast<variable_type>(*y_variable_raw);
    variable_type x_ref = reinterpret_cast<variable_type>(*x_ref_raw);

    Fill_Variable<T, d, Nx, Ny, Nz>(x_ref, 1, rX);
    Fill_Variable<T, d, Nx, Ny, Nz>(y_variable, 2, false);

    Multiply<T, d, Nx, Ny, Nz>(K_matrix, x_ref, y_variable);   // Kx=y
    memcpy(x_variable_raw,y_variable_raw,Nx*Ny*Nz*d*sizeof(T));
    MACROBLOCK::MATRIX_STRUCTURE* matrix_structure;
    MACROBLOCK::MATRIX_DATA* matrix_data;
    MACROBLOCK::DOF_INTERIOR_DATA* x;
    MACROBLOCK::DOF_SCRATCH_DATA* scratch;
    if(posix_memalign((void**)(&(x)), 64,
                      (1<<(hierarchical_cholesky.depth-mb_level))*sizeof( MACROBLOCK::DOF_INTERIOR_DATA))){
        LOG::cout << "Could not allocate aligned memory... exiting..." << std::endl;
        exit(1);
    }
    if(posix_memalign((void**)(&(scratch)), 64,
                      (1<<(hierarchical_cholesky.depth-mb_level))*sizeof( MACROBLOCK::DOF_SCRATCH_DATA))){
        LOG::cout << "Could not allocate aligned memory... exiting..." << std::endl;
        exit(1);
    }
    if(posix_memalign((void**)(&(matrix_data)), 64,
                      (1<<(hierarchical_cholesky.depth-mb_level))*sizeof( MACROBLOCK::MATRIX_DATA))){
        LOG::cout << "Could not allocate aligned memory... exiting..." << std::endl;
        exit(1);
    }
    if(posix_memalign((void**)(&(matrix_structure)), 64,
                      sizeof( MACROBLOCK::MATRIX_STRUCTURE))){
        LOG::cout << "Could not allocate aligned memory... exiting..." << std::endl;
        exit(1);
    }

    MACROBLOCK::Initialize_Matrix_Structure( *matrix_structure );
    // fill in data for matrix_data
    for (int rank=0; rank<1<<(hierarchical_cholesky.depth-mb_level); rank++) {
        Copy_In_Diags(matrix_data+rank, hierarchical_cholesky, rank, matrix_structure, mb_axes);
        Copy_In_Transfer(matrix_data+rank, K_matrix, rank, matrix_structure, mb_axes);
    }

    long long int num = (sizeof( MACROBLOCK::MATRIX_DATA)+sizeof( MACROBLOCK::MATRIX_STRUCTURE))/4+(sizeof( MACROBLOCK::DOF_INTERIOR_DATA)+sizeof( MACROBLOCK::DOF_SCRATCH_DATA))/4;
    std::cout<<"macroblock = "<<num<<" float"<<std::endl;
    std::cout<<"macroblock scratch = "<<sizeof( MACROBLOCK::DOF_SCRATCH_DATA)/4<<" float"<<std::endl;
    std::cout<<"macroblock dof = "<<sizeof( MACROBLOCK::DOF_INTERIOR_DATA)/4<<" float"<<std::endl;
    std::cout<<"macroblock structure = "<<sizeof( MACROBLOCK::MATRIX_STRUCTURE)/4<<" float"<<std::endl;
    std::cout<<"macroblock data = "<<sizeof( MACROBLOCK::MATRIX_DATA)/4<<" float"<<std::endl;
    std::cout<<"temp size = "<<(hierarchical_cholesky.depth-mb_level)*Nx*Ny*Nz*d<<" float"<<std::endl;
    // int i;
    // std::cin>>i;


    {
        LOG::SCOPE scope("Solve");
        Solve(x_variable,y_variable,K_matrix,level_map,hierarchical_cholesky,matrix_data,matrix_structure, x, scratch, mb_axes,hierarchical_cholesky.depth);
    }




#if TEST_RESULT
    //test for correct result
    Compare_To_Reference(x_ref,x_variable);
    std::cout<<"L_inf = "<<L_Inf_Distance(x_ref, x_variable)<<std::endl;
    // compute residual
    Multiply(K_matrix, x_variable, x_ref);
    std::cout<<"L_inf res "<<L_Inf_Distance(x_ref, y_variable)<<std::endl;
    #endif
    LOG::Finish_Logging();
    //garbage collection
    delete[] x_variable_raw;
    delete[] y_variable_raw;
    delete[] x_ref_raw;
    delete[] K_matrix_raw;
    free( matrix_data);
    free( matrix_structure);
    free( x );
    free( scratch );

    #endif
#if 1
    for (int r=0; r<hierarchical_cholesky.sbd_schur.size(); r++) {
        delete[] hierarchical_cholesky.sbd_schur[r].L_diagonal;
        delete[] hierarchical_cholesky.sbd_schur[r].L_lower;
    }
    for (int l=1; l<=depth; l++)
        for (int r=0; r<hierarchical_cholesky.agg_schur[l-1].size(); r++)
            delete hierarchical_cholesky.agg_schur[l-1][r];
#endif

    return 0;
}

std::array<std::array<int,3>,2> Level_Range(int level, int rank, std::vector<int> agg_axes) {
    // mapping from hierarchical rank and its level to the range of a level (interface or sbd)
    // agg_axes - aggregate axis from subdomains
    constexpr int d = 3;
    std::array<std::array<int, d>,2> range{0};
    if (level<=agg_axes.size()) {
        const int axis = agg_axes[level-1];
        std::array<int,d> bit_delta{1,1,1};
        std::array<int,d> bit_cnt{0,0,0};

        for (int l=level, r=rank; l<agg_axes.size(); l++, r>>=1) {
            const int v = agg_axes[l];
            range[0][v]|=r&1?bit_delta[v]:0;
            bit_delta[v]<<=1;
            bit_cnt[v]++;
        }
        // should reverse the bit traverse order
        bit_delta = {1,1,1};
        for (int v=0; v<d; v++) {
            for (;bit_cnt[v];bit_cnt[v]--) {
                if (range[0][v]&bit_delta[v])
                    range[0][v]^=(bit_delta[v]-1);
                bit_delta[v]<<=1;
            }
        }
        //std::cout<<" ax_rank : ["<<range[0][0]<<","<<range[0][1]<<","<<range[0][2]<<"]"<<std::endl;
        if (level) {
            for (int v=0; v<d ;v++) {
                if (v!=axis) {
                    range[1][v] = range[0][v]+1;
                    range[1][v] <<= 2;
                    range[0][v] <<= 2;
                } else {
                    range[0][v] <<= 1;
                    range[0][v] |= 1;
                    range[0][v] <<= 2;
                    range[1][v] = range[0][v];
                }
            }

            for (int l=1; l<level; l++) {
                for (int s=0; s<2; s++)
                    range[s][agg_axes[l-1]]<<=1;
            }
            for (int v=(axis+1)%d; v!=axis; v=(v+1)%d) {
                range[0][v] += 1;
                range[1][v] -= 1;
            }
        } else {
            for (int v=0; v<d ;v++) {
                range[1][v] = range[0][v]+1;
                range[1][v] <<= 2;
                range[0][v] <<= 2;
            }

            for (int l=1; l<=level; l++) {
                for (int s=0; s<2; s++)
                    range[s][agg_axes[l-1]]<<=1;
            }
            for (int v=0; v<d; v++) {
                range[0][v] += 1;
                range[1][v] -= 1;
            }
        }



    }
    return range;
}

template <class T, size_t Nx, size_t Ny, size_t Nz>
void IdentityForInterfaces(T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3], const std::array<std::array<int,3>,2>& range) {
    constexpr int d=3;

    for (size_t i=range[0][0]-1; i<=range[1][0]-1; i++)
        for (size_t j=range[0][1]-1; j<=range[1][1]-1; j++)
            for (size_t k=range[0][2]-1; k<=range[1][2]-1; k++) {
                for (int di=-1; di<=1; di++)
                    for (int dj=-1; dj<=1; dj++)
                        for (int dk=-1; dk<=1; dk++)
                            for (int w=0; w<d; w++)
                                for (int v=0; v<d; v++) {
                                    K_matrix[i][j][k][di+1][dj+1][dk+1][w][v] = T();
                                    if (i+di>=0&&i+di<Nx&&j+dj>=0&&j+dj<Ny&&k+dk>=0&&k+dk<Nz)
                                        K_matrix[i+di][j+dj][k+dk][-di+1][-dj+1][-dk+1][v][w] = T();
                                }
                for (int w=0; w<d; w++)
                    K_matrix[i][j][k][1][1][1][w][w] = (T)1;
            }
}

template <class T, size_t Nx, size_t Ny, size_t Nz>
void IdentityForNode(T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3], const std::array<int,3>& coord) {
    constexpr int d=3;
    const int i = coord[0]-1;
    const int j = coord[1]-1;
    const int k = coord[2]-1;
    //std::cout<<"coord = ["<<i<<" "<<j<<" "<<k<<"]"<<std::endl;
    for (int di=-1; di<=1; di++)
        for (int dj=-1; dj<=1; dj++)
            for (int dk=-1; dk<=1; dk++)
                for (int w=0; w<d; w++)
                    for (int v=0; v<d; v++) {
                        K_matrix[i][j][k][di+1][dj+1][dk+1][w][v] = T();
                        if (i+di>=0&&i+di<Nx&&j+dj>=0&&j+dj<Ny&&k+dk>=0&&k+dk<Nz)
                            K_matrix[i+di][j+dj][k+dk][-di+1][-dj+1][-dk+1][v][w] = T();
                    }
    for (int w=0; w<d; w++)
        K_matrix[i][j][k][1][1][1][w][w] = (T)1;
}

template <class T, size_t Nx, size_t Ny, size_t Nz>
void ClearTransferForInterfaces(T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3], const std::array<std::array<int,3>,2>& range) {
    // the range is assuming the range of the domain starts at [1,1,1]
    constexpr int d=3;
    for (size_t i=range[0][0]-1; i<=range[1][0]-1; i++)
        for (size_t j=range[0][1]-1; j<=range[1][1]-1; j++)
            for (size_t k=range[0][2]-1; k<=range[1][2]-1; k++) {
                for (int di=-1; di<=1; di++)
                    for (int dj=-1; dj<=1; dj++)
                        for (int dk=-1; dk<=1; dk++)
                            if (i+1+di<range[0][0]||i+1+di>range[1][0]||
                                j+1+dj<range[0][1]||j+1+dj>range[1][1]||
                                k+1+dk<range[0][2]||k+1+dk>range[1][2])
                                for (int w=0; w<d; w++)
                                    for (int v=0; v<d; v++) {
                                        K_matrix[i][j][k][di+1][dj+1][dk+1][w][v] = T();
                                        if (i+di>=0&&i+di<Nx&&j+dj>=0&&j+dj<Ny&&k+dk>=0&&k+dk<Nz)
                                            K_matrix[i+di][j+dj][k+dk][-di+1][-dj+1][-dk+1][v][w] = T();
                                    }
            }
}

template <class T, size_t Nx, size_t Ny, size_t Nz>
void ClearTransferForNodes(T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3], const std::array<std::array<int,3>,2>& range) {
    constexpr int d=3;
    for (size_t i=range[0][0]-1; i<=range[1][0]-1; i++)
        for (size_t j=range[0][1]-1; j<=range[1][1]-1; j++)
            for (size_t k=range[0][2]-1; k<=range[1][2]-1; k++) {
                for (int di=-1; di<=1; di++)
                    for (int dj=-1; dj<=1; dj++)
                        for (int dk=-1; dk<=1; dk++)
                            if (di!=0||dj!=0||dk!=0)
                                for (int w=0; w<d; w++)
                                    for (int v=0; v<d; v++) {
                                        K_matrix[i][j][k][di+1][dj+1][dk+1][w][v] = T();
                                        if (i+di>=0&&i+di<Nx&&j+dj>=0&&j+dj<Ny&&k+dk>=0&&k+dk<Nz)
                                            K_matrix[i+di][j+dj][k+dk][-di+1][-dj+1][-dk+1][v][w] = T();
                                    }
            }
}

template <class T, size_t Nx, size_t Ny, size_t Nz>
void ClearTransferForNode(T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3], const std::array<int,3>& coord) {
    constexpr int d=3;
    const int i = coord[0]-1;
    const int j = coord[1]-1;
    const int k = coord[2]-1;
    for (int di=-1; di<=1; di++)
        for (int dj=-1; dj<=1; dj++)
            for (int dk=-1; dk<=1; dk++)
                if (di!=0||dj!=0||dk!=0)
                    for (int w=0; w<d; w++)
                        for (int v=0; v<d; v++) {
                            K_matrix[i][j][k][di+1][dj+1][dk+1][w][v] = T();
                            if (i+di>=0&&i+di<Nx&&j+dj>=0&&j+dj<Ny&&k+dk>=0&&k+dk<Nz)
                                K_matrix[i+di][j+dj][k+dk][-di+1][-dj+1][-dk+1][v][w] = T();
                        }
}


void Print_Indices(const std::unordered_map<std::array<int, 3>, int, Macroblock_Utilities::ArrayHasher<3>>& indices, const std::array<int, 3>& size)
{
    // size is the block size w/o periphery

    const int offset = 5;
    for(int i=0;i<=size[0]+1;i++)
        for(int j=0;j<=size[1]+1;j++) {
            for (int v=0;v<(offset-1)*j;v++)
                std::cout<<" ";
            for(int k=0;k<=size[2]+1;k++)
            {
                if (indices.find({i,j,k}) == indices.end()) {
                    std::cout<<std::setw(offset)<<-1;
                } else
                {
                    std::cout<<std::setw(offset)<<(*indices.find({i,j,k})).second;
                }
            }
            std::cout<<std::endl;
        }
}

int linear_index(const std::array<int,3>& coord, const std::vector<int>& axes) {
    // w/ sbd size of 3x3x3 and interleaven within sbd is z1z0y1y0x1x0
    // axis = {0,1,2,2}
    constexpr int d = 3;
    std::array<int,d> coord_tmp = coord;
    int linear_index = 0;
    int s=1;
    for (int v=0; v<d; v++)
        for (int b=0; b<2; b++,s<<=1,coord_tmp[v]>>=1)
            if (coord_tmp[v]&1)
                linear_index|=s;

    for (int l=0; l<axes.size(); l++, s<<=1) {
        int v = axes[l];
        if (coord_tmp[axes[v]]&1)
            linear_index|=s;
        coord_tmp[v]>>=1;
    }
    return linear_index;
}

template <class T, size_t Nx, size_t Ny, size_t Nz>
void Compare_To_Reference(const T (&x_ref)[Nx][Ny][Nz][3], const T (&x_variable)[Nx][Ny][Nz][3]) {
    constexpr int d=3;
    using namespace PhysBAM;

    // test for correct result
    const double thresh=1e-3; // 1e-13 for double, 1e-4 for float
    std::cout<<"Testing Result with thresh = "<<thresh<<std::endl;
    for(int i=0;i<Nx;i++)
        for(int j=0;j<Ny;j++)
            for(int k=0;k<Nz;k++)
                for(int v=0;v<d;v++) {
                    if ((abs(x_ref[i][j][k][v])<thresh*thresh&&abs(x_variable[i][j][k][v])>thresh*thresh)||(abs(x_ref[i][j][k][v])>thresh*thresh&&(abs((x_variable[i][j][k][v]-x_ref[i][j][k][v])/x_ref[i][j][k][v])>thresh||std::isnan(x_variable[i][j][k][v]))))
                        std::cout<<" "<<VECTOR<int,4>(i+1,j+1,k+1,v)<<" x="<<x_variable[i][j][k][v]<<" x_ref="<<x_ref[i][j][k][v]<<"    diff="<<x_variable[i][j][k][v]-x_ref[i][j][k][v]<<"    %err = "<<(x_variable[i][j][k][v]-x_ref[i][j][k][v])/x_ref[i][j][k][v]<<std::endl;}
}

template <class T, size_t Nx, size_t Ny, size_t Nz>
T L_Inf_Distance(const T (&x_ref)[Nx][Ny][Nz][3], const T (&x_variable)[Nx][Ny][Nz][3]) {
    T norm = 0;
    constexpr int d=3;
    for(int i=0;i<Nx;i++)
    for(int j=0;j<Ny;j++)
    for(int k=0;k<Nz;k++)
        for(int v=0;v<d;v++)
            if (abs(x_variable[i][j][k][v]-x_ref[i][j][k][v])>norm)
                norm = abs(x_variable[i][j][k][v]-x_ref[i][j][k][v]);
    return norm;
}

template <class T, size_t Nx, size_t Ny, size_t Nz>
void Visualize_Error(const T (&x_ref)[Nx][Ny][Nz][3], const T (&x_variable)[Nx][Ny][Nz][3]) {
    constexpr int d=3;
    // test for correct result
    const double thresh=1e-5; // 1e-13 for double, 1e-4 for float
    std::cout<<"Testing Result with thresh = "<<thresh<<std::endl;
    for(int i=0;i<Nx;i++)
        for(int j=0;j<Ny;j++) {
            for (int v=0; v<j; v++)
                std::cout<<" ";
            for(int k=0;k<Nz;k++)
            {
                int v = 0;
                if ((abs(x_ref[i][j][k][v])<thresh*thresh&&abs(x_variable[i][j][k][v])>thresh*thresh)||(abs(x_ref[i][j][k][v])>thresh*thresh&&(abs((x_variable[i][j][k][v]-x_ref[i][j][k][v])/x_ref[i][j][k][v])>thresh||std::isnan(x_variable[i][j][k][v])))){
                    std::cout<<std::setw(5)<<std::setprecision(3)<<(x_variable[i][j][k][v]-x_ref[i][j][k][v])/x_ref[i][j][k][v]/thresh/10;
                } else {
                    std::cout<<std::setw(5)<<0;
                }
            }
            std::cout<<std::endl;
        }
}


template <class T, size_t Nx, size_t Ny, size_t Nz>
void Visualize_Error(const std::array<std::array<int,3>,2>& range, const T (&x_ref)[Nx][Ny][Nz][3], const T (&x_variable)[Nx][Ny][Nz][3]) {
    // the range is assuming the range of the domain starts at [1,1,1]
    constexpr int d=3;
    std::cout<<"visualizing range ["<<range[0][0]<<" "<<range[0][1]<<" "<<range[0][2]<<"] -> ["<<range[1][0]<<" "<<range[1][1]<<" "<<range[1][2]<<"]"<<std::endl;
    // test for correct result
    const double thresh=1e-5; // 1e-13 for double, 1e-4 for float
    std::cout<<"Testing Result with thresh = "<<thresh<<std::endl;



    for (size_t i=range[0][0]-1; i<=range[1][0]-1; i++)
        for (size_t j=range[0][1]-1; j<=range[1][1]-1; j++) {
            for (int v=0; v<j+1-range[0][1]; v++)
                std::cout<<" ";
            for (size_t k=range[0][2]-1; k<=range[1][2]-1; k++) {
                int v = 0;
                if ((abs(x_ref[i][j][k][v])<thresh*thresh&&abs(x_variable[i][j][k][v])>thresh*thresh)||(abs(x_ref[i][j][k][v])>thresh*thresh&&(abs((x_variable[i][j][k][v]-x_ref[i][j][k][v])/x_ref[i][j][k][v])>thresh||std::isnan(x_variable[i][j][k][v])))){
                    std::cout<<std::setw(5)<<std::setprecision(3)<<(x_variable[i][j][k][v]-x_ref[i][j][k][v])/x_ref[i][j][k][v]/thresh/10;
                } else {
                    std::cout<<std::setw(5)<<0;
                }
            }
            std::cout<<std::endl;
        }
}
