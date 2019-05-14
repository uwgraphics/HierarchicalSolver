#include "../Common_Utilities.h"
#include "../Macroblock_Utilities.h"
#include "../Macroblock_Utilities/BlockedCSCSymmetricMatrix3.h"
#include "../Macroblock_Utilities/BlockedCSCSymmetricMatrix3_Reference.h"
#include "../Macroblock_Utilities/BlockedCSCSymmetricMatrix3_Utilities.h"
#include "../Reordering_Utilities/REORDERING.h"
#include "../Reordering_Utilities/Utilities.h"

#include <iostream>
#include <iomanip>
//#include <assert.h>
#include <type_traits>
//#include <set>
#include <vector>
#include "mkl_pardiso.h"
#include "mkl_types.h"

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>

#include <stdio.h>
#include <stdlib.h>
#include <ctime> //time
//#include <math.h>

template<class T>
void printCSR(const MKL_INT* ia, const MKL_INT* ja,const T* a, MKL_INT n, bool onebased);

template <class T, size_t Nx, size_t Ny, size_t Nz>
void Compare_To_Reference(const T (&x_ref)[Nx][Ny][Nz][3], const T (&x_variable)[Nx][Ny][Nz][3]);

template <class T, size_t Nx, size_t Ny, size_t Nz>
T L_Inf_Distance(const T (&x_ref)[Nx][Ny][Nz][3], const T (&x_variable)[Nx][Ny][Nz][3]);

template <class IndexType>
void Block_To_Scalar_Topology(IndexType* offsets, IndexType* rows, const int* offsets_B, const int* rows_B, const int n);

int main(int argc,char *argv[])
{
    using namespace PhysBAM;
    std::cout<<sizeof(size_t)<<std::endl;
    {
        std::cout<<"solver"<<std::endl;
    // ############################################################################
    //   Set up parameters and intitialize the matrix, variables, and reordering
    // ############################################################################
    constexpr size_t Nx=63;
    constexpr size_t Ny=63;
    constexpr size_t Nz=63;
    constexpr size_t d=3;
    using index_type = long long;

    std::srand(0);
    using T = float;
    using T_INDEX = VECTOR<int,d>;
    using T_INDEX_ARRAYS = ARRAY<int,T_INDEX>;
    using T_COORD_ARRAYS = ARRAY<T_INDEX>;
    using T_RANGE = RANGE<T_INDEX>;
    using matrix_type = T (&)[Nx][Ny][Nz][3][3][3][d][d];

    const T_INDEX size(Nx,Ny,Nz);

    LOG::Initialize_Logging();
    bool rK,rX;
    {
        PARSE_ARGS parse_args;
        parse_args.Add_Option_Argument("-rK", "Use random SPD stiffness mat");
        parse_args.Add_Option_Argument("-rX", "Use random solution vector");
        parse_args.Parse(argc,argv);
        rK = parse_args.Get_Option_Value("-rK");
        rX = parse_args.Get_Option_Value("-rX");
    }

    REORDERING<d> reordering(T_INDEX(Nx,Ny,Nz),T_INDEX(3));
    reordering.use_queue = true;
    {
        // setup reordering so it conform with macroblock aggregate order (xyzxyzzyxzyxzyx...)
        char* cut_c = new char[reordering.depth+4];
        cut_c[reordering.depth] = 'z';
        cut_c[reordering.depth+1] = 'y';
        cut_c[reordering.depth+2] = 'x';
        cut_c[reordering.depth+3] = '\0';
        int clamped = reordering.depth>3?3:reordering.depth;
        switch (clamped) {
        case 3:
            cut_c[reordering.depth-3] = 'z';
        case 2:
            cut_c[reordering.depth-2] = 'y';
        case 1:
            cut_c[reordering.depth-1] = 'x';
        }
        for (int v=0; v<reordering.depth-3; v++)
            cut_c[v] ='x'+2-(reordering.depth-1-v)%3;
        std::string cut_string(cut_c);
        LOG::cout<<std::endl<<cut_string<<std::endl;
        reordering.Set_Cut_Array(cut_string);
        delete[] cut_c;
    }
    reordering.Create_Linear_Indices();
    LOG::cout<<"###################### "<<reordering.size<<" ###################"<<std::endl;

    // MATRICES
    T* K_matrix_raw = new T [Nx*Ny*Nz*3*3*3*d*d];
    matrix_type K_matrix = reinterpret_cast<matrix_type>(*K_matrix_raw);
    Populate_Stiffness_Matrix(K_matrix,rK);
    LOG::cout<<"Populate_Stiffness_Matrix"<<std::endl;

    T* x = new T [Nx*Ny*Nz*d];
    T* x_ref = new T [Nx*Ny*Nz*d];
    T* b = new T [Nx*Ny*Nz*d];

    using variable_type = T (&)[Nx][Ny][Nz][d];
    variable_type x_ref_variable = reinterpret_cast<variable_type>(*x_ref);
    variable_type y_variable = reinterpret_cast<variable_type>(*x);
    variable_type x_variable = reinterpret_cast<variable_type>(*b);

    Fill_Variable<T, d, Nx, Ny, Nz>(x_ref_variable, 1, rX);
    Multiply<T, d, Nx, Ny, Nz>(K_matrix, x_ref_variable, y_variable);   // Kx=y

    index_type n;
    index_type *ia_B, *ja_B, *perm;
    T* a_B;
    {
        using namespace Macroblock_Utilities;
        CSC_Topology topology{size.Product()};
        Reordering_To_BSC_Topology(reordering.linear_indices, reordering.coordinates, topology);
        n = topology.n;
        {
            LOG::SCOPE scope("generating scalar topology");
            ia_B = new index_type[n*d+1];
            ja_B = new index_type[(topology.offsets[n]-n)*d*d+n*d*(d+1)/2];
            Block_To_Scalar_Topology(ia_B, ja_B, &topology.offsets[0], &topology.row[0], n);
        }
        // std::cout<<"number of nodes = "<<topology.n<<std::endl;
        // std::cout<<"number of nonzeros = "<<topology.offsets[topology.n]<<std::endl;
    }
    {

        // fill in data
        LOG::SCOPE scope("fill data");
        a_B = new T[(ia_B[n]-n)*d*d+n*d*(d+1)/2];
        n=n*d;
        for (index_type j=0; j<n; j++)
            for (index_type ii=ia_B[j]; ii<ia_B[j+1]; ii++) {
                index_type i=ja_B[ii];
                const T_INDEX& index = reordering.coordinates(j/d+1);
                const T_INDEX& neighbor_index = reordering.coordinates(i/d+1);
                a_B[ii] = K_matrix[index(1)-1][index(2)-1][index(3)-1][neighbor_index(1)-index(1)+1][neighbor_index(2)-index(2)+1][neighbor_index(3)-index(3)+1][j%d][i%d];
            }
        //std::cout<<"n = "<<n<<std::endl;
    }
    perm = new index_type[n];
#if 1
    // reorder the result
    for (int i=0; i<Nx; i++)
    for (int j=0; j<Ny; j++)
    for (int k=0; k<Nz; k++) {
        const T_INDEX index = VECTOR<int,3>(i+1,j+1,k+1);
        b[(reordering.linear_indices(index)-1)*3]=y_variable[index(1)-1][index(2)-1][index(3)-1][0];
        b[(reordering.linear_indices(index)-1)*3+1]=y_variable[index(1)-1][index(2)-1][index(3)-1][1];
        b[(reordering.linear_indices(index)-1)*3+2] = y_variable[index(1)-1][index(2)-1][index(3)-1][2];
    }

    // for (int i=0; i<n; i++)
    //     std::cout<<x[i]<<" "<<b[i]<<std::endl;
    // ############################################################################
    //   PARDISO SOLVE
    // ############################################################################

/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
        index_type mtype = 2;       /* Real symmetric positive definite matrix */
        //index_type mtype = -2;       /* Real symmetric (maybe indefinite) matrix */
        // index_type mtype=11; /* Real not symmetric matrix */


        index_type nrhs = 1;     /* Number of right hand sides. */
        /* Internal solver memory pointer pt, */
        /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
        /* or void *pt[64] should be OK on both architectures */
        void *pt[64];
        /* Pardiso control parameters. */
        index_type iparm[64];
        index_type maxfct, mnum, phase, error, msglvl;
        /* Auxiliary variables. */
        index_type i;
        T ddum;        /* Scalar dummy */
        //index_type idum;         /* Integer dummy. */

    for ( i = 0; i < 64; i++ )
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 3;         /* Use omp */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 2;         /* No user fill-in reducing permutation, return permutaion in idum (could change this, but will ignore iparm[2]) */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not inplace */
    iparm[7] = 0;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 0;        /* not needed for mtype = 2 */
    iparm[10] = 0;        /* not needed for mtype = 2 */
    iparm[11] = 0;        /* Not transpose or conjugate A */
    iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[13] = 0;        /* Output: Number of perturbed pivots (not important for mtype = 2) */
    iparm[14] = 0;        /* Output: Peak memory on symbolic fact */
    iparm[15] = 0;        /* Output: Permanent memory on symbolic fact and solve */
    iparm[16] = 0;        /* Output: Total mem consumption */
    // total peak mem = max(iparm[14],iparm[15])+iparm[16]
    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    iparm[23] = 1;        /* Two-level factorization*/
    iparm[24] = 2;        /*parallel solve*/
    iparm[26] = 1;        /* Check matrix for errors */
    iparm[27] = std::is_same<T,float>::value;
                          /* float or double precision */
    iparm[34] = 1;        /*1 for 0-based index*/
    //iparm[36]=1; /* Use Blocked CSR format for input, need to have iparm[12] = 0 */
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Number of matrix */
    msglvl = 1;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */

    /* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */
    phase = 11;
    PARDISO_64 (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a_B, ia_B, ja_B, perm, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during symbolic factorization: %d", error);
        exit (1);
    }
    printf ("\nReordering completed ... ");
    printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
    printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
    printf ("\nPeak Mem on sym fact = %d kB", iparm[14]);
    printf ("\nPerm Mem on sym fact = %d kB", iparm[15]);
    printf ("\nInternal Mem on num fact and solve = %d kB", iparm[16]);
    printf ("\nTotal peak Mem on sym fact = %f GB",(double)(max(iparm[14],iparm[15])+iparm[16])/1024/1024);

    #if 1
    // total peak mem = max(iparm[14],iparm[15])+iparm[16]
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
    phase = 22;
    PARDISO_64 (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a_B, ia_B, ja_B, perm, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during numerical factorization: %d", error);
        exit (2);
    }
    printf ("\nFactorization completed ... ");
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
#if 1
    phase = 33;
    iparm[7] = 0;         /* Max numbers of iterative refinement steps. */

    PARDISO_64 (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a_B, ia_B, ja_B, perm, &nrhs, iparm, &msglvl, b, x, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during solution: %d", error);
        exit (3);
    }
    printf ("\nSolve completed ... ");
#endif
#endif
#if 1
    // reorder the result
    for (int i=0; i<Nx; i++)
    for (int j=0; j<Ny; j++)
    for (int k=0; k<Nz; k++) {
        const T_INDEX index = VECTOR<int,3>(i+1,j+1,k+1);
        x_variable[index(1)-1][index(2)-1][index(3)-1][0]=x[(reordering.linear_indices(index)-1)*3];
        x_variable[index(1)-1][index(2)-1][index(3)-1][1]=x[(reordering.linear_indices(index)-1)*3+1];
        x_variable[index(1)-1][index(2)-1][index(3)-1][2]=x[(reordering.linear_indices(index)-1)*3+2];
    }
    std::cout<<std::endl;
    Compare_To_Reference(x_ref_variable,x_variable);
    std::cout<<"L_inf = "<<L_Inf_Distance(x_ref_variable, x_variable)<<std::endl;
#endif
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
    phase = -1;           /* Release internal memory. */
    PARDISO_64 (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia_B, ja_B, perm, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);

    // ############################################################################
    // epilogue
    // ############################################################################
    delete[] perm;
    delete[] a_B;
    delete[] ia_B;
    delete[] ja_B;
#endif
    delete[] K_matrix_raw;
    delete[] x;
    delete[] x_ref;
    delete[] b;
    LOG::Finish_Logging();
    }
}


namespace PhysBAM {
    template <class T, int d, size_t Nx, size_t Ny, size_t Nz>
    void Multiply(const T (&K_matrix)[Nx][Ny][Nz][3][3][3][d][d], const T (&x)[Nx][Ny][Nz][d], T (&y)[Nx][Ny][Nz][d]) {
        // y=K*x
        for(int i=0;i<Nx;i++)
            for(int j=0;j<Ny;j++)
                for(int k=0;k<Nz;k++)
                    for(int v=0;v<d;v++)
                    {
                        y[i][j][k][v] = T();
                        for(int di=-1;di<=1;di++)
                            for(int dj=-1;dj<=1;dj++)
                                for(int dk=-1;dk<=1;dk++)
                                    if( (i+di >= 0) && (i+di < Nx) && (j+dj >= 0) && (j+dj < Ny) && (k+dk >= 0) && (k+dk < Nz) )
                                        for(int w=0;w<d;w++)
                                            y[i][j][k][v] += K_matrix[i][j][k][di+1][dj+1][dk+1][v][w] * x[i+di][j+dj][k+dk][w];}
    }

}

template<class T>
void printCSR(const MKL_INT* ia, const MKL_INT* ja,const T* a, MKL_INT n, bool onebased)
{
    const int width = 3;
    for (int i=0; i<n; i++)
    {
        MKL_INT prev=(onebased?1:0);
        for (MKL_INT jjj=ia[i]; jjj<ia[i+1]; jjj++)
        {
            const MKL_INT jj=jjj-(onebased?1:0);
            const MKL_INT j=ja[jj];
            for (int w=0; w<j-prev-1; w++)
                for (int v=0; v<width; v++)
                    std::cout<<" ";
            std::cout<<std::setw(width)<<a[jj];
            prev=j;
        }
        std::cout<<std::endl;
    }
}

template<class T>
T Compute_L_Inf_Norm(const T* data, int N)
{
    T result = abs(data[0]);
    for (int i=0; i<N; i++)
        if (abs(data[i])>result)
            result = abs(data[i]);
    return result;
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

template <class IndexType>
void Block_To_Scalar_Topology(IndexType* offsets, IndexType* rows, const int* offsets_B, const int* rows_B, const int n) {
    constexpr int d = 3;
    //std::cout<<"n = "<<n*d<<std::endl;
    // offsets = new IndexType[n*d+1];
    // rows = new IndexType[(offsets_B[n]-n)*d*d+n*d*(d+1)/2];
    IndexType index = 0;
    for (int j=0; j<n; j++)
        for (int v=0; v<d; v++){
            offsets[j*d+v] = index;
            // diagonal terms
            for (int w=v; w<d; w++)
                rows[index++] = (IndexType)j*d+w;

            for (int ii=offsets_B[j]+1; ii<offsets_B[j+1]; ii++)
                for (int w=0; w<d; w++)
                    rows[index++] = (IndexType)rows_B[ii]*d+w;
        }
    offsets[n*d] = index;
     // std::cout<<"offsets_B[n] = "<<offsets_B[n]<<" n = "<<n<<std::endl;
     // std::cout<<"index = "<<offsets[n*d]<<" (offsets_B[n]-n)*d*d+n*d*(d+1)/2 = "<<(offsets_B[n]-n)*d*d+n*d*(d+1)/2<<std::endl;

}
