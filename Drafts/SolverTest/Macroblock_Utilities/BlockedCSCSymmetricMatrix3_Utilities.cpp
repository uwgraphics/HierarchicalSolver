//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "BlockedCSCSymmetricMatrix3.h"
#include "BlockedCSCSymmetricMatrix3_Reference.h"
#include "BlockedCSCSymmetricMatrix3_Utilities.h"

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_NXN.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <iomanip>
#include <array>
#include <cstring>
#include <algorithm>

using namespace PhysBAM;

namespace{
//################################################
// Function Randomize
//################################################
template<class T>
void Randomize(SYMMETRIC_MATRIX_NXN<T>& S,const int seed)
{RANDOM_NUMBERS<T> random_numbers;random_numbers.Set_Seed(seed);
for(int k=0;k<S.size;k++) S.x[k]=random_numbers.Get_Number();}

template<class T>
void Make_Positive_Definite(SYMMETRIC_MATRIX_NXN<T>& S)
{
    ARRAY<VECTOR<int,2> > givens_pairs;
    ARRAY<VECTOR<T,2> > givens_coefficients;
    S.Jacobi_Solve_Eigenproblem(givens_pairs,givens_coefficients);
    for(int j=1;j<S.n;j++) for(int i=j+1;i<=S.n;i++) S(i,j)=T();   // Clear lower-triangular section
    for(int i=1;i<=S.n;i++) S(i,i)=clamp_min(abs(S(i,i)),(T)1e-4); // Flip sign of negative eigenvalues (and clamp to small number)
    // Reconstruct by conjugating with givens matrices
    for(int k=givens_pairs.m;k>=1;k--){ int i,j; givens_pairs(k).Get(i,j); T c,s; givens_coefficients(k).Get(c,s); S.Givens_Conjugate(i,j,c,-s); }
}
}

//################################################
// Function Print_Matrix
//################################################
template<class T>
void Print_Matrix(const BlockedCSCSymmetricMatrix3_Reference<T>& A,std::ostream& out)
{
    ARRAY<std::string> rows(3*A.n);
    for(int i=0;i<A.n;i++){
        for(int j=0;j<i;j++)
            if(A.L_hash.Contains(VECTOR<int,2>(i,j))){
                const MATRIX<T,3>& L=A.L_hash.Get(VECTOR<int,2>(i,j));
                for(int k=1;k<=3;k++) for(int l=1;l<=3;l++)
                    rows(3*i+k)+=STRING_UTILITIES::string_sprintf("%12.6g% ",L(k,l));}
            else for(int k=1;k<=3;k++) rows(3*i+k)+="                           ";
        const SYMMETRIC_MATRIX<T,3>& D=A.D(i);
        for(int k=1;k<=3;k++) for(int l=1;l<=k;l++)
            rows(3*i+k)+=STRING_UTILITIES::string_sprintf("%12.6g% ",D(k,l));}
    for(int i=1;i<=3*A.n;i++)
        out<<rows(i)<<std::endl;
}
template<class T>
void Print_Matrix(const BlockedCSCSymmetricMatrix3<T>& A,std::ostream& out)
{
    ARRAY<std::string> rows(3*A.n);

    for(int j=0;j<A.n;j++){
        const SYMMETRIC_MATRIX<T,3>& D=reinterpret_cast<const SYMMETRIC_MATRIX<T,3>&>(A.D(j));
        for(int k=1;k<=3;k++) for(int l=1;l<=k;l++)
            rows(3*j+k)+=STRING_UTILITIES::string_sprintf("%12.6g% ",D(k,l));
        int i=j+1;
        for(int ii=A.offsets[j];ii<A.offsets[j+1];ii++,i++){
            for(;i<A.row[ii];i++)
                for(int k=1;k<=3;k++) rows(3*i+k)+="                           ";
            const MATRIX<T,3>& L=reinterpret_cast<const MATRIX<T,3>&>(A.L(ii));
            for(int k=1;k<=3;k++) for(int l=1;l<=3;l++)
                rows(3*i+k)+=STRING_UTILITIES::string_sprintf("%12.6g% ",L(k,l));}
        for(;i<A.n;i++)
            for(int k=1;k<=3;k++) rows(3*i+k)+="                           ";
    }

    for(int i=1;i<=3*A.n;i++)
        out<<rows(i)<<std::endl;
}
//################################################
// Function InitializeBlockedLapacian
//################################################
template< class T >
void InitializeBlockedLapacian( BlockedCSCSymmetricMatrix3_Reference<T> &A )
{
    typedef VECTOR<int, 3> T_INDEX;
    
    HASHTABLE<VECTOR<int,3>,int> map;
    int flat_index = 1;
    for( RANGE_ITERATOR<3> iter( RANGE<T_INDEX>( T_INDEX(0,0,0), T_INDEX(1,1,1) ) ); iter.Valid(); iter.Next() )
        map.Insert(iter.Index(), flat_index++);
    
    SYMMETRIC_MATRIX<T, 3> diag_elem( 4, 0, 0, 4, 0, 4 );
    MATRIX<T,3> offdiag_elem( -1, 0, 0, 0, -1, 0, 0, 0, -1 );

    for( RANGE_ITERATOR<3> X( RANGE<T_INDEX>( T_INDEX(0,0,0), T_INDEX(1,1,1) ) ); X.Valid(); X.Next() )
        for( RANGE_ITERATOR<3> Y( RANGE<T_INDEX>( T_INDEX(0,0,0), T_INDEX(1,1,1) ) ); Y.Valid(); Y.Next() )
            {        
                if( X.Index() == Y.Index() )
                    A.D(map.Get( X.Index())) = diag_elem;
                                
                if( (X.Index()  - Y.Index()).Magnitude_Squared() == 1 &&
                     map.Get( Y.Index() ) < map.Get( X.Index() ) )
                    A.L_Or_Insert( map.Get( X.Index() ), map.Get( Y.Index() ) ) = offdiag_elem;
            }
    
}
//################################################
// Function InitializeBlockedDense
//################################################
template<class T>
void InitializeBlockedDense( BlockedCSCSymmetricMatrix3_Reference<T> &A, int seed, int size )
{
    typedef VECTOR<int, 3> T_INDEX;

    SYMMETRIC_MATRIX<T, 3> diag_elem( 1, 0, 0, 1, 0, 1 );
    MATRIX<T,3> offdiag_elem( 0, 0, 0, 0, 0, 0, 0, 0, 0 );

    for( int i = 0; i < size; i++ ){
        for( int j = 0; j < size; j++ ){
            if( i == j )
                A.D( i ) = diag_elem;               
            
            if( i != j && j < i )
                A.L_Or_Insert( i, j ) = offdiag_elem;
        }
    }              
}

//################################################
// Function InitializeBlockedIdentity
//################################################
template< class T >
void InitializeBlockedIdentity( BlockedCSCSymmetricMatrix3_Reference<T> &A )
{
    typedef VECTOR<int, 3> T_INDEX;
    
    SYMMETRIC_MATRIX<T, 3> diag_elem( 1, 0, 0, 1, 0, 1 );
    MATRIX<T,3> offdiag_elem( 0, 0, 0, 0, 0, 0, 0, 0, 0 );

    for( int i = 0; i < A.n; i++){
        A.D(i) = diag_elem;
    }
    
    for( HASHTABLE_ITERATOR<VECTOR<int, 2>, MATRIX<T, 3> > iter( A.L_hash); iter.Valid(); iter.Next() ){
        iter.Data() = offdiag_elem;
    }
}

//################################################
// Function LoadFromPhysbam
//################################################
template<class T, class T_DATA>
void LoadFromPhysbam( BlockedCSCSymmetricMatrix3<T_DATA> &A, const BlockedCSCSymmetricMatrix3_Reference<T> &S )
{
    assert( A.n == S.n );
    typedef const T (&SYM_TYPE)[6];
    typedef const T (&MAT_TYPE)[9];

    typedef T_DATA (*Lower_Triangular_Matrix3_ptr)[6];
    typedef T_DATA (&Lower_Triangular_Matrix3_type)[6];
    typedef const T_DATA (&Lower_Triangular_Matrix3_const_type)[6];

    typedef T_DATA (*Matrix3_ptr)[9];
    typedef T_DATA (&Matrix3_type)[9];
    typedef const T_DATA (&Matrix3_const_type)[9];

    const int T_DATA_FCOUNT = (sizeof(T_DATA) / sizeof( std::remove_all_extents<T_DATA>::type));

    for( int Bj = 0; Bj < A.n; Bj++)
        for( int Bi = Bj; Bi < A.n; Bi++)
            if( Bj == Bi ){
                const SYMMETRIC_MATRIX<T,3>& diag = S.D(Bi);
                SYM_TYPE diag_casted = reinterpret_cast<SYM_TYPE>( diag );
                T expanded_diag_casted[6*T_DATA_FCOUNT];
                for( int i = 0; i < 6; i++ )
                    for( int j = 0; j < T_DATA_FCOUNT; j++)
                        expanded_diag_casted[ i * T_DATA_FCOUNT + j] = diag_casted[i];                           
                memcpy( A.D(Bi), expanded_diag_casted, 6 * sizeof(T_DATA) );
            }
            else{
                int index = -1;
                for( int r = A.offsets[Bj]; r < A.offsets[Bj+1]; r++)
                    if( A.row[r] == Bi )
                        index = r;
                if( index == -1)
                    continue; // This block was not allocated, skip to preserve sparcity

                const MATRIX<T,3>& offdiag = S.L(Bi,Bj);
                MAT_TYPE offdiag_casted = reinterpret_cast<MAT_TYPE>( offdiag );
                T expanded_offdiag_casted[9*T_DATA_FCOUNT];
                for( int i = 0; i < 9; i++ )
                    for( int j = 0; j < T_DATA_FCOUNT; j++)
                        expanded_offdiag_casted[ i * T_DATA_FCOUNT + j] = offdiag_casted[i];
                memcpy( A.L(index), expanded_offdiag_casted, 9 * sizeof(T_DATA) );
            }   
    
}

template<class T, class T_DATA>
void LoadFromPhysbam( BlockedCSCSymmetricMatrix3<T_DATA> &A, const SYMMETRIC_MATRIX_NXN<T>& S )
{
    assert( 3*A.n == S.n );

    for( int Bj = 0; Bj < A.n; Bj++){
        for( int Bi = Bj; Bi < A.n; Bi++){
            // Do the diagonal copy
            if( Bj == Bi ){
                auto& A_diag_block = A.D(Bj);
                
                int idx = 0;
                for( int j = 0; j<3; j++)
                    for( int i = j; i<3; i++){
                        A_diag_block[idx] = S(Bi*3+i+1, Bj*3+j+1);
                        idx++;
                    }         
            }
            // Do the lower copy
            else{
                int index = -1;
                for( int r = A.offsets[Bj]; r < A.offsets[Bj+1]; r++)
                    if( A.row[r] == Bi )
                        index = r;
                if( index == -1)
                    continue; // This block was not allocated, skip to preserve sparcity

                auto& A_lower_block = A.L(index);

                int idx = 0;
                for( int j = 0; j<3; j++)
                    for( int i = 0; i<3; i++){
                        A_lower_block[idx] = S(Bi*3+i+1, Bj*3+j+1);
                        idx++;
                    }                  
            }
        }
    }    

    //SYMMETRIC_MATRIX_NXN<T> S_debug(S.n);
    //StoreToPhysbam(S_debug,A);
    //PHYSBAM_ASSERT((S-S_debug).Max_Abs()==T());
}
//################################################
// Function StoreFromPhysbam
//################################################
template<class T, class T_DATA>
void StoreToPhysbam( BlockedCSCSymmetricMatrix3_Reference<T> &S, const BlockedCSCSymmetricMatrix3<T_DATA> &A ){
    assert( A.n == S.n );
    typedef T (&SYM_TYPE)[6];
    typedef T (&MAT_TYPE)[9];

    typedef T_DATA (*Lower_Triangular_Matrix3_ptr)[6];
    typedef T_DATA (&Lower_Triangular_Matrix3_type)[6];
    typedef const T_DATA (&Lower_Triangular_Matrix3_const_type)[6];

    typedef T_DATA (*Matrix3_ptr)[9];
    typedef T_DATA (&Matrix3_type)[9];
    typedef const T_DATA (&Matrix3_const_type)[9];

    const int T_DATA_FCOUNT = (sizeof(T_DATA) / sizeof( std::remove_all_extents<T_DATA>::type));   

    for( int Bj = 0; Bj < A.n; Bj++)
        for( int Bi = Bj; Bi < A.n; Bi++)
            if( Bj == Bi ){
                SYMMETRIC_MATRIX<T,3>& diag = S.D(Bi);
                SYM_TYPE diag_casted = reinterpret_cast<SYM_TYPE>( diag );
                T expanded_diag_casted[6*T_DATA_FCOUNT];
                memcpy( expanded_diag_casted, A.D(Bi), 6 * sizeof(T_DATA) );
                for( int i = 0; i < 6; i++ )
                    diag_casted[i] = expanded_diag_casted[ i * T_DATA_FCOUNT ];                           
            }
            else{
                int index = -1;
                for( int r = A.offsets[Bj]; r < A.offsets[Bj+1]; r++)
                    if( A.row[r] == Bi )
                        index = r;
                if( index == -1)
                    continue; // This block was not allocated, skip to preserve sparcity

                MATRIX<T,3>& offdiag = S.L_Or_Insert(Bi,Bj);
                MAT_TYPE offdiag_casted = reinterpret_cast<MAT_TYPE>( offdiag );
                T expanded_offdiag_casted[9*T_DATA_FCOUNT];
                memcpy( expanded_offdiag_casted, A.L(index), 9 * sizeof(T_DATA) );
                for( int i = 0; i < 9; i++ )
                    offdiag_casted[i] = expanded_offdiag_casted[ i * T_DATA_FCOUNT];

            }       
}


template<class T, class T_DATA>
void StoreToPhysbam( SYMMETRIC_MATRIX_NXN<T>& S, const BlockedCSCSymmetricMatrix3<T_DATA> &A )
{
    assert( 3*A.n == S.n );

    for( int Bj = 0; Bj < A.n; Bj++){
        for( int Bi = Bj; Bi < A.n; Bi++){
            // Do the diagonal copy
            if( Bj == Bi ){
                auto& A_diag_block = A.D(Bj);
                
                int idx = 0;
                for( int j = 0; j<3; j++)
                    for( int i = j; i<3; i++){
                        S(Bi*3+i+1, Bj*3+j+1) = A_diag_block[idx];                       
                        idx++;
                    }               
            }
            // Do the lower copy
            else{
                int index = -1;
                for( int r = A.offsets[Bj]; r < A.offsets[Bj+1]; r++)
                    if( A.row[r] == Bi )
                        index = r;
                if( index == -1 ){
                    int idx = 0;
                    for( int j = 0; j<3; j++)
                        for( int i = 0; i<3; i++){
                            S(Bi*3+i+1, Bj*3+j+1) = 0;
                            idx++;
                        }       
                }else{
                    auto& A_lower_block = A.L(index);
                    
                    int idx = 0;
                    for( int j = 0; j<3; j++)
                        for( int i = 0; i<3; i++){
                            S(Bi*3+i+1, Bj*3+j+1) = A_lower_block[idx]; 
                            idx++;
                        }             
                }     
            }
        }
    }    
}

template<class T>
void StoreToPhysbam ( MATRIX_MXN<T>& S, const BlockedCSCSymmetricMatrix3_Reference<T> &A )
{
    assert( 3*A.n <= S.n );
    for( int Bj = 0; Bj < A.n; Bj++){
        auto& A_diag_block = A.D(Bj);
        for( int j = 0; j<3; j++)
            for( int i = j; i<3; i++)
                S(Bj*3+i+1, Bj*3+j+1) = A_diag_block(i+1,j+1);                               
    }
    for( HASHTABLE_ITERATOR<VECTOR<int, 2>, const MATRIX<T, 3> > iter( A.L_hash); iter.Valid(); iter.Next() ){
        auto& A_lower_block = iter.Data();
        for( int j = 0; j<3; j++)
            for( int i = 0; i<3; i++)
                S(iter.Key().x*3+i+1, iter.Key().y*3+j+1) = A_lower_block(i+1,j+1);        
    }
}
//################################################
// Function Dense_Cholesky_Allocate
//################################################
template<class T_DATA>
void Dense_Cholesky_Allocate(const int block_dim,void *&D_raw,void *&L_raw,int *&offsets,int *&row)
{
    std::cout<<"In Dense_Cholesky_Allocate()"<<std::endl;
    std::cout<<"block_dim="<<block_dim<<std::endl;
    std::cout<<"sizeof(T_DATA)="<<sizeof(T_DATA)<<std::endl;
    posix_memalign(&D_raw,64,block_dim*6*sizeof(T_DATA));
    posix_memalign(&L_raw,64,((block_dim*(block_dim-1))/2)*9*sizeof(T_DATA));
    offsets=new int[block_dim+1];
    row=new int[(block_dim*(block_dim-1))/2];
    int index=0;
    for(int j=0;j<block_dim;j++){
        offsets[j]=index;
        for(int i=j+1;i<block_dim;i++)
            row[index++]=i;}
    if(index!=(block_dim*(block_dim-1))/2) abort();
    offsets[block_dim]=index;
}
//################################################
// Function Sparse_Cholesky_Allocate_From_Reference
//################################################
template<class T, class T_DATA>
void Sparse_Cholesky_Allocate_From_Reference(void *&D_raw,void *&L_raw,int *&offsets,int *&row,
                                             const BlockedCSCSymmetricMatrix3_Reference<T> &A)
{
    //std::cout<<"In Sparse_Cholesky_Allocate_From_Reference()"<<std::endl;
    //std::cout<<"block_dim="<<A.n<<std::endl;
    //std::cout<<"sizeof(T_DATA)="<<sizeof(T_DATA)<<std::endl;
    posix_memalign(&D_raw,64,A.n*6*sizeof(T_DATA));
    posix_memalign(&L_raw,64,(A.L_hash.Size())*9*sizeof(T_DATA));
    
    struct VectorCompare {
        bool operator() (VECTOR<int,2> i, VECTOR<int,2> j) {
            if(i.y < j.y) return true;
            else if( i.y == j.y )
                return i.x < j.x;
            return false;
        }
    } vcompare;

    struct myclass {           // function object type:
        void operator() (int i) {std::cout << ' ' << i;}
        void operator() ( const VECTOR<int,2>& v ) {std::cout << ' ' << v << std::endl;}
    } myobject;

    offsets=new int[A.n+1];
    row=new int[A.L_hash.Size()];

    ARRAY< VECTOR<int, 2> > index_list;
    index_list.Preallocate(A.L_hash.Size());
    for( HASHTABLE_ITERATOR< VECTOR<int,2>, const MATRIX<T,3> > hiter( A.L_hash ); hiter.Valid(); hiter.Next() ){
        index_list.Append( hiter.Key() );
    }
    std::sort(index_list.begin(), index_list.end(), vcompare);

    //std::cout << "Entries: "<< std::endl;
    //std::for_each (index_list.begin(), index_list.end(), myobject); std::cout << std::endl;

    int flat_index = 1;
    int sparse_index = 0;
    for( int col = 0; col < A.n; col ++ ){
        bool col_done = false;
        offsets[col] = sparse_index;
        //std::for_each (offsets, offsets+(A.n+1), myobject); std::cout << std::endl;


        if( flat_index > index_list.m ){
            offsets[col] = index_list.m;
            continue;
        }           
        do {
            const VECTOR<int,2> current_index = index_list(flat_index);
            if( current_index.y != col ){
                col_done = true;
                continue;
            }
            else{                
                //std::cout << "Marking location ("<< current_index.x << ", " << col << ") as filled with index " << sparse_index << std::endl;

                row[sparse_index] = current_index.x;
                sparse_index++;
                flat_index++;
            }
        } while( !col_done && flat_index <= index_list.m );
    }
    offsets[A.n] = index_list.m;           
 
    //std::cout << "Final offsets: " << std::endl;
    //std::for_each (offsets, offsets+(A.n+1), myobject); std::cout << std::endl;
       
}
//################################################
// Function Sparse_Cholesky_Initialize
//################################################
template<> void Sparse_Cholesky_Initialize(BlockedCSCSymmetricMatrix3<float> &A,const int seed)
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true,true,true,true,false,false);

    typedef float T;
    SYMMETRIC_MATRIX_NXN<T> S(3*A.n);
    Randomize(S,seed);
    Make_Positive_Definite(S);
    LoadFromPhysbam(A, S);

#if 0
    // Test to confirm copy worked.
    SYMMETRIC_MATRIX_NXN<T> S_copy(n);
    StoreToPhysbam(S_copy, A);
    Print_Formatted(S, std::cout);
    Print_Formatted(S_copy, std::cout);
    Print_Formatted(S-S_copy, std::cout);

    // The following actually goes in "compare"
    // MATRIX_MXN<T> U(S.n),L;
    // for(int i=1;i<=S.n;i++) for(int j=1;j<=S.n;j++) U(i,j)=S(i,j);
    // U.In_Place_LU_Factorization(L);
    // for(int i=1;i<=S.n;i++){ for(int j=1;j<i;j++) S(i,j)=L(i,j); S(i,i)=T(1.f)/U(i,i);
#endif
}
//#####################################################################
template void InitializeBlockedLapacian<double>(BlockedCSCSymmetricMatrix3_Reference<double> &A);
template void InitializeBlockedIdentity<double>(BlockedCSCSymmetricMatrix3_Reference<double> &A);
template void InitializeBlockedDense<double>( BlockedCSCSymmetricMatrix3_Reference<double> &A, int seed, int size );
template void InitializeBlockedLapacian<float>(BlockedCSCSymmetricMatrix3_Reference<float> &A);
template void InitializeBlockedIdentity<float>(BlockedCSCSymmetricMatrix3_Reference<float> &A);
template void InitializeBlockedDense<float>( BlockedCSCSymmetricMatrix3_Reference<float> &A, int seed, int size );


template void Print_Matrix(const BlockedCSCSymmetricMatrix3_Reference<double>& A,std::ostream& out);
template void Print_Matrix(const BlockedCSCSymmetricMatrix3<double>& A,std::ostream& out);

template void Print_Matrix(const BlockedCSCSymmetricMatrix3_Reference<float>& A,std::ostream& out);
template void Print_Matrix(const BlockedCSCSymmetricMatrix3<float>& A,std::ostream& out);


template void LoadFromPhysbam<double,double>(BlockedCSCSymmetricMatrix3<double> &A, const BlockedCSCSymmetricMatrix3_Reference<double> &S);
template void LoadFromPhysbam<float,float>(BlockedCSCSymmetricMatrix3<float> &A, const BlockedCSCSymmetricMatrix3_Reference<float> &S);
template void LoadFromPhysbam<float,float[4]>(BlockedCSCSymmetricMatrix3<float[4]> &A, const BlockedCSCSymmetricMatrix3_Reference<float> &S);
template void LoadFromPhysbam<float,float[8]>(BlockedCSCSymmetricMatrix3<float[8]> &A, const BlockedCSCSymmetricMatrix3_Reference<float> &S);
template void LoadFromPhysbam<float,float[16]>(BlockedCSCSymmetricMatrix3<float[16]> &A, const BlockedCSCSymmetricMatrix3_Reference<float> &S);

template void StoreToPhysbam<double,double>( BlockedCSCSymmetricMatrix3_Reference<double> &S, const BlockedCSCSymmetricMatrix3<double> &A );
template void StoreToPhysbam<float,float>( BlockedCSCSymmetricMatrix3_Reference<float> &S, const BlockedCSCSymmetricMatrix3<float> &A );
template void StoreToPhysbam<float,float[4]>( BlockedCSCSymmetricMatrix3_Reference<float> &S, const BlockedCSCSymmetricMatrix3<float[4]> &A );
template void StoreToPhysbam<float,float[8]>( BlockedCSCSymmetricMatrix3_Reference<float> &S, const BlockedCSCSymmetricMatrix3<float[8]> &A );
template void StoreToPhysbam<float,float[16]>( BlockedCSCSymmetricMatrix3_Reference<float> &S, const BlockedCSCSymmetricMatrix3<float[16]> &A );

template void Sparse_Cholesky_Allocate_From_Reference<double, double> (void *&D_raw,void *&L_raw,int *&offsets,int *&row, const BlockedCSCSymmetricMatrix3_Reference<double> &A);
template void Sparse_Cholesky_Allocate_From_Reference<float, float> (void *&D_raw,void *&L_raw,int *&offsets,int *&row, const BlockedCSCSymmetricMatrix3_Reference<float> &A);
template void Sparse_Cholesky_Allocate_From_Reference<float, float[4]> (void *&D_raw,void *&L_raw,int *&offsets,int *&row, const BlockedCSCSymmetricMatrix3_Reference<float> &A);
template void Sparse_Cholesky_Allocate_From_Reference<float, float[8]> (void *&D_raw,void *&L_raw,int *&offsets,int *&row, const BlockedCSCSymmetricMatrix3_Reference<float> &A);
template void Sparse_Cholesky_Allocate_From_Reference<float, float[16]> (void *&D_raw,void *&L_raw,int *&offsets,int *&row, const BlockedCSCSymmetricMatrix3_Reference<float> &A);
template void Sparse_Cholesky_Allocate_From_Reference<VECTOR<float,16>, float[16]> (void *&D_raw,void *&L_raw,int *&offsets,int *&row, const BlockedCSCSymmetricMatrix3_Reference<VECTOR<float,16> > &A );

template void Dense_Cholesky_Allocate<double>(const int n,void *&D_raw,void *&L_raw,int *&offsets,int *&row);
template void Dense_Cholesky_Allocate<float>(const int n,void *&D_raw,void *&L_raw,int *&offsets,int *&row);
template void Dense_Cholesky_Allocate<float[8]>(const int n,void *&D_raw,void *&L_raw,int *&offsets,int *&row);
template void LoadFromPhysbam<float,float>( BlockedCSCSymmetricMatrix3<float> &A, const SYMMETRIC_MATRIX_NXN<float>& S );
template void StoreToPhysbam<float>( SYMMETRIC_MATRIX_NXN<float>& S, const BlockedCSCSymmetricMatrix3<float> &A );
template void StoreToPhysbam<float>( MATRIX_MXN<float>& S, const BlockedCSCSymmetricMatrix3_Reference<float> &A );
