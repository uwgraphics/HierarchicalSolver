//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "BlockedCSRMatrix3.h"
#include "BlockedCSRMatrix3_Reference.h"
#include "BlockedCSRMatrix3_Utilities.h"

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
//################################################
// Function Randomize
//################################################
template<class T>
void Randomize(SYMMETRIC_MATRIX_NXN<T>& S,const int seed)
{RANDOM_NUMBERS<T> random_numbers;random_numbers.Set_Seed(seed);
for(int k=0;k<S.size;k++) S.x[k]=random_numbers.Get_Number();}

}


//################################################
// Function Sparse_Cholesky_Allocate_From_Reference
//################################################
template<class T, class T_DATA>
void CSR_Allocate_From_Reference(void *&E_raw,int *&offsets,int *&column,
                                 const BlockedCSRMatrix3_Reference<T> &A)
{
    //std::cout<<"In Sparse_Cholesky_Allocate_From_Reference()"<<std::endl;
    //std::cout<<"block_dim="<<A.n<<std::endl;
    //std::cout<<"sizeof(T_DATA)="<<sizeof(T_DATA)<<std::endl;
    posix_memalign(&E_raw,64,(A.L_hash.Size())*9*sizeof(T_DATA));
    
    struct VectorCompare {
       bool operator() (VECTOR<int,2> i, VECTOR<int,2> j) {
            if(i.x < j.x) return true;
            else if( i.x == j.x )
                return i.y < j.y;
            return false;
        }
    } vcompare;

    struct myclass {           // function object type:
        void operator() (int i) {std::cout << ' ' << i;}
        void operator() ( const VECTOR<int,2>& v ) {std::cout << ' ' << v << std::endl;}
    } myobject;

    offsets=new int[A.n+1];
    column=new int[A.L_hash.Size()];

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
    for( int row = 0; row < A.n; row ++ ){
        bool row_done = false;
        offsets[row] = sparse_index;
        //std::for_each (offsets, offsets+(A.n+1), myobject); std::cout << std::endl;


        if( flat_index > index_list.m ){
            offsets[row] = index_list.m;
            continue;
        }           
        do {
            const VECTOR<int,2> current_index = index_list(flat_index);
            if( current_index.x != row ){
                row_done = true;
                continue;
            }
            else{                
                //std::cout << "Marking location ("<< row << ", " << current_index.y << ") as filled with index " << sparse_index << std::endl;

                column[sparse_index] = current_index.y;
                sparse_index++;
                flat_index++;
            }
        } while( !row_done && flat_index <= index_list.m );
    }
    offsets[A.n] = index_list.m;           
 
    //std::cout << "Final offsets: " << std::endl;
    //std::for_each (offsets, offsets+(A.n+1), myobject); std::cout << std::endl;
       
}



//################################################
// Function Load_From_Reference
//################################################
template<class T, class T_DATA>
void Load_From_Reference(BlockedCSRMatrix3<T_DATA> &S, 
                         const BlockedCSRMatrix3_Reference<T> &A){

    typedef const float (&SYM_TYPE)[6];   
    typedef const float (&MAT_TYPE)[9];   
 
    typedef T_DATA (*Lower_Triangular_Matrix3_ptr)[6];
    typedef T_DATA (&Lower_Triangular_Matrix3_type)[6];
    typedef const T_DATA (&Lower_Triangular_Matrix3_const_type)[6];

    typedef T_DATA (*Matrix3_ptr)[9];
    typedef T_DATA (&Matrix3_type)[9];
    typedef const T_DATA (&Matrix3_const_type)[9];
    const int T_DATA_FCOUNT = (sizeof(T_DATA) / sizeof( std::remove_all_extents<T_DATA>::type));   

    for( int Bi = 0; Bi < A.n; Bi++){
        for( int Bj = 0; Bj < A.m; Bj++){
            int index = -1;
            for( int r = S.offsets[Bi]; r < S.offsets[Bi+1]; r++)
                if( S.column[r] == Bj )
                    index = r;
            if( index == -1)
                continue; // This block was not allocated, skip to preserve sparcity
            
            const MATRIX<T,3>& offdiag = A.L(Bi,Bj);
            MAT_TYPE offdiag_casted = reinterpret_cast<MAT_TYPE>( offdiag );
            float expanded_offdiag_casted[9*T_DATA_FCOUNT];
            for( int i = 0; i < 9; i++ )
                for( int j = 0; j < T_DATA_FCOUNT; j++)
                    expanded_offdiag_casted[ i * T_DATA_FCOUNT + j] = offdiag_casted[i];
            memcpy( S.E(index), expanded_offdiag_casted, 9 * sizeof(T_DATA) );
                   
        }
    }  


}


//################################################
// Function StoreFromPhysbam
//################################################
template<class T, class T_DATA>
void StoreToPhysbam( MATRIX_MXN<T>& S, const BlockedCSRMatrix3<T_DATA> &A )
{
    for( int Bi = 0; Bi < A.rows; Bi++){
        for( int Bj = 0; Bj < A.columns; Bj++){
            int index = -1;
            for( int r = A.offsets[Bi]; r < A.offsets[Bi+1]; r++)
                if( A.column[r] == Bj )
                    index = r;
            if( index == -1 ){
                int idx = 0;
                for( int j = 0; j<3; j++)
                    for( int i = 0; i<3; i++){
                        S(Bi*3+i+1, Bj*3+j+1) = 0;
                        idx++;
                    }       
            }else{
                auto& A_lower_block = A.E(index);
                
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


//#################################################################################

template void CSR_Allocate_From_Reference<float, float> (void *&E_raw,int *&offsets,int *&column, const BlockedCSRMatrix3_Reference<float> &A);
template void CSR_Allocate_From_Reference<float, float[4]> (void *&E_raw,int *&offsets,int *&column, const BlockedCSRMatrix3_Reference<float> &A);
template void CSR_Allocate_From_Reference<float, float[8]> (void *&E_raw,int *&offsets,int *&column, const BlockedCSRMatrix3_Reference<float> &A);
template void CSR_Allocate_From_Reference<float, float[16]> (void *&E_raw,int *&offsets,int *&column, const BlockedCSRMatrix3_Reference<float> &A);
template void CSR_Allocate_From_Reference<VECTOR<float,16>, float[16]> (void *&E_raw,int *&offsets,int *&column, const BlockedCSRMatrix3_Reference<VECTOR<float,16> > &A);

template void Load_From_Reference<float,float>(BlockedCSRMatrix3<float> &A, const BlockedCSRMatrix3_Reference<float> &S);
template void Load_From_Reference<float,float[4]>(BlockedCSRMatrix3<float[4]> &A, const BlockedCSRMatrix3_Reference<float> &S);
template void Load_From_Reference<float,float[8]>(BlockedCSRMatrix3<float[8]> &A, const BlockedCSRMatrix3_Reference<float> &S);
template void Load_From_Reference<float,float[16]>(BlockedCSRMatrix3<float[16]> &A, const BlockedCSRMatrix3_Reference<float> &S);

template void StoreToPhysbam<float>( MATRIX_MXN<float>& S, const BlockedCSRMatrix3<float> &A );
