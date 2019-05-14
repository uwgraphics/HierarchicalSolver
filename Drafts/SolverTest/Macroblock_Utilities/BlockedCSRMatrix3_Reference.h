//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef  __BlockedCSRMatrix3_Reference_h__
#define  __BlockedCSRMatrix3_Reference_h__


//################################################
// Class BlockedCSRMatrix3_Reference
//################################################

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>


using namespace PhysBAM;

template<class T>
struct BlockedCSRMatrix3_Reference
// Sparse symmetric matrix in blocked Compressed Sparse Column format, with Matrix3x3 blocks
{
    int n, m;        // Matrix dimension in blocks (row,col)
    HASHTABLE<VECTOR<int, 2>, MATRIX<T, 3> > L_hash;

    BlockedCSRMatrix3_Reference(const int n_input, const int m_input)
    :n(n_input), m(m_input)
    {
    }
    
    MATRIX<T,3,3>& L( const int l, const int j )  {
        PHYSBAM_ASSERT( 0<=l && l<n &&  0<=j && j<m );
        return L_hash.Get( VECTOR<int,2>(l,j));
    }

    MATRIX<T,3,3>& L_Or_Insert( const int l, const int j ) {
                PHYSBAM_ASSERT( 0<=l && l<n &&  0<=j && j<m );
        return L_hash.Get_Or_Insert( VECTOR<int,2>(l,j) );
    }
    
    const MATRIX<T,3,3>& L( const int l, const int j ) const {
                PHYSBAM_ASSERT( 0<=l && l<n &&  0<=j && j<m );
        return L_hash.Get( VECTOR<int,2>(l,j));
    }

    const MATRIX<T,3,3>& L_Or_Insert( const int l, const int j ) const {
                PHYSBAM_ASSERT( 0<=l && l<n &&  0<=j && j<m );
        return L_hash.Get_Or_Insert( VECTOR<int,2>(l,j) );
    }
    

    void RandomInitialize(int seed){
        RANDOM_NUMBERS<T> random_numbers;random_numbers.Set_Seed(seed);
        L_hash.Clean_Memory();
        
        for( int i = 0; i < n; i++ )
            for( int j = 0; j < m; j++ ){
                bool active = random_numbers.Get_Uniform_Integer(0, 2);
                if( active ){
                    random_numbers.Fill_Uniform( L_Or_Insert( i, j ), 0, 1 );
                    L_Or_Insert( i, j )(1,1) *= 5;
                    L_Or_Insert( i, j )(2,2) *= 5;
                    L_Or_Insert( i, j )(3,3) *= 5;
                }
            }
        std::cout << "Initialized matrix with " << L_hash.Size() << " entries." << std::endl;
    }
            
};

//#####################################################################
#endif
