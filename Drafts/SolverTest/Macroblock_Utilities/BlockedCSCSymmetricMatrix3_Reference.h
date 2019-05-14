//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef  __BlockedCSCSymmetricMatrix3_Reference_h__
#define  __BlockedCSCSymmetricMatrix3_Reference_h__


//################################################
// Class BlockedCSCSymmetricMatrix3_Reference
//################################################

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>


using namespace PhysBAM;

template<class T>
struct BlockedCSCSymmetricMatrix3_Reference
// Sparse symmetric matrix in blocked Compressed Sparse Column format, with Matrix3x3 blocks
{
    int n;        // Matrix dimension in blocks
    ARRAY<SYMMETRIC_MATRIX<T, 3>, VECTOR<int,1> > D;
    HASHTABLE<VECTOR<int, 2>, MATRIX<T, 3> > L_hash;

    BlockedCSCSymmetricMatrix3_Reference(const int n_input)
        :n(n_input)
    {
        D.Resize(0,n-1);
    }

    MATRIX<T,3,3>& L( const int l, const int j )  {
        PHYSBAM_ASSERT( 0<=j && j<l && l<n );
        return L_hash.Get( VECTOR<int,2>(l,j));
    }

    MATRIX<T,3,3>& L_Or_Insert( const int l, const int j ) {
        PHYSBAM_ASSERT( 0<=j && j<l && l<n );
        return L_hash.Get_Or_Insert( VECTOR<int,2>(l,j) );
    }

    const MATRIX<T,3,3>& L( const int l, const int j ) const {
        PHYSBAM_ASSERT( 0<=j && j<l && l<n );
        return L_hash.Get( VECTOR<int,2>(l,j));
    }

    const MATRIX<T,3,3>& L_Or_Insert( const int l, const int j ) const {
        PHYSBAM_ASSERT( 0<=j && j<l && l<n );
        return L_hash.Get_Or_Insert( VECTOR<int,2>(l,j) );
    }

    BlockedCSCSymmetricMatrix3_Reference SubMatrix( const ARRAY<int>& indices ) const {
        BlockedCSCSymmetricMatrix3_Reference result(indices.m);
        for( int i = 1; i <= indices.m; i++ ){
            result.D(i-1) = D(indices(i));
        }
        for( HASHTABLE_ITERATOR< VECTOR<int, 2>, const MATRIX<T, 3> > iter(L_hash); iter.Valid(); iter.Next() ){
            int i = indices.Find(iter.Key().x)-1;
            int j = indices.Find(iter.Key().y)-1;
            if( i >=0 && j >=0 )
                result.L_Or_Insert( i, j ) = iter.Data();
        }
        return result;
    }

    void ExtendSparsity(){
        //std::cout << "Extending sparsity for cholesky operations" << std::endl;
        //std::cout << "Total occupency: " << ( (float)(L_hash.Size() + n) / (float)(n*n)) * 100 << "%, with " <<  (L_hash.Size() + n) << " active entries." << std::endl;
        // Run a stripped down version of sparse cholesky to look for missing non-zeros
        int new_entries = 0;

        for( int j = 0; j < n; j++ ){ // For each column
            //std::cout << "Processing column " << j << std::endl;

            for( int l = j+1; l < n; l++ ){ // Examine each potential row

                if( ! L_hash.Contains( VECTOR<int,2>(l,j) ))
                    continue; // Skip over rows without entries

                //std::cout << "Entry found at (" << l << ", " << j << ")" << std::endl;

                for( int r = l+1, k = l+1; k < n; k++ ){ // Now we look at pairs in the columns of j and l

                    if( ! L_hash.Contains( VECTOR<int,2>(k,j) ))
                        continue; // Skip over rows without entries

                    for(;r!=k;r++){ // Real work here.
                        assert(  r < n ); } // We can't go over the size of the matrix...

                    // But we do need to add entries to column l to match those in column j
                    if( ! L_hash.Contains( VECTOR<int,2>(r,l) )) {
                        L_hash.Insert( VECTOR<int,2>(r,l), MATRIX<T, 3>() );
                        //std::cout << "Adding additional entry: ("<<r<<", " << l << ")" << std::endl;
                        new_entries++;
                    }

                }
            }
        }

        //std::cout << "Sparsity was extended by adding " << new_entries << " additional entries." << std::endl;
        //std::cout << "Total occupency is now: " << ( (float)(L_hash.Size() + n) / (float)(n*n)) * 100 << "%, with " <<  (L_hash.Size() + n) << " active entries." << std::endl;
    }

};

//#####################################################################
#endif
