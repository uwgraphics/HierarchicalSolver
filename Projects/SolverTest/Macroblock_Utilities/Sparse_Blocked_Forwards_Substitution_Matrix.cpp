//#####################################################################
// Copyright 2012-2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef SUBROUTINE_Sparse_Blocked_Forwards_Substitution_Matrix
#include "Common/KernelCommon.h"
#include "BlockedCSCSymmetricMatrix3.h"
#include "BlockedCSRMatrix3.h"
#include "DiagonalMatrix3.h"
#include "SymmetricMatrix3.h"
#include "Matrix3.h"
#include "UnitriangularMatrix3.h"
#include <cassert>
#else
namespace {
#endif

template<class Tw,class T_DATA,class I_DATA>
#ifdef SUBROUTINE_Sparse_Blocked_Forwards_Substitution_Matrix
inline
#endif
void Sparse_Blocked_Forwards_Substitution_Matrix(const BlockedCSCSymmetricMatrix3<T_DATA> &U,const BlockedCSRMatrix3<T_DATA> &A,BlockedCSRMatrix3<T_DATA> &B)
{
    static_assert( std::rank<T_DATA>::value == std::rank<I_DATA>::value, "Error: Base datatypes mismatch in width." );
    static_assert( std::rank<T_DATA>::value == 1 || std::rank<T_DATA>::value == 0 , "Error: Base datatypes must be a rank 1 array or a scalar." );

    typedef Number<Tw> Tn;
    typedef typename std::remove_extent<T_DATA>::type T_Base;
    const int T_SIZE = std::rank<T_DATA>::value == 0 ? 1 : std::extent<T_DATA>::value;
    typedef T_Base (*MAT3_PTR_TYPE)[9][T_SIZE];

    // MAT3_PTR_TYPE B=reinterpret_cast<MAT3_PTR_TYPE>(B_raw);
    
    assert(A.rows==B.rows);
    assert(A.columns==B.columns);
    assert(A.columns==U.n);

    for(int i=0;i<B.rows;i++){
        
        // Copy A->B, inserting zeros when necessary
        for(int jj=B.offsets[i],kk=A.offsets[i];jj<B.offsets[i+1];jj++)
            if( (kk==A.offsets[i+1]) || (B.column[jj]<A.column[kk]) ){
                Matrix3<Tn> Zero;
                Zero.Store(B.E(jj));}
            else{
                assert(B.column[jj]==A.column[kk]);
                Matrix3<Tn> Xij;
                Xij.Load_Aligned(A.E(kk++));
                Xij.Store(B.E(jj));}

        // Perform sparse forward substitution, one block row at a time

        for(int jj=B.offsets[i];jj<B.offsets[i+1];jj++){

            int j=B.column[jj];
            Matrix3<Tn> X_ij;
            X_ij.Load_Aligned(B.E(jj));            

            UnitriangularMatrix3<Tn> D_tri;
            D_tri.Load_Aligned(U.D(j));
            D_tri.In_Place_Forward_Substitution(X_ij);

            for(int kk=U.Offsets(j),ll=jj+1;kk<U.Offsets(j+1);kk++){
                int k=U.Row(kk);
                while(B.column[ll]<k){ ll++; assert(ll<B.offsets[i+1]); }
                assert(B.column[ll]==k);
                Matrix3<Tn> X_ik;
                X_ik.Load_Aligned(B.E(ll));
                Matrix3<Tn> U_kj;
                U_kj.Load_Aligned(U.L(kk));
                X_ik-=U_kj*X_ij;
                X_ik.Store(B.E(ll));
            }
        
            X_ij.Store(B.E(jj));            
        }

    }
}

#ifndef SUBROUTINE_Sparse_Blocked_Forwards_Substitution_Matrix
#define INSTANCE_KERNEL_Sparse_Blocked_Forwards_Substitution_Matrix(WIDTH) const BlockedCSCSymmetricMatrix3< WIDETYPE(float,WIDTH) > &U, const BlockedCSRMatrix3< WIDETYPE(float,WIDTH) > &A, BlockedCSRMatrix3< WIDETYPE(float,WIDTH) > &B
INSTANCE_KERNEL(Sparse_Blocked_Forwards_Substitution_Matrix);
#undef INSTANCE_KERNEL_Sparse_Blocked_Forwards_Substitution_Matrix
#else
}
#endif

////////////////
