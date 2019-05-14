//#####################################################################
// Copyright 2012-2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef SUBROUTINE_Sparse_Blocked_Forwards_Substitution
#include "Common/KernelCommon.h"
#include "BlockedCSCSymmetricMatrix3.h"
#include "DiagonalMatrix3.h"
#include "SymmetricMatrix3.h"
#include "Matrix3.h"
#include "UnitriangularMatrix3.h"
#include <cassert>
#else
namespace {
#endif

template<class Tw,class T_DATA,class I_DATA>
#ifdef SUBROUTINE_Sparse_Blocked_Forwards_Substitution
inline
#endif
void Sparse_Blocked_Forwards_Substitution(const BlockedCSCSymmetricMatrix3<T_DATA> &U, typename std::remove_extent<T_DATA>::type *b_raw, bool use_diagonal)
{
    static_assert( std::rank<T_DATA>::value == std::rank<I_DATA>::value, "Error: Base datatypes mismatch in width." );
    static_assert( std::rank<T_DATA>::value == 1 || std::rank<T_DATA>::value == 0 , "Error: Base datatypes must be a rank 1 array or a scalar." );

    typedef Number<Tw> Tn;
    typedef typename std::remove_extent<T_DATA>::type T_Base;
    const int T_SIZE = std::rank<T_DATA>::value == 0 ? 1 : std::extent<T_DATA>::value;
    typedef T_Base (*VEC3_PTR_TYPE)[3][T_SIZE];    

    VEC3_PTR_TYPE b=reinterpret_cast<VEC3_PTR_TYPE>(b_raw);

    for(int j = 0;j<U.n;j++){
        Vector3<Tn> x_j;
        x_j.Load_Aligned(b[j]);

        UnitriangularMatrix3<Tn> D_tri;
        D_tri.Load_Aligned(U.D(j));
        D_tri.In_Place_Forward_Substitution(x_j);       

        for(int ii = U.Offsets(j);ii<U.Offsets(j+1);ii++){
            int i=U.Row(ii);
            Matrix3<Tn> U_ij;
            U_ij.Load_Aligned(U.L(ii));
            Vector3<Tn> x_i;
            x_i.Load_Aligned(b[i]);
            x_i=x_i-U_ij*x_j;
            x_i.Store(b[i]);}       

        if(use_diagonal){
            DiagonalMatrix3<Tn> D_diag;
            D_diag.Load_Aligned(U.D(j));
            x_j=D_diag*x_j;}
        
        x_j.Store(b[j]);
    }
}

#ifndef SUBROUTINE_Sparse_Blocked_Forwards_Substitution
#define INSTANCE_KERNEL_Sparse_Blocked_Forwards_Substitution(WIDTH) const BlockedCSCSymmetricMatrix3< WIDETYPE(float,WIDTH) > &U, float *b, bool use_diagonal
INSTANCE_KERNEL(Sparse_Blocked_Forwards_Substitution);
#undef INSTANCE_KERNEL_Sparse_Blocked_Forwards_Substitution

#define INSTANCE_KERNEL_Sparse_Blocked_Forwards_Substitution(WIDTH) const BlockedCSCSymmetricMatrix3< WIDETYPE(double,WIDTH) > &U, double *b, bool use_diagonal
INSTANCE_KERNEL_DOUBLE(Sparse_Blocked_Forwards_Substitution);
#undef INSTANCE_KERNEL_Sparse_Blocked_Forwards_Substitution

#else
}
#endif

////////////////


