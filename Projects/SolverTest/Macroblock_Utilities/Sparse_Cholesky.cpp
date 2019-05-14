//#####################################################################
// Copyright 2012-2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef SUBROUTINE_Sparse_Cholesky
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
#ifdef SUBROUTINE_Sparse_Cholesky
inline
#endif
void Sparse_Cholesky(BlockedCSCSymmetricMatrix3<T_DATA> &A)
{
    typedef Number<Tw> Tn;

    for(int j=0;j<A.n;j++){

        // Perform LDLt (with inverted diagonal) transformation on pivot

        DiagonalMatrix3<Tn> D;
        {
            UnitriangularMatrix3<Tn> L;
            {
                SymmetricMatrix3<Tn> S;
                S.Load_Aligned(A.D(j));
                S.LDLt(L,D);
            }
            L.Store(A.D(j));
            D.Store(A.D(j));

            // Initially right-scale entire column with L^{-T}

            for(int ii=A.Offsets(j);ii<A.Offsets(j+1);ii++){
                Matrix3<Tn> M;
                M.Load_Aligned(A.L(ii));
                L.In_Place_Forward_Substitution_On_Transpose(M);
                M.Store(A.L(ii));
            }
        }

        for(int ll=A.Offsets(j);ll<A.Offsets(j+1);ll++){

            int l=A.Row(ll);
            Matrix3<Tn> M;
            M.Load_Aligned(A.L(ll));
            {
                SymmetricMatrix3<Tn> S;
                S.Load_Aligned(A.D(l));
                S-=M.Conjugate(D);
                S.Store(A.D(l));
            }
            D.In_Place_Column_Scale(M);
            M.Store(A.L(ll));

            for(int rr=A.Offsets(l),kk=ll+1;kk<A.Offsets(j+1);kk++){

                for(;A.row[rr]!=A.row[kk];rr++) assert(A.row[rr]<A.row[kk] && rr+1<A.Offsets(l+1));

                {
                    Matrix3<Tn> N;
                    N.Load_Aligned(A.L(kk));
                    Matrix3<Tn> Q=N.Times_Transpose(M),P;
                    P.Load_Aligned(A.L(rr));
                    P-=Q;
                    P.Store(A.L(rr));
                }
            }

        }

    }

}


#ifndef SUBROUTINE_Sparse_Cholesky
#define INSTANCE_KERNEL_Sparse_Cholesky(WIDTH) BlockedCSCSymmetricMatrix3< WIDETYPE(float,WIDTH) > &A
INSTANCE_KERNEL(Sparse_Cholesky);
#undef INSTANCE_KERNEL_Sparse_Cholesky

#define INSTANCE_KERNEL_Sparse_Cholesky(WIDTH) BlockedCSCSymmetricMatrix3< WIDETYPE(double,WIDTH) > &A
INSTANCE_KERNEL_DOUBLE(Sparse_Cholesky);
#undef INSTANCE_KERNEL_Sparse_Cholesky
#else
}
#endif
