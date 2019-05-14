//#####################################################################
// Copyright 2012-2018, Qisi Wang, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef SUBROUTINE_Sparse_Blocked_CSC_Multiply
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
// only works for double and float for now
    template<class Arch_Type>
    struct Wide_One_Policy;

    template<>
    struct Wide_One_Policy<double> {
        constexpr static double value[1] = {1.};
    };

    template<>
    struct Wide_One_Policy<float> {
        constexpr static float value[1] = {1.};
    };


template<class Tw,class T_DATA,class I_DATA>
#ifdef SUBROUTINE_Sparse_Blocked_CSC_Multiply
inline
#endif

void Sparse_Blocked_CSC_Multiply(BlockedCSCSymmetricMatrix3<T_DATA> &A, const int end)
{
    typedef Number<Tw> Tn;
    Tn ONE;

    ONE.Load_Aligned(Wide_One_Policy<Tw>::value);

    for(int j=end-1;j>=0;j--){
        DiagonalMatrix3<Tn> D;
        D.Load_Aligned(A.D(j));
        D.Inverse();

        for (int ll=A.Offsets(j); ll<A.Offsets(j+1); ll++)
        {
            int l=A.Row(ll);
            Matrix3<Tn> M;
            M.Load_Aligned(A.L(ll));
            {
                SymmetricMatrix3<Tn> S;
                S.Load_Aligned(A.D(l));
                S+=M.Conjugate(D);
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
                    P+=Q;
                    P.Store(A.L(rr));
                }
            }
        }

        {
            SymmetricMatrix3<Tn> L;
            L.Load_Aligned(A.D(j));
            Matrix3<Tn> L_full(ONE,  L.x21, L.x31,
                               Tn(), ONE,   L.x32,
                               Tn(), Tn(),  ONE   );
            for(int ii=A.Offsets(j);ii<A.Offsets(j+1);ii++){
                Matrix3<Tn> M;
                M.Load_Aligned(A.L(ii));
                Matrix3<Tn> Q=M.Times_Transpose(L_full);
                Q.Store(A.L(ii));
            }
            // compute A(j,j)
            L = L_full.Conjugate(D);
            L.Store(A.D(j));
        }
    }
}


#ifndef SUBROUTINE_Sparse_Blocked_CSC_Multiply
    template void Sparse_Blocked_CSC_Multiply<float, float, int>(BlockedCSCSymmetricMatrix3<float> &A, const int end);
    template void Sparse_Blocked_CSC_Multiply<double, double, int>(BlockedCSCSymmetricMatrix3<double> &A, const int end);
    #if 0
#define INSTANCE_KERNEL_Sparse_Blocked_CSC_Multiply(WIDTH) BlockedCSCSymmetricMatrix3< WIDETYPE(float,WIDTH) > &A, const int end
INSTANCE_KERNEL(Sparse_Blocked_CSC_Multiply);
#undef INSTANCE_KERNEL_Sparse_Blocked_CSC_Multiply
    #endif
#else
}
#endif
