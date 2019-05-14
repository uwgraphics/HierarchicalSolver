#ifndef SUBROUTINE_Sparse_Cholesky_Partial
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
#ifdef SUBROUTINE_Sparse_Cholesky_Partial
inline
#endif
    void Sparse_Cholesky_Partial(BlockedCSCSymmetricMatrix3<T_DATA> &A, const int stop, const int start)
{
    // perform cholesky on a slice of columns
    typedef Number<Tw> Tn;
    assert(stop<=A.n);
    for(int j=start;j<stop;j++){

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

#ifndef SUBROUTINE_Sparse_Cholesky_Partial

#define INSTANCE_KERNEL_Sparse_Cholesky_Partial(WIDTH) BlockedCSCSymmetricMatrix3< WIDETYPE(float,WIDTH) > &A, const int stop, const int
INSTANCE_KERNEL(Sparse_Cholesky_Partial);
#undef INSTANCE_KERNEL_Sparse_Cholesky_Partial

#define INSTANCE_KERNEL_Sparse_Cholesky_Partial(WIDTH) BlockedCSCSymmetricMatrix3< WIDETYPE(double,WIDTH) > &A, const int stop, const int
INSTANCE_KERNEL_DOUBLE(Sparse_Cholesky_Partial);
#undef INSTANCE_KERNEL_Sparse_Cholesky_Partial

#else
}
#endif
