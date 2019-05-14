//#####################################################################
// Copyright 2012-2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
using namespace PhysBAM;

#include "BlockedMatrixNXN_Cholesky.h"
#include "Matrix4.h"
#include "Matrix4_AVX.h"
#include "Matrix4_AVX512.h"

template<class T>
inline void Block_Cholesky(BlockedMatrixNXN<T> &A)
{
    for(int j=0,Ilj=0;j<A.n;j++){

        Matrix4<float,float> Dj;
        Dj.Load(A.D(j));
        Dj.Invert();
        Dj.Store(A.D(j));

        for(int l=j+1;l<A.n;l++,Ilj++){

            Matrix4<float,float> Llj;
            Llj.Load(A.L(Ilj));
            {
                Matrix4<float,float> Dl;
                Dl.Load(A.D(l));
                Dl-=Llj.Conjugate(Dj);
                Dl.Store(A.D(l));

            }
            Llj*=Dj;
            Llj.Store(A.L(Ilj));

            for(int k=l+1,Ikj=Ilj+1,Ikl=((2*A.n-l-1)*l)/2;k<A.n;k++,Ikj++,Ikl++){
                Matrix4<float,float> Lkj,Lkl;
                Lkj.Load(A.L(Ikj));
                Lkl.Load(A.L(Ikl));
                Lkl-=Lkj.Times_Transpose(Llj);
                Lkl.Store(A.L(Ikl));
            }            
        }
    }
}

template<class T>
inline void Block_Cholesky_From_Scalar(BlockedMatrixNXN<T> &A)
{
    for(int j=0,Ilj=0;j<A.n;j++){

        Matrix4<float,float> Dj;
        Dj.Load(A.D(j));

        for(int l=j+1;l<A.n;l++,Ilj++){
            Matrix4<float,float> Llj;
            Llj.Load(A.L(Ilj));
            Dj.Backward_Substitution_On_Rows(Llj);
            Llj.Store(A.L(Ilj));
        }

        Dj.Conjugate_Diagonal_With_Lower_Triangular_Inverse_Transpose();
        Dj.Store(A.D(j));

    }
}

template<class Tw,class T>
void Block_Forward_Substitution(const BlockedMatrixNXN<T> &A,T* x,const bool use_diagonal)
{
    for(int j=0,index=0;j<A.n;j++){

        Vector4<T,Tw,true> xj;
        xj.Broadcast(&x[4*j]);

        for(int i=j+1;i<A.n;i++,index++){
            Matrix4<T,Tw> L;
            L.Load_Prefetch<128>(A.L(index));
            auto xi=L*xj;
            xi.Accumulate_Negative(&x[4*i]);
        }

        if(use_diagonal){
            Matrix4<T,Tw> D;
            D.Load(A.D(j));
            auto xj_scaled=D*xj;
            xj_scaled.Store(&x[4*j]);
        }

    }
}

template<class Tw,class T>
void Block_Backward_Substitution(const BlockedMatrixNXN<T> &A,T* x,const bool use_diagonal)
{
    for(int i=A.n-1,index=(A.n*(A.n-1))/2-1;i>=0;i--){

        Vector4<T,Tw,true> xi;
        if(use_diagonal){
            Vector4<T,Tw,false> xiT;
            xiT.Broadcast(&x[4*i]);            
            Matrix4<T,Tw> D;
            D.Load(A.D(i));
            xi=D*xiT;
        }
        else
            xi.Load(&x[4*i]);

        for(int j=A.n-1;j>i;j--,index--){
            Matrix4<T,Tw> L;
            L.Load_Prefetch<-256>(A.L(index));
            Vector4<T,Tw,false> xj;
            xj.Broadcast(&x[4*j]);
            xi-=L*xj;
        }
        
        xi.Store(&x[4*i]);

    }
}

template void Block_Cholesky(BlockedMatrixNXN<float> &A);
template void Block_Cholesky_From_Scalar(BlockedMatrixNXN<float> &A);
template void Block_Forward_Substitution<float,float>(const BlockedMatrixNXN<float> &A,float* x,const bool use_diagonal);
template void Block_Backward_Substitution<float,float>(const BlockedMatrixNXN<float> &A,float* x,const bool use_diagonal);
template void Block_Forward_Substitution<__m256,float>(const BlockedMatrixNXN<float> &A,float* x,const bool use_diagonal);
template void Block_Backward_Substitution<__m256,float>(const BlockedMatrixNXN<float> &A,float* x,const bool use_diagonal);
template void Block_Forward_Substitution<__m512,float>(const BlockedMatrixNXN<float> &A,float* x,const bool use_diagonal);
template void Block_Backward_Substitution<__m512,float>(const BlockedMatrixNXN<float> &A,float* x,const bool use_diagonal);

