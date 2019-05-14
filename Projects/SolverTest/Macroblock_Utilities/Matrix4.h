//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Matrix4_h__
#define __Matrix4_h__

#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <cstring>

using namespace PhysBAM;

template<class T,class Tw,bool broadcast_in_lane> struct Vector4;
template<class T,class Tw> struct Matrix4;

template<> struct Vector4<float,float,true>
{
    typedef float T;
    typedef float Tw;

    T data[16];

    Vector4& operator-=(const Vector4& v)
    {for(int k=0;k<16;k++) data[k]-=v.data[k];return *this;}

    void Load(const float* x)
    {for(int k=0;k<16;k++) data[k]=0.f;
    for(int k=0;k<16;k+=5) data[k]=x[k/4];}

    void Broadcast(const float* x)
    {for(int k=0;k<16;k++) data[k]=x[k/4];}

    void Store(float* x) const
    {for(int k=0;k<4;k++) x[k]=0.f;
    for(int k=0;k<16;k++) x[k/4]+=data[k];}

    void Accumulate(float* x) const
    {for(int k=0;k<16;k++) x[k/4]+=data[k];}

    void Accumulate_Negative(float* x) const
    {for(int k=0;k<16;k++) x[k/4]-=data[k];}
};

template<> struct Vector4<float,float,false>
{
    typedef float T;
    typedef float Tw;

    T data[16];

    Vector4& operator-=(const Vector4& v)
    {for(int k=0;k<16;k++) data[k]-=v.data[k];return *this;}

    void Broadcast(const float* x)
    {for(int k=0;k<16;k++) data[k]=x[k%4];}

    void Store(float* x) const
    {for(int k=0;k<4;k++) x[k]=0.f;
    for(int k=0;k<16;k++) x[k%4]+=data[k];}

    void Accumulate(float* x) const
    {for(int k=0;k<16;k++) x[k%4]+=data[k];}

    void Accumulate_Negative(float* x) const
    {for(int k=0;k<16;k++) x[k%4]-=data[k];}
};

template<> struct Matrix4<float,float>
{
    typedef float T;
    typedef float Tw;

    T data[16];

    void Load(const float (&A)[16])
    {for(int k=0;k<16;k++) data[k]=A[k];}

    template<int offset>
    void Load_Prefetch(const float (&A)[16])
    {for(int k=0;k<16;k++) data[k]=A[k];}

    void Store(float (&A)[16]) const
    {for(int k=0;k<16;k++) A[k]=data[k];}

    void Invert(){
        MATRIX<T,4> A;
        std::memcpy(A.x,data,16*sizeof(T));
        A.Invert();
        std::memcpy(data,A.x,16*sizeof(T));
    }

    template<bool broadcast_in_lane>
    Vector4<T,Tw,!broadcast_in_lane> operator*(const Vector4<T,Tw,broadcast_in_lane>& x) const
    {
        Vector4<T,Tw,!broadcast_in_lane> y;
        for(int k=0;k<16;k++)
            y.data[k]=data[k]*x.data[k];
        return y;
    }

    Matrix4& operator*=(const Matrix4& M)
    {
        MATRIX<T,4> A,B;
        std::memcpy(A.x,data,16*sizeof(T));
        std::memcpy(B.x,M.data,16*sizeof(T));
        A*=B;
        std::memcpy(data,A.x,16*sizeof(T));
        return *this;
    }

    Matrix4& operator-=(const Matrix4& M)
    {
        MATRIX<T,4> A,B;
        std::memcpy(A.x,data,16*sizeof(T));
        std::memcpy(B.x,M.data,16*sizeof(T));
        A-=B;
        std::memcpy(data,A.x,16*sizeof(T));
        return *this;
    }

    Matrix4 Times_Transpose(const Matrix4& M) const
    {
        MATRIX<T,4> A,B,C;
        std::memcpy(A.x,data,16*sizeof(T));
        std::memcpy(B.x,M.data,16*sizeof(T));
        C=A.Times_Transpose(B);
        Matrix4 result;
        std::memcpy(result.data,C.x,16*sizeof(T));
        return result;
    }

    Matrix4 Conjugate(const Matrix4& D) const
    {
        MATRIX<T,4> A,B;
        std::memcpy(A.x,data,16*sizeof(T));
        std::memcpy(B.x,D.data,16*sizeof(T));
        B=B.Times_Transpose(A);
        A*=B;
        Matrix4 result;
        std::memcpy(result.data,A.x,16*sizeof(T));
        return result;
    }

    void Backward_Substitution_On_Rows(Matrix4& M) const
    {
        MATRIX<T,4> L,A;
        std::memcpy(L.x,data,16*sizeof(T));
        for(int i=1;i<=4;i++){
            L(i,i)=(T)1.f;
            for(int j=i+1;j<=4;j++)
                L(i,j)=T();}
        std::memcpy(A.x,M.data,16*sizeof(T));
        A*=L.Inverse();
        std::memcpy(M.data,A.x,16*sizeof(T));
    }

    void Conjugate_Diagonal_With_Lower_Triangular_Inverse_Transpose()
    {
        MATRIX<T,4> L,D;
        std::memcpy(L.x,data,16*sizeof(T));
        std::memset(D.x,0,16*sizeof(T));
        for(int i=1;i<=4;i++){
            D(i,i)=L(i,i);
            L(i,i)=(T)1.f;
            for(int j=i+1;j<=4;j++)
                L(i,j)=T();}
        L=L.Inverse();
        L=L.Transpose_Times(D)*L;
        std::memcpy(data,L.x,16*sizeof(T));
    }

};

#endif
