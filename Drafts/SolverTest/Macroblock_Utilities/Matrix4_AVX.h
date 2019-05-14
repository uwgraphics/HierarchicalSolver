//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Matrix4_AVX_h__
#define __Matrix4_AVX_h__

#include <immintrin.h>
//using namespace PhysBAM;

template<class T,class Tw,bool broadcast_in_lane> struct Vector4;
template<class T,class Tw> struct Matrix4;

template<> struct Vector4<float,__m256,true>
{
    __m256 data0;
    __m256 data1;

    Vector4& operator-=(const Vector4& v)
    {
        data0=_mm256_sub_ps(data0,v.data0);
        data1=_mm256_sub_ps(data1,v.data1);
    }

    void Load(const float* x)
    {
#if 1
        __m256 tmp=_mm256_broadcast_ps((__m128 const*)x);
        data0=_mm256_setzero_ps();
        data1=_mm256_setzero_ps();
        data0=_mm256_blend_ps(data0,tmp,0x21);
        data1=_mm256_blend_ps(data1,tmp,0x84);
#else
        alignas(64) float data[16];
        for(int k=0;k<16;k++) data[k]=0.f;
        for(int k=0;k<16;k+=5) data[k]=x[k/4];
        data0=_mm256_load_ps(&data[0]);
        data1=_mm256_load_ps(&data[8]);
#endif
    }

    void Broadcast(const float* x)
    {
#if 1
        data0=_mm256_set_m128(_mm_broadcast_ss(x+1),_mm_broadcast_ss(x));
        data1=_mm256_set_m128(_mm_broadcast_ss(x+3),_mm_broadcast_ss(x+2));
#else
        alignas(64) float data[16];
        for(int k=0;k<16;k++) data[k]=x[k/4];
        data0=_mm256_load_ps(&data[0]);
        data1=_mm256_load_ps(&data[8]);
#endif
    }

    void Store(float* x) const
    {
#if 1
        __m256 tmp0=_mm256_permute_ps(data0,0x4e);
        tmp0=_mm256_add_ps(tmp0,data0);
        __m256 tmp1=_mm256_permute_ps(data1,0x4e);
        tmp1=_mm256_add_ps(tmp1,data1);        
        tmp0=_mm256_blend_ps(tmp0,tmp1,0xcc);
        tmp1=_mm256_permute_ps(tmp0,0xb1);
        tmp0=_mm256_add_ps(tmp0,tmp1);
        tmp1=_mm256_permute2f128_ps(tmp0,tmp0,0x81);
        tmp0=_mm256_blend_ps(tmp0,tmp1,0xfa);
        _mm_store_ps(x,_mm256_castps256_ps128(tmp0));
#else
        alignas(64) float data[16];
        _mm256_store_ps(&data[0],data0);
        _mm256_store_ps(&data[8],data1);
        for(int k=0;k<4;k++) x[k]=0.f;
        for(int k=0;k<16;k++) x[k/4]+=data[k];
#endif
    }
};

template<> struct Vector4<float,__m256,false>
{
    __m256 data0;
    __m256 data1;

    void Broadcast(const float* x)
    {
#if 1
        data0=data1=_mm256_broadcast_ps((__m128 const*)x);
#else
        alignas(64) float data[16];
        for(int k=0;k<16;k++) data[k]=x[k%4];        
        data0=_mm256_load_ps(&data[0]);
        data1=_mm256_load_ps(&data[8]);        
#endif
    }

    void Store(float* x) const
    {
#if 1
        __m256 tmp0=_mm256_add_ps(data0,data1);
        __m256 tmp1=_mm256_permute2f128_ps(tmp0,tmp0,0x81);
        tmp0=_mm256_add_ps(tmp0,tmp1);
        _mm_store_ps(x,_mm256_castps256_ps128(tmp0));
#else
        alignas(64) float data[16];
        _mm256_store_ps(&data[0],data0);
        _mm256_store_ps(&data[8],data1);
        for(int k=0;k<4;k++) x[k]=0.f;
        for(int k=0;k<16;k++) x[k%4]+=data[k];
#endif
    }

    void Accumulate_Negative(float* x) const
    {
#if 1
        __m256 tmp0=_mm256_add_ps(data0,data1);
        __m256 tmp1=_mm256_permute2f128_ps(tmp0,tmp0,0x81);
        tmp0=_mm256_add_ps(tmp0,tmp1);
        __m128 dst=_mm_load_ps(x);
        dst=_mm_sub_ps(dst,_mm256_castps256_ps128(tmp0));
        _mm_store_ps(x,dst);
#else
        alignas(64) float data[16];
        _mm256_store_ps(&data[0],data0);
        _mm256_store_ps(&data[8],data1);
        for(int k=0;k<16;k++) x[k%4]-=data[k];
#endif
    }
};

template<> struct Matrix4<float,__m256>
{
    typedef float T;
    typedef __m256 Tw;

    __m256 data0;
    __m256 data1;

    void Load(const float (&A)[16])
    {
        data0=_mm256_load_ps(&A[0]);
        data1=_mm256_load_ps(&A[8]);
    }

    template<int offset>
    void Load_Prefetch(const float (&A)[16])
    {
        data0=_mm256_load_ps(&A[0]);
        data1=_mm256_load_ps(&A[8]);
    }

    template<bool broadcast_in_lane>
    Vector4<T,Tw,!broadcast_in_lane> operator*(const Vector4<T,Tw,broadcast_in_lane>& x) const
    {
        Vector4<T,Tw,!broadcast_in_lane> y;  
        y.data0=_mm256_mul_ps(data0,x.data0);
        y.data1=_mm256_mul_ps(data1,x.data1);
        return y;
    }

};

#endif
