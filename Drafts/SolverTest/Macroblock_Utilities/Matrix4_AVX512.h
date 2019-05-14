//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Matrix4_AVX512_h__
#define __Matrix4_AVX512_h__

#include <immintrin.h>

//using namespace PhysBAM;

template<class T,class Tw,bool broadcast_in_lane> struct Vector4;
template<class T,class Tw> struct Matrix4;

template<> struct Vector4<float,__m512,true>
{
    __m512 data0;

    Vector4& operator-=(const Vector4& v)
    {
        data0=_mm512_sub_ps(data0,v.data0);
    }

    void Load(const float* x)
    {
#if 1
        data0=_mm512_broadcast_f32x4(_mm_load_ps(x));//x0x1x2x3 x0x1x2x3 x0x1x2x3 x0x1x2x3
        data0=_mm512_mask_blend_ps(_mm512_int2mask(0x8421),_mm512_setzero_ps(),data0);
#else
        alignas(64) float data[16];
        for(int k=0;k<16;k++) data[k]=0.f;
        for(int k=0;k<16;k+=5) data[k]=x[k/4];
        data0=_mm512_load_ps(&data[0]);
#endif
    }

    void Broadcast(const float* x)
    {
#if 1
        /*
          data0=_mm512_mask_blend_ps(_mm512_int2mask(0xFF00),
          _mm512_broadcast_f32x8(_mm256_set_m128(_mm_broadcast_ss(x+1),_mm_broadcast_ss(x  ))),
          _mm512_broadcast_f32x8(_mm256_set_m128(_mm_broadcast_ss(x+3),_mm_broadcast_ss(x+2))));
        */
        data0=_mm512_mask_blend_ps(_mm512_int2mask(0xFF00),
                                   _mm512_mask_blend_ps(_mm512_int2mask(0x00F0),
                                                        _mm512_broadcast_f32x4(_mm_broadcast_ss(x)),
                                                        _mm512_broadcast_f32x4(_mm_broadcast_ss(x+1))),
                                   _mm512_mask_blend_ps(_mm512_int2mask(0xF000),
                                                        _mm512_broadcast_f32x4(_mm_broadcast_ss(x+2)),
                                                        _mm512_broadcast_f32x4(_mm_broadcast_ss(x+3))));

#else
        alignas(64) float data[16];
        for(int k=0;k<16;k++) data[k]=x[k/4];
        data0=_mm512_load_ps(&data[0]);
#endif
    }

    void Store(float* x) const
    {
#if 1
        __m512 tmp0=_mm512_add_ps(data0,_mm512_permute_ps(data0,0xB1));
        __m512 tmp1=_mm512_add_ps(tmp0,_mm512_permute_ps(tmp0,0x0A));//now it holds x0x0x0x0 x1x1x1x1 x2x2x2x2 x3x3x3x3
        tmp0=_mm512_shuffle_f32x4(tmp1,tmp1,0x0E);//x2x2x2x2 x3x3x3x3 .... .... 
        tmp1=_mm512_mask_blend_ps(_mm512_int2mask(0x33),tmp0,tmp1);// x0x0x2x2 x1x1x3x3 .... ....
        tmp0=_mm512_shuffle_f32x4(tmp1,tmp1,0x01);//x1x1x3x3 .... .... .... 
        tmp1=_mm512_mask_blend_ps(_mm512_int2mask(0x05),tmp0,tmp1);// x0x1x2x3 .... .... ....
        _mm_store_ps(x,_mm512_castps512_ps128(tmp1));
#else
        alignas(64) float data[16];
        _mm512_store_ps(&data[0],data0);
        for(int k=0;k<4;k++) x[k]=0.f;
        for(int k=0;k<16;k++) x[k/4]+=data[k];
#endif
    }
};

template<> struct Vector4<float,__m512,false>
{
    __m512 data0;

    void Broadcast(const float* x)
    {
#if 1
        data0=_mm512_broadcast_f32x4(_mm_load_ps(x));
#else
        alignas(64) float data[16];
        for(int k=0;k<16;k++) data[k]=x[k%4];        
        data0=_mm512_load_ps(&data[0]);
#endif
    }

    void Store(float* x) const
    {
#if 1
        __m512 tmp0=_mm512_add_ps(data0,_mm512_shuffle_f32x4(data0,data0,0xB1));
        __m512 tmp1=_mm512_add_ps(tmp0,_mm512_shuffle_f32x4(tmp0,tmp0,0xAA));
        _mm_store_ps(x,_mm512_castps512_ps128(tmp1));
#else
        alignas(64) float data[16];
        _mm512_store_ps(&data[0],data0);
        for(int k=0;k<4;k++) x[k]=0.f;
        for(int k=0;k<16;k++) x[k%4]+=data[k];
#endif
    }

    void Accumulate_Negative(float* x) const
    {
#if 1
        __m512 tmp0=_mm512_add_ps(data0,_mm512_shuffle_f32x4(data0,data0,0xB1));
        __m512 tmp1=_mm512_add_ps(tmp0,_mm512_shuffle_f32x4(tmp0,tmp0,0xAA));
        __m128 dst=_mm_load_ps(x);
        dst=_mm_sub_ps(dst,_mm512_castps512_ps128(tmp1));
        _mm_store_ps(x,dst);
#else
        alignas(64) float data[16];
        _mm512_store_ps(&data[0],data0);
        for(int k=0;k<16;k++) x[k%4]-=data[k];
#endif
    }
};

template<> struct Matrix4<float,__m512>
{
    typedef float T;
    typedef __m512 Tw;

    __m512 data0;

    void Load(const float (&A)[16])
    {
        data0=_mm512_load_ps(&A[0]);
    }

    template<int offset>
    void Load_Prefetch(const float (&A)[16])
    {
        data0=_mm512_load_ps(&A[0]);
    }

    template<bool broadcast_in_lane>
    Vector4<T,Tw,!broadcast_in_lane> operator*(const Vector4<T,Tw,broadcast_in_lane>& x) const
    {
        Vector4<T,Tw,!broadcast_in_lane> y;  
        y.data0=_mm512_mul_ps(data0,x.data0);
        return y;
    }

};

#endif
