//#####################################################################
// Copyright 2012-2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef SUBROUTINE_Sparse_Blocked_CSR_Multiply
#include "Common/KernelCommon.h"
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
#ifdef SUBROUTINE_Sparse_Blocked_CSR_Multiply
inline
#endif
void Sparse_Blocked_CSR_Multiply(const BlockedCSRMatrix3<T_DATA> &U,
                                 const typename std::remove_extent<T_DATA>::type *b,
                                 typename std::remove_extent<T_DATA>::type *x,
                                 bool transpose)
{
    typedef Number<Tw> Tn;
    typedef typename std::remove_extent<T_DATA>::type T_Base;
    typedef typename std::remove_extent<I_DATA>::type I_Base;
    static_assert( std::rank<T_DATA>::value == std::rank<I_DATA>::value, "Error: Base datatypes mismatch in width." );
    static_assert( std::rank<T_DATA>::value == 1 ||
                   std::rank<T_DATA>::value == 0 , "Error: Base datatypes must be a rank 1 array or a scalar." );

    const int T_SIZE = std::rank<T_DATA>::value == 0 ? 1 : std::extent<T_DATA>::value;
    typedef T_Base (&VEC3_TYPE)[3][T_SIZE];
    typedef const T_Base (&CVEC3_TYPE)[3][T_SIZE];

    if( !transpose ){
        for( int i = 0; i < U.rows; i++ ){
            Vector3<Tn> x_i;
            for( int jj = U.offsets[i]; jj < U.offsets[i+1]; jj++ ){
                int j = U.column[jj];
                Matrix3<Tn> U_ij;
                U_ij.Load_Aligned(U.E(jj));
                CVEC3_TYPE b_j_raw = reinterpret_cast<CVEC3_TYPE>(*(b+(j*3+0)*T_SIZE));
                Vector3<Tn> b_j;
                b_j.Load_Aligned( b_j_raw );
                x_i = x_i + ( U_ij * b_j );
            }
            VEC3_TYPE x_i_raw = reinterpret_cast<VEC3_TYPE>(*(x+(i*3+0)*T_SIZE));
            x_i.Store( x_i_raw );
        }
    }
    else{
        // Store zeroes first
        for( int j = 0; j < U.columns; j++){
            VEC3_TYPE x_j_raw = reinterpret_cast<VEC3_TYPE>(*(x+(j*3+0)*T_SIZE));
            Vector3<Tn> x_j;
            x_j.Store( x_j_raw );
        }
        
        for( int i = 0; i < U.rows; i++ ){
            for( int jj = U.offsets[i]; jj < U.offsets[i+1]; jj++ ){
                int j = U.column[jj];
                Matrix3<Tn> U_ij;
                U_ij.Load_Aligned_Prefetch(U.E(jj),0);

                CVEC3_TYPE b_i_raw = reinterpret_cast<CVEC3_TYPE>(*(b+(i*3+0)*T_SIZE));
                Vector3<Tn> b_i;
                b_i.Load_Aligned( b_i_raw );

                VEC3_TYPE x_j_raw = reinterpret_cast<VEC3_TYPE>(*(x+(j*3+0)*T_SIZE));
                Vector3<Tn> x_j;
                x_j.Load_Aligned( x_j_raw );

                x_j = x_j + ( U_ij.Times_Transpose(b_i) );

                x_j.Store( x_j_raw );
            }
        }
    }
}



template<class Tw,class T_DATA,class I_DATA>
#ifdef SUBROUTINE_Sparse_Blocked_CSR_Multiply
inline
#endif
void Sparse_Blocked_CSR_Masked_Multiply(const BlockedCSRMatrix3<T_DATA> &U,
                                        const typename std::remove_extent<T_DATA>::type *b,
                                        typename std::remove_extent<T_DATA>::type *x,
                                        const T_DATA (&mask),
                                        bool transpose)
{
    typedef Number<Tw> Tn;
    typedef typename std::remove_extent<T_DATA>::type T_Base;
    typedef typename std::remove_extent<I_DATA>::type I_Base;
    static_assert( std::rank<T_DATA>::value == std::rank<I_DATA>::value, "Error: Base datatypes mismatch in width." );
    static_assert( std::rank<T_DATA>::value == 1 ||
                   std::rank<T_DATA>::value == 0 , "Error: Base datatypes must be a rank 1 array or a scalar." );

    const int T_SIZE = std::rank<T_DATA>::value == 0 ? 1 : std::extent<T_DATA>::value;
    typedef T_Base (&VEC3_TYPE)[3][T_SIZE];
    typedef const T_Base (&CVEC3_TYPE)[3][T_SIZE];

    Tn vMask;
    vMask.Load_Aligned( mask );

    if( !transpose ){
        for( int i = 0; i < U.rows; i++ ){
            Vector3<Tn> x_i;
            for( int jj = U.offsets[i]; jj < U.offsets[i+1]; jj++ ){
                int j = U.column[jj];
                Matrix3<Tn> U_ij;
                U_ij.Load_Aligned(U.E(jj));
                CVEC3_TYPE b_j_raw = reinterpret_cast<CVEC3_TYPE>(*(b+(j*3+0)*T_SIZE));
                Vector3<Tn> b_j;
                b_j.Load_Aligned( b_j_raw );
                x_i = x_i + ( U_ij * b_j );
            }
            VEC3_TYPE x_i_raw = reinterpret_cast<VEC3_TYPE>(*(x+(i*3+0)*T_SIZE));
            x_i *= vMask;
            x_i.Store( x_i_raw );
        }
    }
    else{
        // Store zeroes first
        for( int j = 0; j < U.columns; j++){
            VEC3_TYPE x_j_raw = reinterpret_cast<VEC3_TYPE>(*(x+(j*3+0)*T_SIZE));
            Vector3<Tn> x_j;
            x_j.Store( x_j_raw );
        }
        
        for( int i = 0; i < U.rows; i++ ){
            for( int jj = U.offsets[i]; jj < U.offsets[i+1]; jj++ ){
                int j = U.column[jj];
                Matrix3<Tn> U_ij;
                U_ij.Load_Aligned_Prefetch(U.E(jj),0);

                CVEC3_TYPE b_i_raw = reinterpret_cast<CVEC3_TYPE>(*(b+(i*3+0)*T_SIZE));
                Vector3<Tn> b_i;
                b_i.Load_Aligned( b_i_raw );

                VEC3_TYPE x_j_raw = reinterpret_cast<VEC3_TYPE>(*(x+(j*3+0)*T_SIZE));
                Vector3<Tn> x_j;
                x_j.Load_Aligned( x_j_raw );

                x_j = x_j + ( U_ij.Times_Transpose(b_i) );
                x_j *= vMask;
                x_j.Store( x_j_raw );
            }
        }
    }
}

#ifndef SUBROUTINE_Sparse_Blocked_CSR_Multiply

#define INSTANCE_KERNEL_Sparse_Blocked_CSR_Multiply(WIDTH) const BlockedCSRMatrix3< WIDETYPE(float,WIDTH) > &U, const float *b, float *x, bool transpose
INSTANCE_KERNEL(Sparse_Blocked_CSR_Multiply);
#undef INSTANCE_KERNEL_Sparse_Blocked_CSR_Multiply
    
#define INSTANCE_KERNEL_Sparse_Blocked_CSR_Multiply(WIDTH) const BlockedCSRMatrix3< WIDETYPE(double,WIDTH) > &U, const double *b, double *x, bool transpose
INSTANCE_KERNEL_DOUBLE(Sparse_Blocked_CSR_Multiply);
#undef INSTANCE_KERNEL_Sparse_Blocked_CSR_Multiply
    
#define INSTANCE_KERNEL_Sparse_Blocked_CSR_Masked_Multiply(WIDTH) const BlockedCSRMatrix3< WIDETYPE(float,WIDTH) > &U, const float *b, float *x, const WIDETYPE(float,WIDTH) (&mask), bool transpose
INSTANCE_KERNEL(Sparse_Blocked_CSR_Masked_Multiply);
#undef INSTANCE_KERNEL_Sparse_Blocked_CSR_Masked_Multiply

#define INSTANCE_KERNEL_Sparse_Blocked_CSR_Masked_Multiply(WIDTH) const BlockedCSRMatrix3< WIDETYPE(double,WIDTH) > &U, const double *b, double *x, const WIDETYPE(double,WIDTH) (&mask), bool transpose
INSTANCE_KERNEL_DOUBLE(Sparse_Blocked_CSR_Masked_Multiply);
#undef INSTANCE_KERNEL_Sparse_Blocked_CSR_Masked_Multiply

#else
}
#endif

