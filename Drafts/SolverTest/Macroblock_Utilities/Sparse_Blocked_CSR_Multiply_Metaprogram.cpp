//#####################################################################
// Copyright 2012-2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef SUBROUTINE_Sparse_Blocked_CSR_Multiply_Metaprogram
#include "Common/KernelCommon.h"
#include "BlockedCSRMatrix3.h"
#include "DiagonalMatrix3.h"
#include "SymmetricMatrix3.h"
#include "Matrix3.h"
#include "UnitriangularMatrix3.h"
#include <cassert>
#include <fstream>
#include <iostream>
#else
namespace {
#endif

template<class Tw,class T_DATA,class I_DATA>
#ifdef SUBROUTINE_Sparse_Blocked_CSR_Multiply_Metaprogram
inline
#endif
void Sparse_Blocked_CSR_Multiply_Metaprogram(const BlockedCSRMatrix3<T_DATA> &U,
                                 const typename std::remove_extent<T_DATA>::type *b,
                                 typename std::remove_extent<T_DATA>::type *x,
                                 bool transpose)
{
    abort();
#if 0
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
#endif
}



template<class Tw,class T_DATA,class I_DATA>
#ifdef SUBROUTINE_Sparse_Blocked_CSR_Multiply_Metaprogram
inline
#endif
void Sparse_Blocked_CSR_Masked_Multiply_Metaprogram(const BlockedCSRMatrix3<T_DATA> &U,
                                        const typename std::remove_extent<T_DATA>::type *b,
                                        typename std::remove_extent<T_DATA>::type *x,
                                        const T_DATA (&mask),
                                        bool transpose)
{
    // Print header

    const int T_SIZE = std::rank<T_DATA>::value == 0 ? 1 : std::extent<T_DATA>::value;
    std::cout<<"T_SIZE: "<<T_SIZE<<std::endl; 
    std::ofstream myfile;
    myfile.open("test.cpp");
    
    
    myfile <<
        "#include \"Macroblock_Utilities/Common/KernelCommon.h\""                                               << std::endl <<
        "#include \"Macroblock_Utilities/BlockedCSRMatrix3.h\""                                                 << std::endl <<
        "#include \"Macroblock_Utilities/Matrix3.h\""                                                           << std::endl <<
        "void Sparse_Blocked_CSR_Masked_Multiply_Flat(const BlockedCSRMatrix3<float[16]> &U,"                   << std::endl <<
        "                                             const typename std::remove_extent<float[16]>::type *b,"   << std::endl <<
        "                                             typename std::remove_extent<float[16]>::type *x,"         << std::endl <<
        "                                             const float (&mask)[16],"                                 << std::endl <<
        "                                             bool transpose)"                                          << std::endl <<
        "{"                                                                                                     << std::endl <<
        "    typedef Number<__m512> Tn;"                                                                      << std::endl <<
        "    typedef float (&VEC3_TYPE)[3][16];"                                                              << std::endl <<
        "    typedef const float (&CVEC3_TYPE)[3][16]; "                                                      << std::endl <<
        "    typedef const float (*v3c_ptr)[3][16];"                                                            << std::endl <<
        "    typedef float (*v3_ptr)[3][16];"                                                                   << std::endl;

#if 1
    myfile <<
        "    __m512 vMask = _mm512_load_ps(mask);"                                                              << std::endl;
    for(int i=0;i<U.rows;i++){
        myfile <<
            "    "     << "{"                                                                                     << std::endl <<
            "        " << "__m512 vM11_e = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM21_e = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM31_e = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM12_e = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM22_e = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM32_e = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM13_e = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM23_e = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM33_e = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM11_o = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM21_o = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM31_o = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM12_o = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM22_o = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM32_o = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM13_o = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM23_o = _mm512_setzero_ps();"                                                  << std::endl <<
            "        " << "__m512 vM33_o = _mm512_setzero_ps();"                                                  << std::endl;
        for(int jj=U.offsets[i];jj<U.offsets[i+1];jj++){
            std::string parity = (jj%2) ? "o" : "e";
            int j=U.column[jj];
            myfile << "        " << "{"                                                                                     << std::endl;
            for(int w=0;w<3;w++) myfile <<
                "            " << "__m512 vX_"<<j<<"_"<<w+1<<" = _mm512_load_ps(reinterpret_cast<v3c_ptr>(b)["<<j<<"]["<<w<<"]);"       << std::endl;
            for(int w=0;w<3;w++) for(int v=0;v<3;v++) myfile <<                  
                "            " << "vM" << v+1 << w+1 << "_" << parity << " = _mm512_fmadd_ps(" <<
                                      " _mm512_load_ps(U.E("<<jj<<")["<<v+3*w<<"]) ," <<
                                      " vX_"<<j<<"_"<<w+1<<" ," <<
                                                          " vM" << v+1 << w+1 << "_" << parity << " );" << std::endl;
            myfile << "        " << "}"                                                                                     << std::endl;
            std::cout << "(i,j)=("<<i<<","<<U.column[jj]<<")"<<std::endl;
        }

        myfile << "        __m512 vY_1 = vM11_e + vM12_e + vM13_e + vM11_o + vM12_o + vM13_o;" << std::endl;
        myfile << "        __m512 vY_2 = vM21_e + vM22_e + vM23_e + vM21_o + vM22_o + vM23_o;" << std::endl;
        myfile << "        __m512 vY_3 = vM31_e + vM32_e + vM33_e + vM31_o + vM32_o + vM33_o;" << std::endl;

        myfile <<
            "        " << "vY_1 = _mm512_mul_ps(vMask,vY_1);"                                             << std::endl <<
            "        " << "vY_2 = _mm512_mul_ps(vMask,vY_2);"                                             << std::endl <<
            "        " << "vY_3 = _mm512_mul_ps(vMask,vY_3);"                                             << std::endl <<
            "        " << "_mm512_store_ps(reinterpret_cast<v3_ptr>(x)["<<i<<"][0],vY_1);"                        << std::endl <<
            "        " << "_mm512_store_ps(reinterpret_cast<v3_ptr>(x)["<<i<<"][1],vY_2);"                        << std::endl <<
            "        " << "_mm512_store_ps(reinterpret_cast<v3_ptr>(x)["<<i<<"][2],vY_3);"                        << std::endl <<
            "    " << "}"                                                                                             << std::endl;
    }
#else
    myfile <<
        "    typedef Number<__m512> Tn;"                                                                      << std::endl <<
        "    typedef float (&VEC3_TYPE)[3][16];"                                                              << std::endl <<
        "    typedef const float (&CVEC3_TYPE)[3][16]; "                                                      << std::endl <<
        "    Tn vMask;"                                                                                       << std::endl <<  
        "    vMask.Load_Aligned( mask );"                                                                     << std::endl;
    for(int i=0;i<U.rows;i++){
        myfile <<
            "    {"                                                                                           << std::endl <<
            "    " << "Vector3<Tn> x_"<<i<<";"                                                                << std::endl;
        for(int jj=U.offsets[i];jj<U.offsets[i+1];jj++){
            int j=U.column[jj];
            myfile <<
                "    "     << "{"                                                                                 << std::endl <<
                "        " << "Matrix3<Tn> U_"<<i<<j<<";"                                                         << std::endl <<
                "        " << "U_"<<i<<""<<j<<".Load_Aligned(U.E("<<jj<<"));"                                     << std::endl <<
                "        " << "CVEC3_TYPE b_"<<j<<"_raw = reinterpret_cast<CVEC3_TYPE>(*(b+("<<j<<"*3+0)*16));"   << std::endl <<
                "        " << "Vector3<Tn> b_"<<j<<";"                                                            << std::endl <<
                "        " << "b_"<<j<<".Load_Aligned( b_"<<j<<"_raw );"                                          << std::endl <<
                "        " << "x_"<<i<<" = x_"<<i<<" + ( U_"<<i<<j<<" * b_"<<j<<" );"                             << std::endl <<
                "    "     << "}"                                                                                 << std::endl;           
        }
        myfile <<
            "    " << "VEC3_TYPE x_"<<i<<"_raw = reinterpret_cast<VEC3_TYPE>(*(x+("<<i<<"*3+0)*16));"             << std::endl <<
            "    " << "x_"<<i<<" *= vMask;"                                                                       << std::endl <<
            "    " << "x_"<<i<<".Store( x_"<<i<<"_raw );"                                                         << std::endl <<
            "    }"                                                                                               << std::endl;
    }
#endif
    // Print footer
    myfile <<
        "}"                                                                                                  << std::endl;

    myfile.close();

    exit(0);

#if 0
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
#endif
}

#ifndef SUBROUTINE_Sparse_Blocked_CSR_Multiply_Metaprogram

#define INSTANCE_KERNEL_Sparse_Blocked_CSR_Multiply_Metaprogram(WIDTH) const BlockedCSRMatrix3< WIDETYPE(float,WIDTH) > &U, const float *b, float *x, bool transpose
INSTANCE_KERNEL(Sparse_Blocked_CSR_Multiply_Metaprogram);
#undef INSTANCE_KERNEL_Sparse_Blocked_CSR_Multiply_Metaprogram
    
#define INSTANCE_KERNEL_Sparse_Blocked_CSR_Multiply_Metaprogram(WIDTH) const BlockedCSRMatrix3< WIDETYPE(double,WIDTH) > &U, const double *b, double *x, bool transpose
INSTANCE_KERNEL_DOUBLE(Sparse_Blocked_CSR_Multiply_Metaprogram);
#undef INSTANCE_KERNEL_Sparse_Blocked_CSR_Multiply_Metaprogram
    
#define INSTANCE_KERNEL_Sparse_Blocked_CSR_Masked_Multiply_Metaprogram(WIDTH) const BlockedCSRMatrix3< WIDETYPE(float,WIDTH) > &U, const float *b, float *x, const WIDETYPE(float,WIDTH) (&mask), bool transpose
INSTANCE_KERNEL(Sparse_Blocked_CSR_Masked_Multiply_Metaprogram);
#undef INSTANCE_KERNEL_Sparse_Blocked_CSR_Masked_Multiply_Metaprogram

#define INSTANCE_KERNEL_Sparse_Blocked_CSR_Masked_Multiply_Metaprogram(WIDTH) const BlockedCSRMatrix3< WIDETYPE(double,WIDTH) > &U, const double *b, double *x, const WIDETYPE(double,WIDTH) (&mask), bool transpose
INSTANCE_KERNEL_DOUBLE(Sparse_Blocked_CSR_Masked_Multiply_Metaprogram);
#undef INSTANCE_KERNEL_Sparse_Blocked_CSR_Masked_Multiply_Metaprogram

#else
}
#endif

