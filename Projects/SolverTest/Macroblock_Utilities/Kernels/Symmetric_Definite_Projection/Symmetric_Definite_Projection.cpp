//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Symmetric_Definite_Projection
//#include <algorithm>
//#include <cmath>
#include "KernelCommon.h"
#else
namespace{
#endif

#include "Symmetric_Definite_Projection.h"

#define COMPUTE_V_AS_MATRIX

namespace{
    //float rsqrt(const float f)
    //{__m128 val=_mm_load_ss(&f);val=_mm_rsqrt_ss(val);float result;_mm_store_ss(&result,val);return result;}
}

 template<class Tw, class T_DATA, class I_DATA>
#ifdef SUBROUTINE_Symmetric_Definite_Projection
inline
#endif
void Symmetric_Definite_Projection(const T_DATA (&A)[6], T_DATA (&Apd)[6])
{

typedef enum { x11=0, x21, x31, x22, x32, x33 } Entry;

#include "Symmetric_Definite_Projection_Kernel_Declarations.hpp"

    Vs11.Load(A[x11]);
    Vs21.Load(A[x21]);
    Vs31.Load(A[x31]);
    Vs22.Load(A[x22]);
    Vs32.Load(A[x32]);
    Vs33.Load(A[x33]);

#include "Symmetric_Definite_Projection_Main_Kernel_Body.hpp"    

    Store(Apd[x11],Vs11);
    Store(Apd[x21],Vs21);
    Store(Apd[x31],Vs31);
    Store(Apd[x22],Vs22);    
    Store(Apd[x32],Vs32);
    Store(Apd[x33],Vs33);
}

#undef COMPUTE_V_AS_MATRIX

#ifndef SUBROUTINE_Symmetric_Definite_Projection
#define INSTANCE_KERNEL_Symmetric_Definite_Projection(WIDTH) const WIDETYPE(float,WIDTH) (&A)[6], WIDETYPE(float,WIDTH) (&Apd)[6]
INSTANCE_KERNEL(Symmetric_Definite_Projection);
#undef INSTANCE_KERNEL_Symmetric_Definite_Projection
#else
}
#endif
