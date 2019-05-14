//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_RSD_Positive_Definite_Part
#include "KernelCommon.h"
#else
namespace{
#endif

#ifdef FORCE_INLINE
#include "../Symmetric_Definite_Projection/Symmetric_Definite_Projection.h"
#else
#define SUBROUTINE_Symmetric_Definite_Projection
#include "../Symmetric_Definite_Projection/Symmetric_Definite_Projection.cpp"
#undef SUBROUTINE_Symmetric_Definite_Projection
#endif

namespace nm_RSD_Positive_Definite_Part {
    BUILD_CONSTANT(one_half,0.5f);
}

 using namespace nm_RSD_Positive_Definite_Part;

template<class Tw,class T_DATA=void,class I_DATA=void> inline void
RSD_Positive_Definite_Part(const T_DATA (&RSD)[12], T_DATA (&RSDpd)[12])
{
    typedef enum { x1111=0,x1122,x1133,x2222,x2233,x3333,x1212,x1221,x1313,x1331,x2323,x2332 } RSD_Entry;
    typedef Number<Tw> Tn;
    typedef const T_DATA (&refConstSymmMatrix)[6];
    typedef T_DATA (&refSymmMatrix)[6];

    Symmetric_Definite_Projection<Tw,T_DATA,I_DATA>(reinterpret_cast<refConstSymmMatrix>(RSD), reinterpret_cast<refSymmMatrix>(RSDpd));

    Tn val1,val2,val3,val_one_half;
    val_one_half.Load_Aligned(one_half);

    for(int i=0;i<5;i+=2){
        val1.Load(RSD[x1212+i]);
        val3.Load(RSD[x1221+i]);
        val2=val1-val3;
        val1=val1+val3;
        val3=Tn();
        val1=max(val1,val3);
        val3=max(val2,val3);
        val2=val1-val3;
        val1=val1+val3;
        val1=val1*val_one_half;
        val2=val2*val_one_half;
        Store(RSDpd[x1212+i],val1);
        Store(RSDpd[x1221+i],val2);
    }
}

#ifndef SUBROUTINE_RSD_Positive_Definite_Part
#define INSTANCE_KERNEL_RSD_Positive_Definite_Part(WIDTH) const WIDETYPE(float,WIDTH) (&RSD)[12], WIDETYPE(float,WIDTH) (&RSDpd)[12]
INSTANCE_KERNEL(RSD_Positive_Definite_Part);
#undef INSTANCE_KERNEL_RSD_Positive_Definite_Part
#else
}
#endif

