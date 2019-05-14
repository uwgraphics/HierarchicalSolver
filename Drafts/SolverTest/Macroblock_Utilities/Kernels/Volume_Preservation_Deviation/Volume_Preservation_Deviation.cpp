//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Volume_Preservation_Deviation
#include <assert.h>
#include <Common/KernelCommon.h>
#else
namespace{
#endif

#ifndef __Volume_Preservation_Deviation_defined__
#define __Volume_Preservation_Deviation_defined__
template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Volume_Preservation_Deviation
{
static inline void Run(T_DATA (&M), const T_DATA (&Sigma)[3]);
};


 namespace {
     BUILD_CONSTANT(three,3.f);
 }

template<class Tw,class T_DATA,class I_DATA>
struct Volume_Preservation_Deviation<COROTATED_TAG,Tw,T_DATA,I_DATA>
{
static inline void Run(T_DATA (&M), const T_DATA (&Sigma)[3])
{
    typedef Number<Tw> Tn;

    Tn rSigma1;
    Tn rSigma2;
    Tn rSigma3;
    Tn rThree;

    rSigma1.Load(Sigma[0]);
    rSigma2.Load(Sigma[1]);
    rSigma3.Load(Sigma[2]);
    rThree.Load_Aligned(three);

    Tn rM;
    
    rM = rSigma1 + rSigma2 + rSigma3 - rThree;  // (Sigma-1).Trace();
    
    Store(M, rM);
}
};

template<class Tw,class T_DATA,class I_DATA>
struct Volume_Preservation_Deviation<NEOHOOKEAN_TAG,Tw,T_DATA,I_DATA>
{
static inline void Run(T_DATA (&M), const T_DATA (&Sigma)[3])
{
    typedef Number<Tw> Tn;

    Tn rSigma1;
    Tn rSigma2;
    Tn rSigma3;

    rSigma1.Load(Sigma[0]);
    rSigma2.Load(Sigma[1]);
    rSigma3.Load(Sigma[2]);

    Tn rM;
    
    rM = (rSigma1 * rSigma2 * rSigma3).log(); // log( Sigma.Determinate() )
    
    Store(M, rM);
}
};

template<class Tw,class T_DATA,class I_DATA>
struct Volume_Preservation_Deviation<BIPHASIC_TAG,Tw,T_DATA,I_DATA>
{
static inline void Run(T_DATA (&M), const T_DATA (&Sigma)[3])
{
    typedef Number<Tw> Tn;

    Tn rZero;
    Tn rSigma1;
    Tn rSigma2;
    Tn rSigma3;
    Tn rThree;
    Tn rBiphasicThreshold;

    rBiphasicThreshold.Load(nm_biphasic::biphasic_threshold);
    rSigma1.Load(Sigma[0]);
    rSigma2.Load(Sigma[1]);
    rSigma3.Load(Sigma[2]);
    rThree.Load_Aligned(three);

    Tn rM;
    
	rM = rSigma1 + rSigma2 + rSigma3 - rThree;  // (Sigma-1).Trace();
    //rM = max( rSigma1-rBiphasicThreshold, rZero ) +
    //     max( rSigma2-rBiphasicThreshold, rZero ) +
    //     max( rSigma3-rBiphasicThreshold, rZero );
    
    Store(M, rM);
}
};
#endif

#ifndef SUBROUTINE_Volume_Preservation_Deviation
#define INSTANCE_KERNEL_Volume_Preservation_Deviation(WIDTH) WIDETYPE(float,WIDTH) (&M), const WIDETYPE(float,WIDTH) (&Sigma)[3]
INSTANCE_KERNEL_MATERIAL(Volume_Preservation_Deviation,COROTATED_TAG);
INSTANCE_KERNEL_MATERIAL(Volume_Preservation_Deviation,NEOHOOKEAN_TAG);
#undef INSTANCE_KERNEL_Volume_Preservation_Deviation
#else
}
#endif
