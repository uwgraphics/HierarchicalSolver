//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Penalty_Measure_Gradient
#include <assert.h>
#include <Common/KernelCommon.h>
#else
namespace{
#endif

    namespace nm_Penalty_Measure_Gradient{
        BUILD_CONSTANT(one,1.f);
    }

    using namespace nm_Penalty_Measure_Gradient;

template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Penalty_Measure_Gradient
{
static inline void Run(const T_DATA (&Sigma)[3], T_DATA (&Q_hat)[3]);
};

template<class Tw,class T_DATA,class I_DATA>
struct Penalty_Measure_Gradient<COROTATED_TAG,Tw,T_DATA,I_DATA>
{
static inline void Run(const T_DATA (&Sigma)[3], T_DATA (&Q_hat)[3])
{
    typedef Number<Tw> Tn;

    Tn rone;
    rone.Load_Aligned(nm_Penalty_Measure_Gradient::one);
    Store(Q_hat[0],rone);
    Store(Q_hat[1],rone);
    Store(Q_hat[2],rone);
}
};

template<class Tw,class T_DATA,class I_DATA>
struct Penalty_Measure_Gradient<NEOHOOKEAN_TAG,Tw,T_DATA,I_DATA>
{
static inline void Run(const T_DATA (&Sigma)[3], T_DATA (&Q_hat)[3])
{
    typedef Number<Tw> Tn;

    Tn rSigma1_inverse;
    Tn rSigma2_inverse;
    Tn rSigma3_inverse;

    rSigma1_inverse.Load(Sigma[0]);
    rSigma2_inverse.Load(Sigma[1]);
    rSigma3_inverse.Load(Sigma[2]);

    rSigma1_inverse=rSigma1_inverse.inverse();
    rSigma2_inverse=rSigma2_inverse.inverse();
    rSigma3_inverse=rSigma3_inverse.inverse();

    Store(Q_hat[0],rSigma1_inverse);
    Store(Q_hat[1],rSigma2_inverse);
    Store(Q_hat[2],rSigma3_inverse);
}
};

template<class Tw,class T_DATA,class I_DATA>
struct Penalty_Measure_Gradient<BIPHASIC_TAG,Tw,T_DATA,I_DATA>
{
static inline void Run(const T_DATA (&Sigma)[3], T_DATA (&Q_hat)[3])
{
    typedef Number<Tw> Tn;

    Tn rone;
    rone.Load_Aligned(nm_Penalty_Measure_Gradient::one);
    Store(Q_hat[0],rone);
    Store(Q_hat[1],rone);
    Store(Q_hat[2],rone);
}
};

#ifndef SUBROUTINE_Penalty_Measure_Gradient
#define INSTANCE_KERNEL_Penalty_Measure_Gradient(WIDTH)  const WIDETYPE(float,WIDTH) (&Sigma)[3],WIDETYPE(float,WIDTH) (&Q_hat)[3]
INSTANCE_KERNEL_MATERIAL(Penalty_Measure_Gradient,COROTATED_TAG);
INSTANCE_KERNEL_MATERIAL(Penalty_Measure_Gradient,NEOHOOKEAN_TAG);
#undef INSTANCE_KERNEL_Penalty_Measure_Gradient
#else
}
#endif
