//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Deviatoric_Piola_Kirchhoff_Stress_Tensor
#include <assert.h>
#include <Common/KernelCommon.h>
#else
namespace{
#endif

#ifndef __Deviatoric_Piola_Kirchhoff_Stress_Tensor_defined__
#define __Deviatoric_Piola_Kirchhoff_Stress_Tensor_defined__
template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Deviatoric_Piola_Kirchhoff_Stress_Tensor
{
    static inline void Run(T_DATA (&DevP_Hat)[3], const T_DATA (&Sigma)[3], const T_DATA (&mu));
};
#endif

namespace {
    BUILD_CONSTANT(one,1.f);
    BUILD_CONSTANT(point_five,0.5f);
 }

template<class Tw,class T_DATA,class I_DATA>
struct Deviatoric_Piola_Kirchhoff_Stress_Tensor<COROTATED_TAG,Tw,T_DATA,I_DATA>
{
    static inline void Run(T_DATA (&DevP_Hat)[3], const T_DATA (&Sigma)[3], const T_DATA (&mu))
{
    typedef Number<Tw> Tn;

    Tn rMu;
    Tn rOne;
    Tn rSigma1;
    Tn rSigma2;
    Tn rSigma3;

    rMu.Load(mu);
    rOne.Load_Aligned(one);
    rSigma1.Load(Sigma[0]);
    rSigma2.Load(Sigma[1]);
    rSigma3.Load(Sigma[2]);  

    Tn rDevPHat1;
    Tn rDevPHat2;
    Tn rDevPHat3;
    
    rDevPHat1 = ( rSigma1 - rOne ) * ( rMu + rMu );   // (Sigma-1.)*(2.*mu)
    rDevPHat2 = ( rSigma2 - rOne ) * ( rMu + rMu );
    rDevPHat3 = ( rSigma3 - rOne ) * ( rMu + rMu );
    
    Store(DevP_Hat[0], rDevPHat1);
    Store(DevP_Hat[1], rDevPHat2);
    Store(DevP_Hat[2], rDevPHat3);
}
};

template<class Tw,class T_DATA,class I_DATA>
struct Deviatoric_Piola_Kirchhoff_Stress_Tensor<NEOHOOKEAN_TAG,Tw,T_DATA,I_DATA>
{
    static inline void Run(T_DATA (&DevP_Hat)[3], const T_DATA (&Sigma)[3], const T_DATA (&mu))
{
    typedef Number<Tw> Tn;

    Tn rMu;
    Tn rSigma1;
    Tn rSigma2;
    Tn rSigma3;

    rMu.Load(mu);
    rSigma1.Load(Sigma[0]);
    rSigma2.Load(Sigma[1]);
    rSigma3.Load(Sigma[2]);  

    Tn rDevPHat1;
    Tn rDevPHat2;
    Tn rDevPHat3;
    
    rDevPHat1 = rMu * ( rSigma1 - rSigma1.inverse() ); // mu*(Sigma-Sigma.Inverse());
    rDevPHat2 = rMu * ( rSigma2 - rSigma2.inverse() );
    rDevPHat3 = rMu * ( rSigma3 - rSigma3.inverse() );
    
    Store(DevP_Hat[0], rDevPHat1);
    Store(DevP_Hat[1], rDevPHat2);
    Store(DevP_Hat[2], rDevPHat3);
}
};

template<class Tw,class T_DATA,class I_DATA>
struct Deviatoric_Piola_Kirchhoff_Stress_Tensor<BIPHASIC_TAG,Tw,T_DATA,I_DATA>
{
    static inline void Run(T_DATA (&DevP_Hat)[3], const T_DATA (&Sigma)[3], const T_DATA (&mu))
{
    typedef Number<Tw> Tn;

    Tn rMu;
    Tn rZero;
    Tn rOne;
    Tn rBiphasicThreshold;
    Tn rSigma1;
    Tn rSigma2;
    Tn rSigma3;
    Tn rKbp;
    Tn rPointFive;

    rPointFive.Load(point_five);
    rKbp.Load(nm_biphasic::biphasic_factor);
    rMu.Load(mu);
    rOne.Load_Aligned(one);
    rBiphasicThreshold.Load(nm_biphasic::biphasic_threshold);
    rSigma1.Load(Sigma[0]);
    rSigma2.Load(Sigma[1]);
    rSigma3.Load(Sigma[2]);  

    Tn rDevPHat1;
    Tn rDevPHat2;
    Tn rDevPHat3;

    //rDevPHat1 = ( rSigma1 - rOne ) * ( rMu + rMu );   // (Sigma-1.)*(2.*mu)
    //rDevPHat2 = ( rSigma2 - rOne ) * ( rMu + rMu );
    //rDevPHat3 = ( rSigma3 - rOne ) * ( rMu + rMu );

    //rDevPHat1 = (max( ( rSigma1 - (rBiphasicThreshold+rOne) ), rZero ) * ( rKbp ));   // max((Sigma-gamma.),0)*(2.*mu)
    //rDevPHat2 = (max( ( rSigma2 - (rBiphasicThreshold+rOne) ), rZero ) * ( rKbp ));
    //rDevPHat3 = (max( ( rSigma3 - (rBiphasicThreshold+rOne) ), rZero ) * ( rKbp ));

    rDevPHat1 = (max( ( rSigma1 - (rBiphasicThreshold) ), rZero ));   // max((Sigma-gamma.),0)*(2.*mu)
    rDevPHat2 = (max( ( rSigma2 - (rBiphasicThreshold) ), rZero ));
    rDevPHat3 = (max( ( rSigma3 - (rBiphasicThreshold) ), rZero ));

    // Square the previous results, multiply by 1/2 Kbp, and add to standard Corotated energy term
    rDevPHat1 = (( rSigma1 - rOne ) * ( rMu + rMu )) + ( rDevPHat1 * rKbp );
    rDevPHat2 = (( rSigma2 - rOne ) * ( rMu + rMu )) + ( rDevPHat2 * rKbp );
    rDevPHat3 = (( rSigma3 - rOne ) * ( rMu + rMu )) + ( rDevPHat3 * rKbp );



    Store(DevP_Hat[0], rDevPHat1);
    Store(DevP_Hat[1], rDevPHat2);
    Store(DevP_Hat[2], rDevPHat3);
}
};

#ifndef SUBROUTINE_Deviatoric_Piola_Kirchhoff_Stress_Tensor
#define INSTANCE_KERNEL_Deviatoric_Piola_Kirchhoff_Stress_Tensor(WIDTH) WIDETYPE(float,WIDTH) (&DevP_Hat)[3],  const WIDETYPE(float,WIDTH) (&Sigma)[3], const WIDETYPE(float,WIDTH) (&mu) 
INSTANCE_KERNEL_MATERIAL(Deviatoric_Piola_Kirchhoff_Stress_Tensor,COROTATED_TAG);
INSTANCE_KERNEL_MATERIAL(Deviatoric_Piola_Kirchhoff_Stress_Tensor,NEOHOOKEAN_TAG);
#undef INSTANCE_KERNEL_Deviatoric_Piola_Kirchhoff_Stress_Tensor
#else
}
#endif
