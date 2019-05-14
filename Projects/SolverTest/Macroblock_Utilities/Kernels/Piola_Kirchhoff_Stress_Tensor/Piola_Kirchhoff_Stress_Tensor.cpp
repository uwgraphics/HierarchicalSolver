//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Piola_Kirchhoff_Stress_Tensor
#include <assert.h>
#include <Common/KernelCommon.h>
#else
namespace {
#endif

#ifdef FORCE_INLINE
#include "../Deviatoric_Piola_Kirchhoff_Stress_Tensor/Deviatoric_Piola_Kirchhoff_Stress_Tensor.h"
#include "../Penalty_Measure_Gradient/Penalty_Measure_Gradient.h"
#ifdef USE_NONMIXED_FORMULAS
#include "../Volume_Preservation_Deviation/Volume_Preservation_Deviation.h"
#endif
#else
#define SUBROUTINE_Deviatoric_Piola_Kirchhoff_Stress_Tensor
#include "../Deviatoric_Piola_Kirchhoff_Stress_Tensor/Deviatoric_Piola_Kirchhoff_Stress_Tensor.cpp"
#undef SUBROUTINE_Deviatoric_Piola_Kirchhoff_Stress_Tensor

#define SUBROUTINE_Penalty_Measure_Gradient
#include "../Penalty_Measure_Gradient/Penalty_Measure_Gradient.cpp"
#undef SUBROUTINE_Penalty_Measure_Gradient

#ifdef USE_NONMIXED_FORMULAS
#define SUBROUTINE_Volume_Preservation_Deviation
#include "../Volume_Preservation_Deviation/Volume_Preservation_Deviation.cpp"
#undef SUBROUTINE_Volume_Preservation_Deviation
#endif

#endif

template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Piola_Kirchhoff_Stress_Tensor
{
static 
#ifdef SUBROUTINE_Piola_Kirchhoff_Stress_Tensor
inline
#endif
void Run(T_DATA (&P_Hat)[3], const T_DATA (&Sigma)[3], const T_DATA (&Q_Hat)[3], const T_DATA (&p), const T_DATA (&mu), const T_DATA (&alpha), const T_DATA (&kappa))
{
    typedef Number<Tw> Tn;

    T_DATA DevP_Hat[3];
    Deviatoric_Piola_Kirchhoff_Stress_Tensor<T_MATERIAL,Tw,T_DATA,I_DATA>::Run(DevP_Hat,Sigma,mu);
    
    #ifdef USE_NONMIXED_FORMULAS
    T_DATA M;
    Volume_Preservation_Deviation<T_MATERIAL,Tw,T_DATA,I_DATA>::Run(M, Sigma);

    Tn rDevP_Hat;
    Tn rM;
    Tn rKappa;
    //Tn rQ_Hat;
    Tn rP_Hat;
    rKappa.Load( kappa );
    rM.Load( M ); 

    for(int i=0; i<3; i++){
        rDevP_Hat.Load(DevP_Hat[i]);
        //rQ_Hat.Load(Q_Hat[i]);
        rP_Hat = rDevP_Hat + (rKappa * rM ); // TODO: Use a real QHat later for other non-corotated materials * rQ_Hat);
        Store(P_Hat[i], rP_Hat);
    }
    #else
    Tn rp;
    Tn rAlpha;
    Tn rDevP_Hat;
    Tn rQ_Hat;
    Tn rP_Hat;
    rp.Load( p );
    rAlpha.Load( alpha );

    for(int i=0; i<3; i++){
        rDevP_Hat.Load(DevP_Hat[i]);
        rQ_Hat.Load(Q_Hat[i]);
        rP_Hat = rDevP_Hat + rAlpha * rp * rQ_Hat;
        Store(P_Hat[i], rP_Hat);
    }
    #endif
}
};

#ifndef SUBROUTINE_Piola_Kirchhoff_Stress_Tensor
#define INSTANCE_KERNEL_Piola_Kirchhoff_Stress_Tensor(WIDTH) WIDETYPE(float,WIDTH) (&P_Hat)[3],    const WIDETYPE(float,WIDTH) (&Sigma)[3],    const WIDETYPE(float,WIDTH) (&Q_Hat)[3],    const WIDETYPE(float,WIDTH) (&p),    const WIDETYPE(float,WIDTH) (&mu),    const WIDETYPE(float,WIDTH) (&alpha),    const WIDETYPE(float,WIDTH) (&kappa)
INSTANCE_KERNEL_MATERIAL(Piola_Kirchhoff_Stress_Tensor, COROTATED_TAG);
INSTANCE_KERNEL_MATERIAL(Piola_Kirchhoff_Stress_Tensor, NEOHOOKEAN_TAG);
#undef INSTANCE_KERNEL_Piola_Kirchhoff_Stress_Tensor
#else
}
#endif
