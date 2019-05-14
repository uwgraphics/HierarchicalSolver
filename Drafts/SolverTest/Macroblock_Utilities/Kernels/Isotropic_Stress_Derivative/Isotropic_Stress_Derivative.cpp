//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Isotropic_Stress_Derivative
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif

#ifdef FORCE_INLINE
#include "../Augmented_Rotated_Stress_Derivative/Augmented_Rotated_Stress_Derivative.h"
#include "../Rotated_Stress_Derivative/Rotated_Stress_Derivative.h"
#include "../RSD_Positive_Definite_Part/RSD_Positive_Definite_Part.h"
#else
#define SUBROUTINE_Augmented_Rotated_Stress_Derivative
#include "../Augmented_Rotated_Stress_Derivative/Augmented_Rotated_Stress_Derivative.cpp"
#undef SUBROUTINE_Augmented_Rotated_Stress_Derivative

#define SUBROUTINE_Rotated_Stress_Derivative
#include "../Rotated_Stress_Derivative/Rotated_Stress_Derivative.cpp"
#undef SUBROUTINE_Rotated_Stress_Derivative

#define SUBROUTINE_RSD_Positive_Definite_Part
#include "../RSD_Positive_Definite_Part/RSD_Positive_Definite_Part.cpp"
#undef SUBROUTINE_RSD_Positive_Definite_Part
#endif


template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Isotropic_Stress_Derivative
{
static 
#ifdef SUBROUTINE_Isotropic_Stress_Derivative
inline
#endif
void Run(T_DATA (&dPdF)[12], const T_DATA (&Sigma)[3], const T_DATA (&p), const T_DATA (&mu), const T_DATA (&kappa), const T_DATA (&alpha), const bool apply_definiteness_fix)
{
    typedef Number<Tw> Tn;

#ifdef USE_NONMIXED_FORMULAS
    Rotated_Stress_Derivative<T_MATERIAL,Tw,T_DATA,I_DATA>::Run(dPdF,Sigma,mu,kappa);
    RSD_Positive_Definite_Part<Tw,T_DATA,I_DATA>(dPdF,dPdF);
#else
    Augmented_Rotated_Stress_Derivative<T_MATERIAL,Tw,T_DATA,I_DATA>::Run(dPdF,Sigma,p,mu,alpha);

    if(apply_definiteness_fix){
        BUILD_TDATA(dPdF_unaugmented,[12]);
        Rotated_Stress_Derivative<T_MATERIAL,Tw,T_DATA,I_DATA>::Run(dPdF_unaugmented,Sigma,mu,kappa);
        BUILD_TDATA(dPdF_spd_unaugmented,[12]);
        RSD_Positive_Definite_Part<Tw,T_DATA,I_DATA>(dPdF_unaugmented,dPdF_spd_unaugmented);

        Tn val1,val2,val3;
        for(int i=0;i<12;i++){
            val1.Load(dPdF[i]);
            val2.Load(dPdF_unaugmented[i]);
            val3.Load(dPdF_spd_unaugmented[i]);
            val1=val1-val2+val3;
            Store(dPdF[i],val1);
        }
    }
#endif

}
};


#ifndef SUBROUTINE_Isotropic_Stress_Derivative
#define INSTANCE_KERNEL_Isotropic_Stress_Derivative(WIDTH) WIDETYPE(float,WIDTH) (&dPdF)[12], const WIDETYPE(float,WIDTH) (&Sigma)[3],    const WIDETYPE(float,WIDTH) (&p),    const WIDETYPE(float,WIDTH) (&mu),    const WIDETYPE(float,WIDTH) (&kappa),    const WIDETYPE(float,WIDTH) (&alpha), const bool apply_definiteness_fix
INSTANCE_KERNEL_MATERIAL(Isotropic_Stress_Derivative,COROTATED_TAG);
INSTANCE_KERNEL_MATERIAL(Isotropic_Stress_Derivative,NEOHOOKEAN_TAG);
INSTANCE_KERNEL_MATERIAL(Isotropic_Stress_Derivative,BIPHASIC_TAG);
#undef INSTANCE_KERNEL_Isotropic_Stress_Derivative
#else
}
#endif
