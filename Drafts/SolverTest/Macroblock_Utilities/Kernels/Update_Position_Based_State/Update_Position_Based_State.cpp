//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <assert.h>

#include "../../Common/KernelCommon.h"

#ifdef FORCE_INLINE
#include "../Isotropic_Stress_Derivative/Isotropic_Stress_Derivative.h"
#include "../Weighted_Gradient/Weighted_Gradient.h"
#include "../Singular_Value_Decomposition/Singular_Value_Decomposition.h"
namespace {
    BUILD_CONSTANT(one,1.f);
}
#else

#define SUBROUTINE_Isotropic_Stress_Derivative
#include "../Isotropic_Stress_Derivative/Isotropic_Stress_Derivative.cpp"
#undef SUBROUTINE_Isotropic_Stress_Derivative

#define SUBROUTINE_Unweighted_Gradient
#include "../Weighted_Gradient/Weighted_Gradient.cpp"
#undef SUBROUTINE_Unweighted_Gradient

#define SUBROUTINE_Singular_Value_Decomposition
#include "../Singular_Value_Decomposition/Singular_Value_Decomposition.cpp"
#undef SUBROUTINE_Singular_Value_Decomposition

#endif

template<class T_MATERIAL,class Tw,class T_DATA=void, class I_DATA=void> 
struct Update_Position_Based_State
{
    static void Run(const T_DATA (&u)[3][8], 
                    const T_DATA (&mu),
                    const T_DATA (&lambda),
                    const T_DATA (&weights)[3],
                    const T_DATA (&cutoff),
                    const T_DATA (&one_over_h)[3],
                    const T_DATA (&cell_volume),
                    
                    T_DATA (&U)[9],
                    T_DATA (&V)[9],
                    T_DATA (&Sigma)[3],
                    T_DATA (&dPdF)[12])
    {
        typedef Number<Tw> Tn;

        BUILD_TDATA(F,[9]);
        T_DATA alpha;
        T_DATA p;

        Weighted_Gradient<Tw,T_DATA,I_DATA>(u,F,weights,one_over_h);
        
        
        Tn val;
        /*
        Tn rone;
        rone.Load_Aligned(one);
        for(int i=0;i<9;i+=4){
            val.Load(F[i]);
            val = val + rone;
            Store(F[i],val);
        }
        */

        Singular_Value_Decomposition<Tw,T_DATA,I_DATA>(F,U,Sigma,V);

        Tn co;
        co.Load(cutoff);
        for(int i=0;i<3;i++){
            val.Load(Sigma[i]);
            val = max(val,co);
            Store(Sigma[i],val);
        }       

        Isotropic_Stress_Derivative<T_MATERIAL,Tw,T_DATA,I_DATA>::Run(dPdF,Sigma,p,mu,lambda,alpha,true);            
    }
    
};
                    


#define INSTANCE_KERNEL_Update_Position_Based_State(WIDTH)   \
    const WIDETYPE(float,WIDTH) (&u)[3][8],                  \
    const WIDETYPE(float,WIDTH) (&mu),                       \
    const WIDETYPE(float,WIDTH) (&lambda),                    \
    const WIDETYPE(float,WIDTH) (&weights)[3],                    \
    const WIDETYPE(float,WIDTH) (&cutoff),                   \
    const WIDETYPE(float,WIDTH) (&one_over_h)[3],               \
    const WIDETYPE(float,WIDTH) (&cell_volume),              \
    WIDETYPE(float,WIDTH) (&U)[9],                           \
    WIDETYPE(float,WIDTH) (&V)[9],                           \
    WIDETYPE(float,WIDTH) (&Sigma)[3],                       \
    WIDETYPE(float,WIDTH) (&dPdF)[12]                       



INSTANCE_KERNEL_MATERIAL(Update_Position_Based_State,COROTATED_TAG);
    //INSTANCE_KERNEL_MATERIAL(Update_Position_Based_State,NEOHOOKEAN_TAG);

#undef INSTANCE_KERNEL_Update_Position_Based_State
