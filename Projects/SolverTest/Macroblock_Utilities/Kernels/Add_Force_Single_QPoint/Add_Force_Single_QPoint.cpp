//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <assert.h>

#include "../../Common/KernelCommon.h"

#ifdef FORCE_INLINE
#include "../Matrix_Times_Transpose/Matrix_Times_Transpose.h"
#include "../Piola_Kirchhoff_Stress_Tensor/Piola_Kirchhoff_Stress_Tensor.h"
namespace {
    BUILD_CONSTANT(one,1.f);
}
#else


#define SUBROUTINE_Matrix_Times_Transpose
#include "../Matrix_Times_Transpose/Matrix_Times_Transpose.cpp"
#undef SUBROUTINE_Matrix_Times_Transpose

#define SUBROUTINE_Piola_Kirchhoff_Stress_Tensor
#include "../Piola_Kirchhoff_Stress_Tensor/Piola_Kirchhoff_Stress_Tensor.cpp"
#undef SUBROUTINE_Piola_Kirchhoff_Stress_Tensor

#endif

template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Add_Force_Single_QPoint
{
    static void Run(const T_DATA (&mu),
                    const T_DATA (&lambda),
                    const T_DATA (&U)[9],
                    const T_DATA (&V)[9],
                    const T_DATA (&Sigma)[3],
                    T_DATA (&P)[9])
    {
        typedef Number<Tw> Tn;
        typedef enum { x11=0,x21,x31,x12,x22,x32,x13,x23,x33} Matrix_Entry;
       
        T_DATA P_Hat[3];
        T_DATA Q_Hat[3];
        T_DATA p;
        T_DATA alpha;
        
        Piola_Kirchhoff_Stress_Tensor<T_MATERIAL,Tw,T_DATA,I_DATA>::Run(P_Hat, Sigma, 
                                                                        Q_Hat, p, mu, alpha, lambda);
        
        Tn rP;
        Store(P[x21], rP);
        Store(P[x31], rP);
        Store(P[x12], rP);
        Store(P[x32], rP);
        Store(P[x13], rP);
        Store(P[x23], rP);
        rP.Load(P_Hat[0]);
        Store(P[x11], rP);
        rP.Load(P_Hat[1]);
        Store(P[x22], rP);
        rP.Load(P_Hat[2]);
        Store(P[x33], rP);

        Matrix_Times_Transpose<Tw,T_DATA,I_DATA>(V,P,P);
        Matrix_Times_Transpose<Tw,T_DATA,I_DATA>(U,P,P);
    }

};


#define INSTANCE_KERNEL_Add_Force_Single_QPoint(WIDTH) \
                       const WIDETYPE(float,WIDTH) (&mu), \
                       const WIDETYPE(float,WIDTH) (&lambda), \
                       const WIDETYPE(float,WIDTH) (&U)[9], \
                       const WIDETYPE(float,WIDTH) (&V)[9], \
                       const WIDETYPE(float,WIDTH) (&Sigma)[3], \
                       WIDETYPE(float,WIDTH) (&P)[9]
INSTANCE_KERNEL_MATERIAL(Add_Force_Single_QPoint,COROTATED_TAG);
//INSTANCE_KERNEL_MATERIAL(Add_Force_Single_QPoint,NEOHOOKEAN_TAG);
//INSTANCE_KERNEL_MATERIAL(Add_Force_Single_QPoint,BIPHASIC_TAG);
#undef INSTANCE_KERNEL_Add_Force_Single_QPoint
