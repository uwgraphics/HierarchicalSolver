//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __DEVIATORIC_PIOLA_KIRCHHOFF_STRESS_TENSOR_H__
#define __DEVIATORIC_PIOLA_KIRCHHOFF_STRESS_TENSOR_H__
template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Deviatoric_Piola_Kirchhoff_Stress_Tensor
{
    static void Run(T_DATA (&DevP_Hat)[3], const T_DATA (&Sigma)[3], const T_DATA (&mu));
};
#endif
