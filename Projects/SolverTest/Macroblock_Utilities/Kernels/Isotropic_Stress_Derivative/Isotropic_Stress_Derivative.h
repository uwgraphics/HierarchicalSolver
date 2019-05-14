//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Isotropic_Stress_Derivative
{
static void Run(T_DATA (&dPdF)[12], const T_DATA (&Sigma)[3], const T_DATA (&p), const T_DATA (&mu), const T_DATA (&kappa), const T_DATA (&alpha), const bool apply_definiteness_fix);
};
