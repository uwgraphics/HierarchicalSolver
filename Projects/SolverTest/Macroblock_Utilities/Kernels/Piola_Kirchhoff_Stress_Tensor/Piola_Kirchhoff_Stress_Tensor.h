//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Piola_Kirchhoff_Stress_Tensor
{
    static void Run(T_DATA (&P_Hat)[3], const T_DATA (&Sigma)[3], const T_DATA (&Q_Hat)[3], const T_DATA (&p), const T_DATA (&mu), const T_DATA (&alpha), const T_DATA (&kappa) );
};
