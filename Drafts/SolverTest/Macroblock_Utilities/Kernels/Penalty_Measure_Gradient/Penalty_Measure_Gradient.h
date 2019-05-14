//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Penalty_Measure_Gradient
{
static void Run(const T_DATA (&Sigma)[3], T_DATA (&Q_hat)[3]);
};
