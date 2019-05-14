//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_RAW,class T_DATA=void,class I_DATA=void>
void Weighted_Gradient(const T_DATA (&u)[3][8], T_DATA (&F)[9],const T_DATA (&W)[3], const T_DATA (&one_over_h)[3]);
