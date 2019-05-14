//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_MATERIAL,class Tw,class T_DATA=void, class I_DATA=void>
struct Add_Force_Single_QPoint
{
    static void Run(const T_DATA (&mu),
                    const T_DATA (&lambda), 
                    const T_DATA (&U)[9],
                    const T_DATA (&V)[9],
                    const T_DATA (&Sigma)[3],
                    T_DATA (&P)[9]);


};


