//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Weighted_Gradient
#include <assert.h>
#include "../../Common/KernelCommon.h"
#else
namespace {
#endif

template<class Tw,class T_DATA=void,class I_DATA=void>
#ifdef SUBROUTINE_Weighted_Gradient
inline
#endif
void Weighted_Gradient(const T_DATA (&u)[3][8], T_DATA (&F)[9],const T_DATA (&W)[3], const T_DATA (&one_over_h)[3])
{
    #define M_I(x,y) (y * 3 + x )

    typedef Number<Tw> Tn;

    Tn Ru000;  Tn Ru001;  Tn Ru010;  Tn Ru011;  
    Tn Ru100;  Tn Ru101;  Tn Ru110;  Tn Ru111;  

    // Interpolated along a single axes (located on edges)
    Tn RuI00;  Tn RuI01;  Tn RuI10;  Tn RuI11;  
    Tn Ru0I0;  Tn Ru0I1;  Tn Ru1I0;  Tn Ru1I1;  

    // Differentiated along one axis, interpolated along one (located on faces)
    Tn RuID0;  Tn RuID1;  
    Tn RuI0D;  Tn RuI1D;  
    Tn RuDI0;  Tn RuDI1;  

    // Differentiated along one axis, interpolated along two (located at cell interior)
    Tn RuDII;  Tn RuIDI;  Tn RuIID;  
	
    // Weight Values
    Tn w0;  
    Tn w1;  
    Tn w2;  
	
    // ONE_OVER_H
    Tn ONE_OVER_DX;
    Tn ONE_OVER_DY;
    Tn ONE_OVER_DZ;
    
    ONE_OVER_DX.Load(one_over_h[0]);
    ONE_OVER_DY.Load(one_over_h[1]);
    ONE_OVER_DZ.Load(one_over_h[2]);
    w0.Load(W[0]);
    w1.Load(W[1]);
    w2.Load(W[2]);

    // V = 0
    for( int v =0; v < 3; v++)
        {
            Ru000.Load(u[v][0]);
            Ru001.Load(u[v][1]);
            Ru010.Load(u[v][2]);
            Ru011.Load(u[v][3]);
            Ru100.Load(u[v][4]);
            Ru101.Load(u[v][5]);
            Ru110.Load(u[v][6]);
            Ru111.Load(u[v][7]);

            RuI00 = ((Ru100 - Ru000) * w0) + Ru000;
            RuI01 = ((Ru101 - Ru001) * w0) + Ru001;
            RuI10 = ((Ru110 - Ru010) * w0) + Ru010;
            RuI11 = ((Ru111 - Ru011) * w0) + Ru011;
            
            Ru0I0 = ((Ru010 - Ru000) * w1) + Ru000;
            Ru0I1 = ((Ru011 - Ru001) * w1) + Ru001;
            Ru1I0 = ((Ru110 - Ru100) * w1) + Ru100;
            Ru1I1 = ((Ru111 - Ru101) * w1) + Ru101;
            
            RuID0 = (RuI10 - RuI00) * ONE_OVER_DY;
            RuID1 = (RuI11 - RuI01) * ONE_OVER_DY;
            
            RuI0D = (RuI01 - RuI00) * ONE_OVER_DZ;
            RuI1D = (RuI11 - RuI10) * ONE_OVER_DZ;
            
            RuDI0 = (Ru1I0 - Ru0I0) * ONE_OVER_DX;
            RuDI1 = (Ru1I1 - Ru0I1) * ONE_OVER_DX;
            
            RuDII = ((RuDI1 - RuDI0) * w2) + RuDI0;
            RuIDI = ((RuID1 - RuID0) * w2) + RuID0;
            RuIID = ((RuI1D - RuI0D) * w1) + RuI0D;
            
            Store(F[M_I(v,0)], RuDII);
            Store(F[M_I(v,1)], RuIDI);
            Store(F[M_I(v,2)], RuIID);
        }
}

#ifdef SUBROUTINE_Weighted_Gradient
}
#else
#define INSTANCE_KERNEL_Weighted_Gradient(WIDTH) const WIDETYPE(float,WIDTH) (&u)[3][8], WIDETYPE(float,WIDTH) (&F)[9], const WIDETYPE(float,WIDTH) (&W)[3], const WIDETYPE(float,WIDTH) (&one_over_h)[3]
INSTANCE_KERNEL(Weighted_Gradient);
#undef INSTANCE_KERNEL_Weighted_Gradient
#endif
