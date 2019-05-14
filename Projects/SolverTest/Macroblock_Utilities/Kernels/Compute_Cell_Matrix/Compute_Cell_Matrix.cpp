//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################

#ifndef SUBROUTINE_Compute_Cell_Matrix
#include <assert.h>
#include "../../Common/KernelCommon.h"
#endif


#ifdef FORCE_INLINE
#include "../Weighted_Accumulation/Weighted_Accumulation.cpp"
#else
#define SUBROUTINE_Unweighted_Accumulation
#include "../Weighted_Accumulation/Weighted_Accumulation.cpp"
#undef SUBROUTINE_Unweighted_Accumulation
#endif


template<class Tw,class T_DATA=void, class I_DATA=void>
void Transpose(T_DATA (&At)[9], const T_DATA (&A)[9])
{
    typedef enum { x11=0,x21,x31,x12,x22,x32,x13,x23,x33} Matrix_Entry;
    typedef enum { t11=0,t21=3,t31=6,t12=1,t22=4,t32=7,t13=2,t23=5,t33=8} TransMatrix_Entry;
    typedef Number<Tw> Tn;
    
    Tn Entry;
    Entry.Load(A[x11]);
    Store(At[t11],Entry);
    
    Entry.Load(A[x21]);
    Store(At[t21],Entry);
    
    Entry.Load(A[x31]);
    Store(At[t31],Entry);
    
    Entry.Load(A[x12]);
    Store(At[t12],Entry);
    
    Entry.Load(A[x22]);
    Store(At[t22],Entry);
    
    Entry.Load(A[x32]);
    Store(At[t32],Entry);
    
    Entry.Load(A[x13]);
    Store(At[t13],Entry);
    
    Entry.Load(A[x23]);
    Store(At[t23],Entry);
    
    Entry.Load(A[x33]);
    Store(At[t33],Entry);
}

template<class Tw,class T_DATA=void, class I_DATA=void>
void Build_M(T_DATA (&M)[9], const T_DATA (&dPdF)[12], const T_DATA (&H)[3][8], const T_DATA (&scale), int i, int k )
{
    
    typedef enum { x1111=0,x1122,x1133,x2222,x2233,x3333,x1212,x1221,x1313,x1331,x2323,x2332 } RSD_Entry;
    typedef enum { x11=0,x21,x31,x12,x22,x32,x13,x23,x33} Matrix_Entry;
    typedef Number<Tw> Tn;
    
    Tn ra1111, ra1122, ra1133, ra2222, ra2233, ra3333, ra1212, ra1221, ra1313, ra1331, ra2323, ra2332;
    Tn rM[9];
    Tn rhi[3];
    Tn rhk[3];
    Tn rscale;
    
    ra1111.Load(dPdF[x1111]);
    ra1122.Load(dPdF[x1122]);
    ra1133.Load(dPdF[x1133]);
    ra2222.Load(dPdF[x2222]);
    ra2233.Load(dPdF[x2233]);
    ra3333.Load(dPdF[x3333]);
    ra1212.Load(dPdF[x1212]);
    ra1221.Load(dPdF[x1221]);
    ra1313.Load(dPdF[x1313]);
    ra1331.Load(dPdF[x1331]);
    ra2323.Load(dPdF[x2323]);
    ra2332.Load(dPdF[x2332]);
    
    rhi[0].Load(H[0][i]);
    rhi[1].Load(H[1][i]);
    rhi[2].Load(H[2][i]);
    
    rhk[0].Load(H[0][k]);
    rhk[1].Load(H[1][k]);
    rhk[2].Load(H[2][k]);
    
    rscale.Load(scale);
    
    /* What goes here??? */
    rM[0] = rscale*(rhi[0]*rhk[0]*ra1111 + rhi[1]*rhk[1]*ra1212 + rhi[2]*rhk[2]*ra1313);
    rM[1] = rscale*(rhi[0]*rhk[0]*ra1212 + rhi[1]*rhk[1]*ra2222 + rhi[2]*rhk[2]*ra2323); 
    rM[2] = rscale*(rhi[0]*rhk[0]*ra1313 + rhi[1]*rhk[1]*ra2323 + rhi[2]*rhk[2]*ra3333);       
    rM[3] = rscale*(rhi[0]*rhk[1]*ra1122 + rhi[1]*rhk[0]*ra1221);
    rM[4] = rscale*(rhi[0]*rhk[1]*ra1221 + rhi[1]*rhk[0]*ra1122);
    rM[5] = rscale*(rhi[0]*rhk[2]*ra1133 + rhi[2]*rhk[0]*ra1331);        
    rM[6] = rscale*(rhi[0]*rhk[2]*ra1331 + rhi[2]*rhk[0]*ra1133);
    rM[7] = rscale*(rhi[2]*rhk[1]*ra2233 + rhi[1]*rhk[2]*ra2332);
    rM[8] = rscale*(rhi[2]*rhk[1]*ra2332 + rhi[1]*rhk[2]*ra2233);
    
    Store(M[x11], rM[0]);
    Store(M[x22], rM[1]);
    Store(M[x33], rM[2]);        
    Store(M[x12], rM[3]);
    Store(M[x21], rM[4]);
    Store(M[x13], rM[5]);        
    Store(M[x31], rM[6]);
    Store(M[x32], rM[7]);
    Store(M[x23], rM[8]);
}



template<class Tw,class T_DATA=void, class I_DATA=void>
#ifdef SUBROUTINE_Compute_Cell_Matrix
inline
#endif
    void Compute_Cell_Matrix(const int v1, const int v2,
                             const T_DATA (&one_over_h),
                             const T_DATA (&cell_volume),
                             const T_DATA (&Weights)[3],
                             const T_DATA (&U)[9],
                             const T_DATA (&V)[9],
                             const T_DATA (&dPdF)[12],
                             T_DATA (&ij_matrix)[9])
{
    #if 0
    typedef Number<Tw> Tn;
    typedef enum { x11=0,x21,x31,x12,x22,x32,x13,x23,x33} Matrix_Entry;

    Tn Zero;
    //Tn rSev;
    //Tn rmustab;
    Tn rcV;

    //Tn rminusseven;
    //Tn rone;
    //Tn rthree;
    //Tn rminusfive;
    BUILD_TDATA(_one, );

    //rSev.Load_Aligned( CDC::seven );
    //rminusseven = Tn() - rSev;
    //rmustab.Load(mu_stab);
    rcV.Load(cell_volume);
    //rone.Load_Aligned( CDC::one );
    //Store(_one, rone);
    //rthree.Load_Aligned(CDC::three);
    //rminusfive.Load_Aligned(CDC::minusfive);

    // K_Stab Matrix // This should work, even with wacky Number classes
    //const Tn* const K_Stab[8][8] = {
    //    /*0*/{&rminusseven,  &rthree,  &rthree,  &rone,  &rthree,  &rone,  &rone, &rminusfive},
    //    /*1*/{ &rthree, &rminusseven,  &rone,  &rthree,  &rone,  &rthree, &rminusfive,  &rone},
    //    /*2*/{ &rthree,  &rone, &rminusseven,  &rthree,  &rone, &rminusfive,  &rthree,  &rone},
    //    /*3*/{ &rone,  &rthree,  &rthree, &rminusseven, &rminusfive,  &rone,  &rone,  &rthree},
    //    /*4*/{ &rthree,  &rone,  &rone, &rminusfive, &rminusseven,  &rthree,  &rthree,  &rone},
    //    /*5*/{ &rone,  &rthree, &rminusfive,  &rone,  &rthree, &rminusseven,  &rone,  &rthree},
    //    /*6*/{ &rone, &rminusfive,  &rthree,  &rone,  &rthree,  &rone, &rminusseven,  &rthree},
    //    /*7*/{&rminusfive,  &rone,  &rone,  &rthree,  &rone,  &rthree,  &rthree, &rminusseven}
    //};


    // BUILD H
    BUILD_TDATA(H,[3][8]);
    for( int k=0; k<3; k++)
        for( int l=0; l<8; l++)
            Store(H[k][l], Zero);
    
    BUILD_TDATA(Vt,[9]);
    Transpose<Tw,T_DATA,I_DATA>(Vt, V);
    Unweighted_Accumulation<Tw,T_DATA,I_DATA>(H,Vt,one_over_h,_one);
    // -------

    for( int i = 0; i < 8; i++){
        for( int k = 0; k < 8; k++){
            // Build M here ???
            Build_M<Tw,T_DATA,I_DATA>(Vt, dPdF, H, i, k);
            for( int j=0; j<3; j++)
                for( int l=0; l<3; l++){
                    Tn df_du_value;
                    for( int p=0; p<3; p++)
                        for( int q=0; q<3; q++){
                            Tn rM_pq, U_jp, U_lq;
                            rM_pq.Load( Vt[p + q*3] );
                            U_jp.Load( U[j + p*3] );
                            U_lq.Load( U[l + q*3] );
                            df_du_value = df_du_value + (rM_pq*U_jp*U_lq);
                        }
                    df_du_value = Tn() - ((df_du_value * rcV) + ((*K_Stab[i][k]) * rmustab));

                    if( !((j>l) || (j==l && k < i )) ){
                        int I = i+j*8;
                        int J = k+l*8;
                        int I_SUM = 0;
                        for( int _j=0; _j<=I; I_SUM+=_j, _j++);
                        int flat_index = J - I_SUM + I*24;
                        Store( matrix[flat_index], df_du_value );
                    }
                }

        }
    }
    #endif
}

#ifndef SUBROUTINE_Compute_Cell_Matrix

#define INSTANCE_KERNEL_Transpose(WIDTH) WIDETYPE(float,WIDTH) (&At)[9], const WIDETYPE(float,WIDTH) (&A)[9]
INSTANCE_KERNEL(Transpose);
#undef INSTANCE_KERNEL_Transpose

#define INSTANCE_KERNEL_Build_M(WIDTH) WIDETYPE(float,WIDTH) (&M)[9], const WIDETYPE(float,WIDTH) (&dPdF)[12], const WIDETYPE(float,WIDTH) (&H)[3][8], const WIDETYPE(float,WIDTH) (&scale), int i, int k
INSTANCE_KERNEL(Build_M);
#undef INSTANCE_KERNEL_Build_M

#define INSTANCE_KERNEL_Compute_Cell_Matrix(WIDTH) const int v1, const int v2, const WIDETYPE(float,WIDTH) (&one_over_h), const WIDETYPE(float,WIDTH) (&cell_volume), const WIDETYPE(float,WIDTH) (&Weights)[3], const WIDETYPE(float,WIDTH) (&U)[9], const WIDETYPE(float,WIDTH) (&V)[9], const WIDETYPE(float,WIDTH) (&dPdF)[12], WIDETYPE(float,WIDTH) (&ij_matrix)[9]
INSTANCE_KERNEL(Compute_Cell_Matrix);
#undef INSTANCE_KERNEL_Compute_Cell_Matrix
#endif
