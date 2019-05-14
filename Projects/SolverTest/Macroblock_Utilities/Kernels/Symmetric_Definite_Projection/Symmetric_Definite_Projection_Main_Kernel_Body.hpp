//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifdef __INTEL_COMPILER
#pragma warning( disable : 592 )
#endif

#define USE_ACCURATE_RSQRT_IN_JACOBI_CONJUGATION
// #define PERFORM_STRICT_QUATERNION_RENORMALIZATION

{ // Begin block : Scope of Sdelta

Tn Vdelta1;
Tn Vdelta2;
Tn Vdelta3;

{ // Begin block : Scope of qV

Tn Vqvs;
Tn Vqvvx;
Tn Vqvvy;
Tn Vqvvz;

{ // Begin block : Scope of Stemp

Tn Vstemp11;
Tn Vstemp21;
Tn Vstemp31;
Tn Vstemp22;
Tn Vstemp32;
Tn Vstemp33;

Vstemp11=Vs11;
Vstemp21=Vs21;
Vstemp31=Vs31;
Vstemp22=Vs22;
Vstemp32=Vs32;
Vstemp33=Vs33;

{ // Begin block : Symmetric eigenanalysis

    Vqvs=Vone;
    Vqvvx = Vqvvx ^ Vqvvx;
    Vqvvy = Vqvvy ^ Vqvvy;
    Vqvvz = Vqvvz ^ Vqvvz;

    //###########################################################
    // Solve symmetric eigenproblem using Jacobi iteration
    //###########################################################

    for(int sweep=1;sweep<=4;sweep++){

        // First Jacobi conjugation

#define SS11 Sstemp11
#define SS21 Sstemp21
#define SS31 Sstemp31
#define SS22 Sstemp22
#define SS32 Sstemp32
#define SS33 Sstemp33
#define SQVVX Sqvvx
#define SQVVY Sqvvy
#define SQVVZ Sqvvz
#define STMP1 Stmp1
#define STMP2 Stmp2
#define STMP3 Stmp3

#define VS11 Vstemp11
#define VS21 Vstemp21
#define VS31 Vstemp31
#define VS22 Vstemp22
#define VS32 Vstemp32
#define VS33 Vstemp33
#define VQVVX Vqvvx
#define VQVVY Vqvvy
#define VQVVZ Vqvvz
#define VTMP1 Vtmp1
#define VTMP2 Vtmp2
#define VTMP3 Vtmp3

#include "Symmetric_Definite_Projection_Jacobi_Conjugation_Kernel.hpp"

#undef SS11
#undef SS21
#undef SS31
#undef SS22
#undef SS32
#undef SS33
#undef SQVVX
#undef SQVVY
#undef SQVVZ
#undef STMP1
#undef STMP2
#undef STMP3

#undef VS11
#undef VS21
#undef VS31
#undef VS22
#undef VS32
#undef VS33
#undef VQVVX
#undef VQVVY
#undef VQVVZ
#undef VTMP1
#undef VTMP2
#undef VTMP3

        // Second Jacobi conjugation

#define SS11 Sstemp22
#define SS21 Sstemp32
#define SS31 Sstemp21
#define SS22 Sstemp33
#define SS32 Sstemp31
#define SS33 Sstemp11
#define SQVVX Sqvvy
#define SQVVY Sqvvz
#define SQVVZ Sqvvx
#define STMP1 Stmp2
#define STMP2 Stmp3
#define STMP3 Stmp1

#define VS11 Vstemp22
#define VS21 Vstemp32
#define VS31 Vstemp21
#define VS22 Vstemp33
#define VS32 Vstemp31
#define VS33 Vstemp11
#define VQVVX Vqvvy
#define VQVVY Vqvvz
#define VQVVZ Vqvvx
#define VTMP1 Vtmp2
#define VTMP2 Vtmp3
#define VTMP3 Vtmp1

#include "Symmetric_Definite_Projection_Jacobi_Conjugation_Kernel.hpp"

#undef SS11
#undef SS21
#undef SS31
#undef SS22
#undef SS32
#undef SS33
#undef SQVVX
#undef SQVVY
#undef SQVVZ
#undef STMP1
#undef STMP2
#undef STMP3

#undef VS11
#undef VS21
#undef VS31
#undef VS22
#undef VS32
#undef VS33
#undef VQVVX
#undef VQVVY
#undef VQVVZ
#undef VTMP1
#undef VTMP2
#undef VTMP3

        // Third Jacobi conjugation

#define SS32 Sstemp21
#define SS33 Sstemp22
#define SQVVX Sqvvz
#define SQVVY Sqvvx
#define SQVVZ Sqvvy
#define STMP1 Stmp3
#define STMP2 Stmp1
#define STMP3 Stmp2

#define VS11 Vstemp33
#define VS21 Vstemp31
#define VS31 Vstemp32
#define VS22 Vstemp11
#define VS32 Vstemp21
#define VS33 Vstemp22
#define VQVVX Vqvvz
#define VQVVY Vqvvx
#define VQVVZ Vqvvy
#define VTMP1 Vtmp3
#define VTMP2 Vtmp1
#define VTMP3 Vtmp2

#include "Symmetric_Definite_Projection_Jacobi_Conjugation_Kernel.hpp"

#undef SS11
#undef SS21
#undef SS31
#undef SS22
#undef SS32
#undef SS33
#undef SQVVX
#undef SQVVY
#undef SQVVZ
#undef STMP1
#undef STMP2
#undef STMP3

#undef VS11
#undef VS21
#undef VS31
#undef VS22
#undef VS32
#undef VS33
#undef VQVVX
#undef VQVVY
#undef VQVVZ
#undef VTMP1
#undef VTMP2
#undef VTMP3
    }

} // End block : Symmetric eigenanalysis

Tn Vsabs21;
Tn Vsabs31;
Tn Vsabs32;

Vtmp1 = Vtmp1 ^ Vtmp1;

Vsabs21 = Vtmp1 - Vstemp21;
Vsabs21 = max(Vsabs21, Vstemp21);
Vsabs31 = Vtmp1 - Vstemp31;
Vsabs31 = max(Vsabs31, Vstemp31);
Vsabs32 = Vtmp1 - Vstemp32;
Vsabs32 = max(Vsabs32, Vstemp32);

Vdelta1 = Vsabs21 + Vsabs31;
Vdelta2 = Vsabs21 + Vsabs32;
Vdelta3 = Vsabs31 + Vsabs32;

Vdelta1 = max(Vstemp11, Vdelta1);
Vdelta2 = max(Vstemp22, Vdelta2);
Vdelta3 = max(Vstemp33, Vdelta3);

Vdelta1 = Vdelta1 - Vstemp11;
Vdelta2 = Vdelta2 - Vstemp22;
Vdelta3 = Vdelta3 - Vstemp33;

} // End block : Scope of Stemp

    //###########################################################
    // Normalize quaternion for matrix V
    //###########################################################

#if !defined(USE_ACCURATE_RSQRT_IN_JACOBI_CONJUGATION)
 || defined(PERFORM_STRICT_QUATERNION_RENORMALIZATION)

    Vtmp2 = Vqvs * Vqvs;
    Vtmp1 = Vqvvx * Vqvvx;
    Vtmp2 = Vtmp1 + Vtmp2;
    Vtmp1 = Vqvvy * Vqvvy;
    Vtmp2 = Vtmp1 + Vtmp2;
    Vtmp1 = Vqvvz * Vqvvz;
    Vtmp2 = Vtmp1 + Vtmp2;

    Vtmp1 = Vtmp2.rsqrt();
    Vtmp4 = Vtmp1 * Vone_half;
    Vtmp3 = Vtmp1 * Vtmp4;
    Vtmp3 = Vtmp1 * Vtmp3;
    Vtmp3 = Vtmp2 * Vtmp3;
    Vtmp1 = Vtmp1 + Vtmp4;
    Vtmp1 = Vtmp1 - Vtmp3;

    Vqvs = Vqvs * Vtmp1;
    Vqvvx = Vqvvx * Vtmp1;
    Vqvvy = Vqvvy * Vtmp1;
    Vqvvz = Vqvvz * Vtmp1;
#endif

{ // Begin block : Convert V to matrix

    //###########################################################
    // Transform quaternion to matrix V
    //###########################################################

#ifndef COMPUTE_V_AS_MATRIX
    Tn Vv11;
    Tn Vv21;
    Tn Vv31;
    Tn Vv12;
    Tn Vv22;
    Tn Vv32;
    Tn Vv13;
    Tn Vv23;
    Tn Vv33;
#endif

    Vtmp1 = Vqvvx * Vqvvx;
    Vtmp2 = Vqvvy * Vqvvy;
    Vtmp3 = Vqvvz * Vqvvz;
    Vv11 = Vqvs * Vqvs;
    Vv22 = Vv11 - Vtmp1;
    Vv33 = Vv22 - Vtmp2;
    Vv33 = Vv33 + Vtmp3;
    Vv22 = Vv22 + Vtmp2;
    Vv22 = Vv22 - Vtmp3;
    Vv11 = Vv11 + Vtmp1;
    Vv11 = Vv11 - Vtmp2;
    Vv11 = Vv11 - Vtmp3;
    Vtmp1 = Vqvvx + Vqvvx;
    Vtmp2 = Vqvvy + Vqvvy;
    Vtmp3 = Vqvvz + Vqvvz;
    Vv32 = Vqvs * Vtmp1;
    Vv13 = Vqvs * Vtmp2;
    Vv21 = Vqvs * Vtmp3;
    Vtmp1 = Vqvvy * Vtmp1;
    Vtmp2 = Vqvvz * Vtmp2;
    Vtmp3 = Vqvvx * Vtmp3;
    Vv12 = Vtmp1 - Vv21;
    Vv23 = Vtmp2 - Vv32;
    Vv31 = Vtmp3 - Vv13;
    Vv21 = Vtmp1 + Vv21;
    Vv32 = Vtmp2 + Vv32;
    Vv13 = Vtmp3 + Vv13;

} // End block : Convert V to matrix

} // End block : Scope of qV

// Project along first eigenvector

Vtmp1 = Vdelta1 * Vv11;
Vtmp2 = Vdelta1 * Vv21;
Vtmp3 = Vdelta1 * Vv31;
Vtmp4 = Vv11 * Vtmp1;
Vs11 = Vs11 + Vtmp4;
Vtmp4 = Vv21 * Vtmp1;
Vs21 = Vs21 + Vtmp4;
Vtmp4 = Vv31 * Vtmp1;
Vs31 = Vs31 + Vtmp4;
Vtmp4 = Vv21 * Vtmp2;
Vs22 = Vs22 + Vtmp4;
Vtmp4 = Vv31 * Vtmp2;
Vs32 = Vs32 + Vtmp4;
Vtmp4 = Vv31 * Vtmp3;
Vs33 = Vs33 + Vtmp4;

// Project along second eigenvector

Vtmp1 = Vdelta2 * Vv12;
Vtmp2 = Vdelta2 * Vv22;
Vtmp3 = Vdelta2 * Vv32;
Vtmp4 = Vv12 * Vtmp1;
Vs11 = Vs11 + Vtmp4;
Vtmp4 = Vv22 * Vtmp1;
Vs21 = Vs21 + Vtmp4;
Vtmp4 = Vv32 * Vtmp1;
Vs31 = Vs31 + Vtmp4;
Vtmp4 = Vv22 * Vtmp2;
Vs22 = Vs22 + Vtmp4;
Vtmp4 = Vv32 * Vtmp2;
Vs32 = Vs32 + Vtmp4;
Vtmp4 = Vv32 * Vtmp3;
Vs33 = Vs33 + Vtmp4;

// Project along third eigenvector

Vtmp1 = Vdelta3 * Vv13;
Vtmp2 = Vdelta3 * Vv23;
Vtmp3 = Vdelta3 * Vv33;
Vtmp4 = Vv13 * Vtmp1;
Vs11 = Vs11 + Vtmp4;
Vtmp4 = Vv23 * Vtmp1;
Vs21 = Vs21 + Vtmp4;
Vtmp4 = Vv33 * Vtmp1;
Vs31 = Vs31 + Vtmp4;
Vtmp4 = Vv23 * Vtmp2;
Vs22 = Vs22 + Vtmp4;
Vtmp4 = Vv33 * Vtmp2;
Vs32 = Vs32 + Vtmp4;
Vtmp4 = Vv33 * Vtmp3;
Vs33 = Vs33 + Vtmp4;

} // End block : Scope of delta

#ifdef __INTEL_COMPILER
#pragma warning( default : 592 )
#endif
