//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Augmented_Rotated_Stress_Derivative
#include <assert.h>
#include "KernelCommon.h"
#else
namespace{
#endif

#ifndef __Augmented_Rotated_Stress_Derivative_defined__
#define __Augmented_Rotated_Stress_Derivative_defined__
template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Augmented_Rotated_Stress_Derivative
{
static inline void Run(T_DATA (&dPdF)[12], const T_DATA (&Sigma)[3], const T_DATA (&p), const T_DATA (&mu), const T_DATA (&alpha));
};
#endif

namespace {
 }

template<class Tw,class T_DATA,class I_DATA>
struct Augmented_Rotated_Stress_Derivative<COROTATED_TAG,Tw,T_DATA,I_DATA>
{
static inline void Run(T_DATA (&dPdF)[12], const T_DATA (&Sigma)[3], const T_DATA (&p), const T_DATA (&mu), const T_DATA (&alpha))
{
    typedef enum { x1111=0,x1122,x1133,x2222,x2233,x3333,x1212,x1221,x1313,x1331,x2323,x2332 } RSD_Entry;

    typedef Number<Tw> Tn;

    Tn rp,rmu,ralpha;
    rp.Load(p);rmu.Load(mu);ralpha.Load(alpha);    

    Tn rtwo_mu=rmu+rmu;
    Tn rtwo_mu_minus_alpha_p=rtwo_mu-ralpha*rp;

    Tn rSigma1;
    Tn rSigma2;
    Tn rSigma3;

    rSigma1.Load(Sigma[0]);
    rSigma2.Load(Sigma[1]);
    rSigma3.Load(Sigma[2]);

    Tn ra1111=rtwo_mu;
    Tn ra1122;
    Tn ra1133;
    Tn ra2222=rtwo_mu;
    Tn ra2233;
    Tn ra3333=rtwo_mu;

    Tn ra1212=rtwo_mu-rtwo_mu_minus_alpha_p*(rSigma1+rSigma2).inverse();
    Tn ra1221=rtwo_mu_minus_alpha_p*(rSigma1+rSigma2).inverse();     
    Tn ra1313=rtwo_mu-rtwo_mu_minus_alpha_p*(rSigma1+rSigma3).inverse();
    Tn ra1331=rtwo_mu_minus_alpha_p*(rSigma1+rSigma3).inverse();     
    Tn ra2323=rtwo_mu-rtwo_mu_minus_alpha_p*(rSigma2+rSigma3).inverse();
    Tn ra2332=rtwo_mu_minus_alpha_p*(rSigma2+rSigma3).inverse();     
        
    Store(dPdF[x1111],ra1111);
    Store(dPdF[x1122],ra1122);
    Store(dPdF[x1133],ra1133);
    Store(dPdF[x2222],ra2222);
    Store(dPdF[x2233],ra2233);
    Store(dPdF[x3333],ra3333);
    Store(dPdF[x1212],ra1212);
    Store(dPdF[x1221],ra1221);
    Store(dPdF[x1313],ra1313);
    Store(dPdF[x1331],ra1331);
    Store(dPdF[x2323],ra2323);
    Store(dPdF[x2332],ra2332);
}
};

template<class Tw,class T_DATA,class I_DATA>
struct Augmented_Rotated_Stress_Derivative<NEOHOOKEAN_TAG,Tw,T_DATA,I_DATA>
{
static inline void Run(T_DATA (&dPdF)[12], const T_DATA (&Sigma)[3], const T_DATA (&p), const T_DATA (&mu), const T_DATA (&alpha))
{
    typedef enum { x1111=0,x1122,x1133,x2222,x2233,x3333,x1212,x1221,x1313,x1331,x2323,x2332 } RSD_Entry;

    typedef Number<Tw> Tn;

    Tn rp,rmu,ralpha;
    rp.Load(p);rmu.Load(mu);ralpha.Load(alpha);    

    Tn rmu_minus_alpha_p=rmu-ralpha*rp;

    Tn rSigma1_inverse;
    Tn rSigma2_inverse;
    Tn rSigma3_inverse;

    rSigma1_inverse.Load(Sigma[0]);
    rSigma2_inverse.Load(Sigma[1]);
    rSigma3_inverse.Load(Sigma[2]);

    rSigma1_inverse=rSigma1_inverse.inverse();
    rSigma2_inverse=rSigma2_inverse.inverse();
    rSigma3_inverse=rSigma3_inverse.inverse();

    Tn rbeta11=rmu_minus_alpha_p*rSigma1_inverse*rSigma1_inverse;
    Tn rbeta12=rmu_minus_alpha_p*rSigma1_inverse*rSigma2_inverse;
    Tn rbeta13=rmu_minus_alpha_p*rSigma1_inverse*rSigma3_inverse;
    Tn rbeta22=rmu_minus_alpha_p*rSigma2_inverse*rSigma2_inverse;
    Tn rbeta23=rmu_minus_alpha_p*rSigma2_inverse*rSigma3_inverse;
    Tn rbeta33=rmu_minus_alpha_p*rSigma3_inverse*rSigma3_inverse;

    Tn ra1111=rmu+rbeta11;
    Tn ra1122;
    Tn ra1133;
    Tn ra2222=rmu+rbeta22;
    Tn ra2233;
    Tn ra3333=rmu+rbeta33;

    Tn ra1212=rmu;
    Tn ra1221=rbeta12;
    Tn ra1313=rmu;
    Tn ra1331=rbeta13;
    Tn ra2323=rmu;
    Tn ra2332=rbeta23;
        
    Store(dPdF[x1111],ra1111);
    Store(dPdF[x1122],ra1122);
    Store(dPdF[x1133],ra1133);
    Store(dPdF[x2222],ra2222);
    Store(dPdF[x2233],ra2233);
    Store(dPdF[x3333],ra3333);
    Store(dPdF[x1212],ra1212);
    Store(dPdF[x1221],ra1221);
    Store(dPdF[x1313],ra1313);
    Store(dPdF[x1331],ra1331);
    Store(dPdF[x2323],ra2323);
    Store(dPdF[x2332],ra2332);
}
};

template<class Tw,class T_DATA,class I_DATA>
struct Augmented_Rotated_Stress_Derivative<BIPHASIC_TAG,Tw,T_DATA,I_DATA>
{
static inline void Run(T_DATA (&dPdF)[12], const T_DATA (&Sigma)[3], const T_DATA (&p), const T_DATA (&mu), const T_DATA (&alpha))
{
    typedef enum { x1111=0,x1122,x1133,x2222,x2233,x3333,x1212,x1221,x1313,x1331,x2323,x2332 } RSD_Entry;

    typedef Number<Tw> Tn;

    Tn rZero;
    Tn rp,rmu,ralpha;
    Tn rBiphasicThreshold;
    rp.Load(p);rmu.Load(mu);ralpha.Load(alpha);    
    rBiphasicThreshold.Load(nm_biphasic::biphasic_threshold);
    Tn rKbp;
    rKbp.Load(nm_biphasic::biphasic_factor);

    Tn rtwo_mu=rmu+rmu;
    Tn rtwo_mu_minus_alpha_p=rtwo_mu-ralpha*rp;

    Tn rSigma1;
    Tn rSigma2;
    Tn rSigma3;

    rSigma1.Load(Sigma[0]);
    rSigma2.Load(Sigma[1]);
    rSigma3.Load(Sigma[2]);

    Tn ra1111=rtwo_mu + blend(rSigma1 > rBiphasicThreshold, rZero, rKbp);
    Tn ra1122;
    Tn ra1133;
    Tn ra2222=rtwo_mu + blend(rSigma2 > rBiphasicThreshold, rZero, rKbp);
    Tn ra2233;
    Tn ra3333=rtwo_mu + blend(rSigma3 > rBiphasicThreshold, rZero, rKbp);

    Tn ra1212=rtwo_mu-rtwo_mu_minus_alpha_p*(rSigma1+rSigma2).inverse();
    Tn ra1221=rtwo_mu_minus_alpha_p*(rSigma1+rSigma2).inverse();     
    Tn ra1313=rtwo_mu-rtwo_mu_minus_alpha_p*(rSigma1+rSigma3).inverse();
    Tn ra1331=rtwo_mu_minus_alpha_p*(rSigma1+rSigma3).inverse();     
    Tn ra2323=rtwo_mu-rtwo_mu_minus_alpha_p*(rSigma2+rSigma3).inverse();
    Tn ra2332=rtwo_mu_minus_alpha_p*(rSigma2+rSigma3).inverse();     
        
    Store(dPdF[x1111],ra1111);
    Store(dPdF[x1122],ra1122);
    Store(dPdF[x1133],ra1133);
    Store(dPdF[x2222],ra2222);
    Store(dPdF[x2233],ra2233);
    Store(dPdF[x3333],ra3333);
    Store(dPdF[x1212],ra1212);
    Store(dPdF[x1221],ra1221);
    Store(dPdF[x1313],ra1313);
    Store(dPdF[x1331],ra1331);
    Store(dPdF[x2323],ra2323);
    Store(dPdF[x2332],ra2332);
}
};

#ifndef SUBROUTINE_Augmented_Rotated_Stress_Derivative
#define INSTANCE_KERNEL_Augmented_Rotated_Stress_Derivative(WIDTH) WIDETYPE(float,WIDTH) (&dPdF)[12],    const WIDETYPE(float,WIDTH) (&Sigma)[3],    const WIDETYPE(float,WIDTH) (&p),    const WIDETYPE(float,WIDTH) (&mu),    const WIDETYPE(float,WIDTH) (&alpha)
INSTANCE_KERNEL_MATERIAL(Augmented_Rotated_Stress_Derivative, COROTATED_TAG);
INSTANCE_KERNEL_MATERIAL(Augmented_Rotated_Stress_Derivative, NEOHOOKEAN_TAG);
INSTANCE_KERNEL_MATERIAL(Augmented_Rotated_Stress_Derivative, BIPHASIC_TAG);
#undef INSTANCE_KERNEL_Augmented_Rotated_Stress_Derivative
#else
}
#endif
