//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Rotated_Stress_Derivative
#include <assert.h>
#include "KernelCommon.h"
#else
namespace{
#endif

    namespace nm_Rotated_Stress_Derivative {
        BUILD_CONSTANT(three,3.0f);
    }

    using namespace nm_Rotated_Stress_Derivative;

#ifndef __Rotated_Stress_Derivative_defined__
#define __Rotated_Stress_Derivative_defined__
template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Rotated_Stress_Derivative
{
static inline void Run(T_DATA (&dPdF)[12], const T_DATA (&Sigma)[3], const T_DATA (&mu), const T_DATA (&kappa));
};
#endif

template<class Tw,class T_DATA,class I_DATA>
struct Rotated_Stress_Derivative<COROTATED_TAG,Tw,T_DATA,I_DATA>
{
static inline void Run(T_DATA (&dPdF)[12], const T_DATA (&Sigma)[3], const T_DATA (&mu), const T_DATA (&kappa))
{
    typedef enum { x1111=0,x1122,x1133,x2222,x2233,x3333,x1212,x1221,x1313,x1331,x2323,x2332 } RSD_Entry;

    typedef Number<Tw> Tn;

    Tn rthree;
    rthree.Load_Aligned(three);

    Tn rmu,rkappa;
    rmu.Load(mu);rkappa.Load(kappa);

    Tn rtwo_mu=rmu+rmu;

    Tn rSigma1;
    Tn rSigma2;
    Tn rSigma3;

    rSigma1.Load(Sigma[0]);
    rSigma2.Load(Sigma[1]);
    rSigma3.Load(Sigma[2]);

    Tn rtwo_mu_minus_kappa_strain_trace=rtwo_mu-rkappa*(rSigma1+rSigma2+rSigma3-rthree);

    Tn ra1111=rtwo_mu+rkappa;
    Tn ra1122=rkappa;
    Tn ra1133=rkappa;
    Tn ra2222=rtwo_mu+rkappa;
    Tn ra2233=rkappa;
    Tn ra3333=rtwo_mu+rkappa;

    Tn ra1212=rtwo_mu-rtwo_mu_minus_kappa_strain_trace*(rSigma1+rSigma2).inverse();
    Tn ra1221=rtwo_mu_minus_kappa_strain_trace*(rSigma1+rSigma2).inverse();     
    Tn ra1313=rtwo_mu-rtwo_mu_minus_kappa_strain_trace*(rSigma1+rSigma3).inverse();
    Tn ra1331=rtwo_mu_minus_kappa_strain_trace*(rSigma1+rSigma3).inverse();     
    Tn ra2323=rtwo_mu-rtwo_mu_minus_kappa_strain_trace*(rSigma2+rSigma3).inverse();
    Tn ra2332=rtwo_mu_minus_kappa_strain_trace*(rSigma2+rSigma3).inverse();

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
struct Rotated_Stress_Derivative<NEOHOOKEAN_TAG,Tw,T_DATA,I_DATA>
{
static inline void Run(T_DATA (&dPdF)[12], const T_DATA (&Sigma)[3], const T_DATA (&mu), const T_DATA (&kappa))
{
    typedef enum { x1111=0,x1122,x1133,x2222,x2233,x3333,x1212,x1221,x1313,x1331,x2323,x2332 } RSD_Entry;

    typedef Number<Tw> Tn;

    Tn rmu,rkappa;
    rmu.Load(mu);rkappa.Load(kappa);

    Tn rSigma1_inverse;
    Tn rSigma2_inverse;
    Tn rSigma3_inverse;

    rSigma1_inverse.Load(Sigma[0]);
    rSigma2_inverse.Load(Sigma[1]);
    rSigma3_inverse.Load(Sigma[2]);

    Tn rlogJ=(rSigma1_inverse*rSigma2_inverse*rSigma3_inverse).log();

    rSigma1_inverse=rSigma1_inverse.inverse();
    rSigma2_inverse=rSigma2_inverse.inverse();
    rSigma3_inverse=rSigma3_inverse.inverse();

    Tn ralpha11=rmu;
    Tn ralpha12=rmu;
    Tn ralpha13=rmu;
    Tn ralpha22=rmu;
    Tn ralpha23=rmu;
    Tn ralpha33=rmu;

    Tn rmu_minus_kappa_logJ=rmu-rkappa*rlogJ;
    Tn rbeta11=rmu_minus_kappa_logJ*rSigma1_inverse*rSigma1_inverse;
    Tn rbeta12=rmu_minus_kappa_logJ*rSigma1_inverse*rSigma2_inverse;
    Tn rbeta13=rmu_minus_kappa_logJ*rSigma1_inverse*rSigma3_inverse;
    Tn rbeta22=rmu_minus_kappa_logJ*rSigma2_inverse*rSigma2_inverse;
    Tn rbeta23=rmu_minus_kappa_logJ*rSigma2_inverse*rSigma3_inverse;
    Tn rbeta33=rmu_minus_kappa_logJ*rSigma3_inverse*rSigma3_inverse;

    Tn rgamma11=rkappa*rSigma1_inverse*rSigma1_inverse;
    Tn rgamma12=rkappa*rSigma1_inverse*rSigma2_inverse;
    Tn rgamma13=rkappa*rSigma1_inverse*rSigma3_inverse;
    Tn rgamma22=rkappa*rSigma2_inverse*rSigma2_inverse;
    Tn rgamma23=rkappa*rSigma2_inverse*rSigma3_inverse;
    Tn rgamma33=rkappa*rSigma3_inverse*rSigma3_inverse;
    
    Tn ra1111=ralpha11+rbeta11+rgamma11;
    Tn ra1122=rgamma12;
    Tn ra1133=rgamma13;
    Tn ra2222=ralpha22+rbeta22+rgamma22;
    Tn ra2233=rgamma23;
    Tn ra3333=ralpha33+rbeta33+rgamma33;

    Tn ra1212=ralpha12;
    Tn ra1221=rbeta12;
    Tn ra1313=ralpha13;
    Tn ra1331=rbeta13;
    Tn ra2323=ralpha23;
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
struct Rotated_Stress_Derivative<BIPHASIC_TAG,Tw,T_DATA,I_DATA>
{
static inline void Run(T_DATA (&dPdF)[12], const T_DATA (&Sigma)[3], const T_DATA (&mu), const T_DATA (&kappa))
{
    typedef enum { x1111=0,x1122,x1133,x2222,x2233,x3333,x1212,x1221,x1313,x1331,x2323,x2332 } RSD_Entry;

    typedef Number<Tw> Tn;

    Tn rZero;
    Tn rthree;
    rthree.Load_Aligned(three);    

    Tn rBiphasicThreshold;
    rBiphasicThreshold.Load(nm_biphasic::biphasic_threshold);
    Tn rKbp;
    rKbp.Load(nm_biphasic::biphasic_factor);

    Tn rmu,rkappa;
    rmu.Load(mu);rkappa.Load(kappa);

    Tn rtwo_mu=rmu+rmu;

    Tn rSigma1;
    Tn rSigma2;
    Tn rSigma3;

    rSigma1.Load(Sigma[0]);
    rSigma2.Load(Sigma[1]);
    rSigma3.Load(Sigma[2]);

    Tn rtwo_mu_minus_kappa_strain_trace=rtwo_mu-rkappa*(rSigma1+rSigma2+rSigma3-rthree);

    Tn ra1111=rtwo_mu+rkappa + blend(rSigma1 > rBiphasicThreshold, rZero, rKbp);
    Tn ra1122=rkappa;
    Tn ra1133=rkappa;
    Tn ra2222=rtwo_mu+rkappa + blend(rSigma2 > rBiphasicThreshold, rZero, rKbp);
    Tn ra2233=rkappa;
    Tn ra3333=rtwo_mu+rkappa + blend(rSigma3 > rBiphasicThreshold, rZero, rKbp);

    Tn ra1212=rtwo_mu-rtwo_mu_minus_kappa_strain_trace*(rSigma1+rSigma2).inverse();
    Tn ra1221=rtwo_mu_minus_kappa_strain_trace*(rSigma1+rSigma2).inverse();     
    Tn ra1313=rtwo_mu-rtwo_mu_minus_kappa_strain_trace*(rSigma1+rSigma3).inverse();
    Tn ra1331=rtwo_mu_minus_kappa_strain_trace*(rSigma1+rSigma3).inverse();     
    Tn ra2323=rtwo_mu-rtwo_mu_minus_kappa_strain_trace*(rSigma2+rSigma3).inverse();
    Tn ra2332=rtwo_mu_minus_kappa_strain_trace*(rSigma2+rSigma3).inverse();

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

#ifndef SUBROUTINE_Rotated_Stress_Derivative
#define INSTANCE_KERNEL_Rotated_Stress_Derivative(WIDTH) WIDETYPE(float,WIDTH) (&dPdF)[12],    const WIDETYPE(float,WIDTH) (&Sigma)[3],    const WIDETYPE(float,WIDTH) (&mu),    const WIDETYPE(float,WIDTH) (&kappa)
INSTANCE_KERNEL_MATERIAL(Rotated_Stress_Derivative,COROTATED_TAG);
INSTANCE_KERNEL_MATERIAL(Rotated_Stress_Derivative,NEOHOOKEAN_TAG);
INSTANCE_KERNEL_MATERIAL(Rotated_Stress_Derivative,BIPHASIC_TAG);
#undef INSTANCE_KERNEL_Rotated_Stress_Derivative
#else
}
#endif
