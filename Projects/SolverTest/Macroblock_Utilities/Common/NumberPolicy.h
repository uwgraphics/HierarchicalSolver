//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################


#ifndef __NUMBER_POLICY_H__
#define __NUMBER_POLICY_H__

template<class Tw> class Number;
template<class Tw> class Discrete;
template<class Twi> class Mask;


//======================================================
//
//      NUMBER POLICY: Links Number to a Mask
//
//======================================================


struct ERROR_NO_MASK_TYPE;
struct ERROR_NO_DISCRETE_TYPE;
struct ERROR_NO_NUMBER_TYPE;

template<class NumberType> struct NumberPolicy{
    typedef Mask<ERROR_NO_MASK_TYPE> MASK_TYPE;
    typedef Discrete<ERROR_NO_DISCRETE_TYPE> DISCRETE_TYPE;
    typedef Number<ERROR_NO_NUMBER_TYPE> NUMBER_TYPE;
    static const int width=0;
};

template<> struct NumberPolicy<Number<float> >{
    typedef Mask<bool> MASK_TYPE;
    typedef Discrete<int> DISCRETE_TYPE;
    typedef Number<float> NUMBER_TYPE;
    static const int width=1;
};

// for double
template<> struct NumberPolicy<Number<double> >{
    typedef Mask<bool> MASK_TYPE;
    typedef Discrete<int> DISCRETE_TYPE;
    typedef Number<double> NUMBER_TYPE;
    static const int width=1;
};

#if defined(ENABLE_SSE_INSTRUCTION_SET)
template<> struct NumberPolicy<Number<__m128> >{
    typedef Mask<__m128> MASK_TYPE;
    typedef Discrete<__m128i> DISCRETE_TYPE;
    typedef Number<__m128> NUMBER_TYPE;
    static const int width=4;
};
#if defined(ENABLE_DOUBLE_SUPPORT)
template<> struct NumberPolicy<Number<__m128d> >{
    typedef Mask<__m128d> MASK_TYPE;
    typedef Discrete<__m128i> DISCRETE_TYPE;
    typedef Number<__m128d> NUMBER_TYPE;
    static const int width=2;
};
#endif
#endif

#if defined(ENABLE_AVX_INSTRUCTION_SET)
template<> struct NumberPolicy<Number<__m256> >{
    typedef Mask<__m256> MASK_TYPE;
    typedef Discrete<__m256i> DISCRETE_TYPE;
    typedef Number<__m256> NUMBER_TYPE;
    static const int width=8;
};
#if defined(ENABLE_DOUBLE_SUPPORT)
template<> struct NumberPolicy<Number<__m256d> >{
    typedef Mask<__m256d> MASK_TYPE;
    typedef Discrete<__m256i> DISCRETE_TYPE;
    typedef Number<__m256d> NUMBER_TYPE;
    static const int width=4;
};
#endif
#endif

#if defined(ENABLE_AVX512_INSTRUCTION_SET)
template<> struct NumberPolicy<Number<__m512> >{
    typedef Mask<__mmask16> MASK_TYPE;
    typedef Discrete<__m512i> DISCRETE_TYPE;
    typedef Number<__m512> NUMBER_TYPE;
    static const int width=16;
};
#if defined(ENABLE_DOUBLE_SUPPORT)
template<> struct NumberPolicy<Number<__m512d> >{
    typedef Mask<__mmask8> MASK_TYPE;
    typedef Discrete<__m512i> DISCRETE_TYPE;
    typedef Number<__m512d> NUMBER_TYPE;
    static const int width=8;
};
#endif
#endif

#if defined(ENABLE_MIC_INSTRUCTION_SET)
template<> struct NumberPolicy<Number<__m512> >{
    typedef Mask<__mmask16> MASK_TYPE;
    typedef Discrete<__m512i> DISCRETE_TYPE;
    typedef Number<__m512> NUMBER_TYPE;
    static const int width=16;
};
#endif

#if defined(ENABLE_NEON_INSTRUCTION_SET)
template<> struct NumberPolicy<Number<float32x4_t> >{
    typedef Mask<float32x4_t> MASK_TYPE;
    typedef Discrete<int32x4_t> DISCRETE_TYPE;
    typedef Number<float32x4_t> NUMBER_TYPE;
    static const int width=4;
};
#endif

//======================================================
//
//      DISCRETE POLICY: Links Discrete to a Mask
//
//      WARNING: FLOAT AND DOUBLE HAVE THE SAME DISCRETE POLICY.
//               DO NOT USE THIS TO LINK BACK TO DOUBLE
//
//======================================================

template<> struct NumberPolicy<Discrete<int> >{
    typedef Mask<bool> MASK_TYPE;
    typedef Discrete<int> DISCRETE_TYPE;
    typedef Number<float> NUMBER_TYPE;
    static const int width=1;
};

#if defined(ENABLE_SSE_INSTRUCTION_SET)
template<> struct NumberPolicy<Discrete<__m128i> >{
    typedef Mask<__m128> MASK_TYPE;
    typedef Discrete<__m128i> DISCRETE_TYPE;
    typedef Number<__m128> NUMBER_TYPE;
    static const int width=4;
};
#endif

#if defined(ENABLE_AVX_INSTRUCTION_SET)
template<> struct NumberPolicy<Discrete<__m256i> >{
    typedef Mask<__m256> MASK_TYPE;
    typedef Discrete<__m256i> DISCRETE_TYPE;
    typedef Number<__m256> NUMBER_TYPE;
    static const int width=8;
};
#endif

#if defined(ENABLE_AVX512_INSTRUCTION_SET)
template<> struct NumberPolicy<Discrete<__m512i> >{
    typedef Mask<__mmask16> MASK_TYPE;
    typedef Discrete<__m512i> DISCRETE_TYPE;
    typedef Number<__m512> NUMBER_TYPE;
    static const int width=16;
};
#endif

#if defined(ENABLE_MIC_INSTRUCTION_SET)
template<> struct NumberPolicy<Discrete<__m512i> >{
    typedef Mask<__mmask16> MASK_TYPE;
    typedef Discrete<__m512i> DISCRETE_TYPE;
    typedef Number<__m512> NUMBER_TYPE;
    static const int width=16;
};
#endif

#if defined(ENABLE_NEON_INSTRUCTION_SET)
template<> struct NumberPolicy<Discrete<int32x4_t> >{
    typedef Mask<float32x4_t> MASK_TYPE;
    typedef Discrete<int32x4_t> DISCRETE_TYPE;
    typedef Number<float32x4_t> NUMBER_TYPE;
    static const int width=4;
};
#endif

//======================================================
//
//      MASK POLICY: Links Mask to a common type
//
//======================================================

struct ERROR_NO_COMMON_TYPE;

template<class MaskType> struct MaskPolicy{
    typedef ERROR_NO_COMMON_TYPE MASK_EXTERNAL_TYPE;
    static const int width=0;
};

template<> struct MaskPolicy<Mask<bool> >{
    typedef float MASK_EXTERNAL_TYPE;
    static const int width=1;
};

#if defined(ENABLE_SSE_INSTRUCTION_SET)
template<> struct MaskPolicy<Mask<__m128> >{
    typedef float MASK_EXTERNAL_TYPE;
    static const int width=4;
};
#if defined(ENABLE_DOUBLE_SUPPORT)
template<> struct MaskPolicy<Mask<__m128d> >{
    typedef double MASK_EXTERNAL_TYPE;
    static const int width=2;
};
#endif
#endif

#if defined(ENABLE_AVX_INSTRUCTION_SET)
template<> struct MaskPolicy<Mask<__m256> >{
    typedef float MASK_EXTERNAL_TYPE;
    static const int width=8;
};
#if defined(ENABLE_DOUBLE_SUPPORT)
template<> struct MaskPolicy<Mask<__m256d> >{
    typedef double MASK_EXTERNAL_TYPE;
    static const int width=4;
};
#endif
#endif

#if defined(ENABLE_AVX512_INSTRUCTION_SET)
template<> struct MaskPolicy<Mask<__mmask16> >{
    typedef int MASK_EXTERNAL_TYPE;
    static const int width=16;
};
#if defined(ENABLE_DOUBLE_SUPPORT)
template<> struct MaskPolicy<Mask<__mmask8> >{
    typedef int MASK_EXTERNAL_TYPE;
    static const int width=8;
};
#endif
#endif

#if defined(ENABLE_SSE_INSTRUCTION_SET)
template<> struct MaskPolicy<Mask<__m128i> >{
    typedef int MASK_EXTERNAL_TYPE;
    static const int width=4;
};
#endif

#if defined(ENABLE_MIC_INSTRUCTION_SET)
template<> struct MaskPolicy<Mask<__mmask16> >{
    typedef int MASK_EXTERNAL_TYPE;
    static const int width=16;
};
#endif

#if defined(ENABLE_NEON_INSTRUCTION_SET)
template<> struct MaskPolicy<Mask<float32x4_t> >{
    typedef int MASK_EXTERNAL_TYPE;
    static const int width=4;
};
#endif


#endif
