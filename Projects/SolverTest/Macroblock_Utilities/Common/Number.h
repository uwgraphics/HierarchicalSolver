//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_NUMBER_H__
#define __KERNEL_NUMBER_H__

//#include <math.h>

#ifdef ENABLE_IO_SUPPORT
#include <iostream>
#endif


namespace {
    typedef union {
        int i;
        float f;
    } floatConverter;

    template<typename Tw> struct ScalarType;
    template<> struct ScalarType<float>{ using type=float; };
    template<> struct ScalarType<double>{ using type=double; };
#if defined(ENABLE_SSE_INSTRUCTION_SET)
    template<> struct ScalarType<__m128>{ using type=float; };
#if defined(ENABLE_DOUBLE_SUPPORT)
    template<> struct ScalarType<__m128d>{ using type=double; };
#endif
#endif

#if defined(ENABLE_AVX_INSTRUCTION_SET)
    template<> struct ScalarType<__m256>{ using type=float; };
#if defined(ENABLE_DOUBLE_SUPPORT)
    template<> struct ScalarType<__m256d>{ using type=double; };
#endif
#endif

#if defined(ENABLE_AVX512_INSTRUCTION_SET) || defined(ENABLE_MIC_INSTRUCTION_SET)
    template<> struct ScalarType<__m512>{ using type=float; };
#if defined(ENABLE_DOUBLE_SUPPORT)
    template<> struct ScalarType<__m512d>{ using type=double; };
#endif
#endif

#if defined(ENABLE_NEON_INSTRUCTION_SET)
    template<> struct ScalarType<float32x4_t>{ using type =float; };
#endif
}



template<class Tw> class Number;
template<class Tw> Number<Tw> min(const Number<Tw>& A, const Number<Tw>& B);
template<class Tw> Number<Tw> max(const Number<Tw>& A, const Number<Tw>& B);
template<class Tw> Number<Tw> blend(const typename NumberPolicy<Number<Tw> >::MASK_TYPE& mask, const Number<Tw>& A, const Number<Tw>& B);
template<class Tw> void Store(typename ScalarType<Tw>::type* data,const Number<Tw>& number);
template<class Tw> void Store(typename ScalarType<Tw>::type& data,const Number<Tw>& number);
template<class Tw> void StoreQuadIn16(typename ScalarType<Tw>::type* data,const Number<Tw>& number,int quad);
template<class Tw> void StoreQuadIn16(typename ScalarType<Tw>::type& data,const Number<Tw>& number,int quad);

//============================================================================//
//============================================================================//

#ifdef ENABLE_IO_SUPPORT
template<class Tw> std::ostream& operator<<( std::ostream& os, const Number<Tw>& number);
#endif

template<class Tw>
struct Number
{
    using T = typename ScalarType<Tw>::type;

    Tw value;
    typedef typename NumberPolicy<Number<Tw> >::MASK_TYPE Mask;
    Number();
    explicit Number(const Mask& mask);

    Number operator+(const Number& other) const;
    Number operator*(const Number& other) const;
    Number operator-(const Number& other) const;
    Number operator/(const Number& other) const;

    Mask operator<(const Number& other) const;
    Mask operator>(const Number& other) const;
    Mask operator<=(const Number& other) const;
    Mask operator>=(const Number& other) const;
    Mask operator==(const Number& other) const;

    Number operator&(const Number& other) const;
    Number operator|(const Number& other) const;
    Number operator^(const Number& other) const;
    Number andnot(const Number& other) const;
    Number operator~() const;

    Number sqrt() const;
    Number rsqrt() const;
    Number log() const;
    Number exp() const;
    Number inverse() const;
    Number abs() const;
    Number sign() const;

    friend Number min<>(const Number& A, const Number& B);
    friend Number max<>(const Number& A, const Number& B);
    friend Number blend<>(const Mask& mask, const Number& A, const Number& B);
    Number mask(const Mask& mask) const;

    void Load_Aligned(const T* data);
    void Load_Aligned(const T& data);

    void Load(const T* data);
    void Load(const T& data);

    void Gather(const T* data,const int* offsets);
    void Gather(const T* data,const int& offsets);

    friend void Store<>(T* data,const Number& number);
    friend void Store<>(T& data,const Number& number);
    friend void StoreQuadIn16<>(T* data,const Number& number,int quad);
    friend void StoreQuadIn16<>(T& data,const Number& number,int quad);

//============================================================================//
//============================================================================//

    Number Spread(int i);
    Number Distribute(int i );
    Number SwizzleAdd(int i, const Number& other);
    Number Horizontal_Quad_Add();
    Number Quad_Mask(int i);

#ifdef ENABLE_IO_SUPPORT
    friend std::ostream& operator<< <>( std:: ostream & os, const Number& number);
#endif

    enum NUMBER_CONSTANTS  {NEG_TEN=0,NEG_NINE,NEG_EIGHT,NEG_SEVEN,NEG_SIX,NEG_FIVE,NEG_FOUR,NEG_THREE,
                            NEG_TWO,NEG_ONE,NEG_ONE_OVER_FOUR,NEG_ONE_HALF,ZERO,ONE_HALF,ONE_OVER_FOUR,
                            ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT,NINE,TEN};



};


//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//

template<> inline
Number<float>::Number()
{value=0.;}

template<> inline
Number<float>::Number(const Mask& mask)
{
    if(mask.value){
        floatConverter ALL_ONES;
        ALL_ONES.i = (0xFFFFFFFF);
        value = ALL_ONES.f;
    }
    else
        value = 0;
}

//==============================================================//
//==============================================================//

//==============================================================//
//                                                              //
//                      BASIC OPERATIONS                        //
//                                                              //
//==============================================================//


//------------------------------------//
//             ADDITION               //
//------------------------------------//

template<> inline
Number<float> Number<float>::operator+(const Number& other) const
{

    Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
    result.value=value+other.value;
#else
    __m128 val1=_mm_load_ss(&value);
    __m128 val2=_mm_load_ss(&other.value);
    __m128 val3=_mm_add_ps(val1,val2);
    _mm_store_ss(&result.value,val3);
#endif
    return result;
}

//------------------------------------//
//           MULTIPLICATION           //
//------------------------------------//

template<> inline
Number<float> Number<float>::operator*(const Number& other) const
{Number result;result.value=value*other.value;return result;}


//------------------------------------//
//           SUBTRACTION              //
//------------------------------------//

template<> inline
Number<float> Number<float>::operator-(const Number& other) const
{Number result;result.value=value-other.value;return result;}

//------------------------------------//
//            DIVISION                //
//------------------------------------//

template<> inline
Number<float> Number<float>::operator/(const Number& other) const
{Number result;result.value=value/other.value;return result;}


//==============================================================//
//                                                              //
//                      COMPARISON OPERATIONS                   //
//                                                              //
//==============================================================//

//------------------------------------//
//             LESS THAN              //
//------------------------------------//

template<> inline
Number<float>::Mask Number<float>::operator<(const Number& other) const
{
    Mask result;
    floatConverter ALL_ONES;
    ALL_ONES.i = (0xFFFFFFFF);
#ifndef FORCE_IDENTICAL_BEHAVIOR
    result.value = value < other.value ? true : false;
#else
    float f_result;
    __m128 val1=_mm_load_ss(&value);
    __m128 val2=_mm_load_ss(&other.value);
    __m128 val3=_mm_cmp_ss(val1,val2,_CMP_LT_OS);
    _mm_store_ss(&f_result,val3);
    ALL_ONES.f=f_result;
    result.value = ALL_ONES.i==0 ? false : true;
#endif
    return result;
}


//------------------------------------//
//             GREATER THAN              //
//------------------------------------//

template<> inline
Number<float>::Mask Number<float>::operator>(const Number& other) const
{
    Mask result;
    floatConverter ALL_ONES;
    ALL_ONES.i = (0xFFFFFFFF);
#ifndef FORCE_IDENTICAL_BEHAVIOR
    result.value = value > other.value ? true : false;
#else
    float f_result;
    __m128 val1=_mm_load_ss(&value);
    __m128 val2=_mm_load_ss(&other.value);
    __m128 val3=_mm_cmp_ss(val1,val2,_CMP_GT_OS);
    _mm_store_ss(&f_result,val3);
    ALL_ONES.f=f_result;
    result.value = ALL_ONES.i==0 ? false : true;
#endif
    return result;
}

//------------------------------------//
//             LESS EQUALS              //
//------------------------------------//

template<> inline
Number<float>::Mask Number<float>::operator<=(const Number& other) const
{
    Mask result;
    floatConverter ALL_ONES;
    ALL_ONES.i = (0xFFFFFFFF);
#ifndef FORCE_IDENTICAL_BEHAVIOR
    result.value = value <= other.value ? true : false;
#else
    float f_result;
    __m128 val1=_mm_load_ss(&value);
    __m128 val2=_mm_load_ss(&other.value);
    __m128 val3=_mm_cmp_ss(val1,val2,_CMP_LE_OS);
    _mm_store_ss(&f_result,val3);
    ALL_ONES.f=f_result;
    result.value = ALL_ONES.i==0 ? false : true;
#endif
    return result;
}


//------------------------------------//
//             GREATER EQUALS              //
//------------------------------------//

template<> inline
Number<float>::Mask Number<float>::operator>=(const Number& other) const
{
    Mask result;
    floatConverter ALL_ONES;
    ALL_ONES.i = (0xFFFFFFFF);
#ifndef FORCE_IDENTICAL_BEHAVIOR
    result.value = value >= other.value ? true : false;
#else
    float f_result;
    __m128 val1=_mm_load_ss(&value);
    __m128 val2=_mm_load_ss(&other.value);
    __m128 val3=_mm_cmp_ss(val1,val2,_CMP_GE_OS);
    _mm_store_ss(&f_result,val3);
    ALL_ONES.f=f_result;
    result.value = ALL_ONES.i==0 ? false : true;
#endif
    return result;
}


//------------------------------------//
//               EQUALS               //
//------------------------------------//

template<> inline
Number<float>::Mask Number<float>::operator==(const Number& other) const
{
    Mask result;
    floatConverter ALL_ONES;
    ALL_ONES.i = (0xFFFFFFFF);
#ifndef FORCE_IDENTICAL_BEHAVIOR
    result.value = value == other.value ? true : false;
#else
    float f_result;
    __m128 val1=_mm_load_ss(&value);
    __m128 val2=_mm_load_ss(&other.value);
    __m128 val3=_mm_cmp_ss(val1,val2,_CMP_EQ_OS);
    _mm_store_ss(&f_result,val3);
    ALL_ONES.f=f_result;
    result.value = ALL_ONES.i==0 ? false : true;
#endif
    return result;
}

//==============================================================//
//                                                              //
//                      BITWISE  OPERATIONS                     //
//                                                              //
//==============================================================//

//------------------------------------//
//             BITWISE AND            //
//------------------------------------//

template<> inline
Number<float> Number<float>::operator&(const Number& other) const
{
    Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
    floatConverter VAL1;
    floatConverter VAL2;
    floatConverter RESULT;
    VAL1.f = value;
    VAL2.f = other.value;
    RESULT.i = VAL1.i & VAL2.i;
    result.value = RESULT.f;
#else
    __m128 val1=_mm_load_ss(&value);
    __m128 val2=_mm_load_ss(&other.value);
    __m128 val3=_mm_and_ps(val1,val2);
    _mm_store_ss(&result.value,val3);
#endif
    return result;
}


//------------------------------------//
//             BITWISE OR             //
//------------------------------------//

template<> inline
Number<float> Number<float>::operator|(const Number& other) const
{
    Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
    floatConverter VAL1;
    floatConverter VAL2;
    floatConverter RESULT;
    VAL1.f = value;
    VAL2.f = other.value;
    RESULT.i = VAL1.i | VAL2.i;
    result.value = RESULT.f;
#else
    __m128 val1=_mm_load_ss(&value);
    __m128 val2=_mm_load_ss(&other.value);
    __m128 val3=_mm_or_ps(val1,val2);
    _mm_store_ss(&result.value,val3);
#endif
    return result;
}


//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//

template<> inline
Number<float> Number<float>::operator^(const Number& other) const
{
    Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
    floatConverter VAL1;
    floatConverter VAL2;
    floatConverter RESULT;
    VAL1.f = value;
    VAL2.f = other.value;
    RESULT.i = VAL1.i ^ VAL2.i;
    result.value = RESULT.f;
#else
    __m128 val1=_mm_load_ss(&value);
    __m128 val2=_mm_load_ss(&other.value);
    __m128 val3=_mm_xor_ps(val1,val2);
    _mm_store_ss(&result.value,val3);
#endif
    return result;
}


//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//

template<> inline
Number<float> Number<float>::operator~() const
{
    Number result;
    floatConverter ALL_ONES;
    ALL_ONES.i = (0xFFFFFFFF);
#ifndef FORCE_IDENTICAL_BEHAVIOR
    floatConverter VAL1;
    floatConverter RESULT;
    VAL1.f = value;
    RESULT.i = ~VAL1.i & ALL_ONES.i;
    result.value = RESULT.f;
#else
    const float __one = ALL_ONES.f;
    __m128 val1=_mm_load_ss(&value);
    __m128 val2=_mm_load_ss(&__one);
    __m128 val3=_mm_andnot_ps(val1,val2);
    _mm_store_ss(&result.value,val3);
#endif
    return result;
}

//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//

template<> inline
Number<float> Number<float>::andnot(const Number& other) const
{
    Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
    floatConverter VAL1;
    floatConverter VAL2;
    floatConverter RESULT;
    VAL1.f = value;
    VAL2.f = other.value;
    RESULT.i = ~VAL1.i & VAL2.i;
    result.value = RESULT.f;
#else
    const float __one = 1.0f;
    __m128 val1=_mm_load_ss(&value);
    __m128 val2=_mm_load_ss(&other.value);
    __m128 val3=_mm_andnot_ps(val1,val2);
    _mm_store_ss(&result.value,val3);
#endif
    return result;
}


//==============================================================//
//==============================================================//

//==============================================================//
//                                                              //
//                  ADVANCED OPERATIONS                         //
//                                                              //
//==============================================================//

//------------------------------------//
//              SQUARE ROOT           //
//------------------------------------//

template<> inline
Number<float> Number<float>::sqrt() const
{Number result;result.value=::sqrt(value);return result;}


//------------------------------------//
//       RECIPROCAL SQUARE ROOT       //
//------------------------------------//

template<> inline
Number<float> Number<float>::rsqrt() const
{
Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
 result.value=1.0/::sqrt(value);
#else
    __m128 val1=_mm_broadcast_ss(&value);
    val1=_mm_rsqrt_ss(val1);
     _mm_store_ss(&result.value,val1);
#endif
 return result;
}

//------------------------------------//
//                 LOG                //
//------------------------------------//

template<> inline
Number<float> Number<float>::log() const
{
    Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
    result.value=::log(value);
#else
    __m128 val1=_mm_broadcast_ss(&value);
    val1=_mm_log_ps(val1);
     _mm_store_ss(&result.value,val1);
#endif
    return result;
}



//------------------------------------//
//                 EXP                //
//------------------------------------//

template<> inline
Number<float> Number<float>::exp() const
{
    Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
    result.value=::exp(value);
#else
    __m128 val1=_mm_broadcast_ss(&value);
    val1=_mm_exp_ps(val1);
     _mm_store_ss(&result.value,val1);
#endif
    return result;
}


//------------------------------------//
//        INVERSE (RECIPROCAL)        //
//------------------------------------//

template<> inline
Number<float> Number<float>::inverse() const
{
    Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
    result.value=1.f/value;
#else
    // 22-bit precision -- See Intel IA32 Architecture Optimization Manual
    __m128 val1=_mm_broadcast_ss(&value);
    __m128 val3=_mm_rcp_ps(val1);
    __m128 val2=_mm_add_ss(val3,val3);
    val3=_mm_mul_ps(val3,val3);
    val3=_mm_mul_ps(val3,val1);
    val2=_mm_sub_ps(val2,val3);
    _mm_store_ss(&result.value,val2);
#endif
    return result;
}



//------------------------------------//
//        Absolulte Value             //
//------------------------------------//

template<> inline
Number<float> Number<float>::abs() const
{
    Number result;
    result.value=(value>=0.f)?value:-value;
    return result;
}

//------------------------------------//
//               SIGN                 //
//------------------------------------//

template<> inline
Number<float> Number<float>::sign() const
{
    Number result;
    result.value=(value<0.f)?-1.0f:1.0f;
    return result;
}

//------------------------------------//
//               MINIMUM              //
//------------------------------------//

template<> inline
Number<float> min(const Number<float>& A, const Number<float>& B)
{Number<float> result;result.value=A.value < B.value ? A.value : B.value;return result;}

//------------------------------------//
//               MAXIMUM              //
//------------------------------------//

template<> inline
Number<float> max(const Number<float>& A, const Number<float>& B)
{Number<float> result;result.value=A.value > B.value ? A.value : B.value;return result;}


//------------------------------------//
//      Masked assignment (blend)     //
//------------------------------------//

template<> inline
Number<float> blend(const NumberPolicy<Number<float> >::MASK_TYPE& mask, const Number<float>& A, const Number<float>& B)
{
#ifndef FORCE_IDENTICAL_BEHAVIOR
    Number<float> result;
    result.value = mask.value ? B.value : A.value;
    return result;
#else
    Number<float> result;
    floatConverter ALL_ONES,ALL_ZEROS;
    ALL_ONES.i = (0xFFFFFFFF);
    ALL_ZEROS.i = (0x0);
    __m128 val1=_mm_load_ss(&A.value);
    __m128 val2=_mm_load_ss(&B.value);
    __m128 val3;
    if(mask.value)
        val3 =_mm_load_ss(&ALL_ONES.f);
    else
        val3=_mm_load_ss(&ALL_ZEROS.f);
    val1=_mm_blendv_ps(val1,val2,val3);
    _mm_store_ss(&result.value,val1);
    return result;
#endif
}

//------------------------------------//
//      Masked assignment (mask )     //
//------------------------------------//
template<> inline
Number<float> Number<float>::mask(const Mask& mask) const
{
#ifndef FORCE_IDENTICAL_BEHAVIOR
    Number<float> result;
    result.value = mask.value ? value : (float)(0.0);
    return result;
#else
    Number<float> result;
    floatConverter ALL_ONES,ALL_ZEROS;
    ALL_ONES.i = (0xFFFFFFFF);
    ALL_ZEROS.i = (0x0);
    __m128 val1=_mm_load_ss(&value);
    __m128 val2;
    if(mask.value)
        val2 =_mm_load_ss(&ALL_ONES.f);
    else
        val2=_mm_load_ss(&ALL_ZEROS.f);
    val1=_mm_and_ps(val1,val2);
    _mm_store_ss(&result.value,val1);
    return result;
#endif
}



//==============================================================//
//==============================================================//

#ifdef ENABLE_IO_SUPPORT
//==============================================================//
//                                                              //
//                  DEBUG PRINT SUPPORT                         //
//                                                              //
//==============================================================//

template<> inline
std:: ostream & operator<<( std:: ostream & os, const Number<float>& number)
{
    os << "[ " << number.value << " ]";
    return os;
}

//==============================================================//
//==============================================================//
#endif


//==============================================================//
//                                                              //
//                  LOADS AND STORES                            //
//                                                              //
//==============================================================//

//------------------------------------//
//              LOADS                 //
//------------------------------------//

template<> inline
void Number<float>::Load(const float* data)
{value=*data;}


template<> inline
void Number<float>::Load(const float& data)
{value=data;}


//------------------------------------//
//              LOADS                 //
//------------------------------------//

template<> inline
void Number<float>::Gather(const float* data,const int* offsets)
{value=data[*offsets];}

template<> inline
void Number<float>::Gather(const float* data,const int& offsets)
{value=data[offsets];}


//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//

template<> inline
void Number<float>::Load_Aligned(const float* data)
{value=*data;}

template<> inline
void Number<float>::Load_Aligned(const float& data)
{value=data;}

//------------------------------------//
//             STORES                 //
//------------------------------------//

template<> inline
void Store(float* data,const Number<float>& number)
{*data=number.value;}

template<> inline
void Store(float& data,const Number<float>& number)
{data=number.value;}
//==============================================================//
//                                                              //
//                  HORIZONTAL OPERATIONS                       //
//                                                              //
//==============================================================//

template<> inline
Number<float> Number<float>::Spread(int i){
    return *this;
};


template<> inline
Number<float> Number<float>::Distribute(int i){
    return *this;
};


template<> inline
Number<float> Number<float>::SwizzleAdd(int i, const Number& other){
    return *this;
};

template<> inline
Number<float> Number<float>::Horizontal_Quad_Add(){
    return *this;
};

template<> inline
Number<float> Number<float>::Quad_Mask(int i){
    return *this;
};

template<> inline
void StoreQuadIn16(float& data, const Number<float>& number, int quad){
    data=number.value;
}

template<> inline
void StoreQuadIn16(float* data, const Number<float>& number, int quad){
    *data=number.value;
}

//==============================================================//
//==============================================================//
#if 0

#endif

#endif
