//#####################################################################
//  Copyright (c) 2018 Haixiang Liu.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_NUMBER_AVX_DOUBLE_H__
#define __KERNEL_NUMBER_AVX_DOUBLE_H__


//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//


template<> inline
Number<__m256d>::Number()
{value=_mm256_xor_pd(value,value);}



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
Number<__m256d> Number<__m256d>::operator+(const Number& other) const
{Number<__m256d> result;result.value=_mm256_add_pd(value,other.value);return result;}



//------------------------------------//
//           MULTIPLICATION           //
//------------------------------------//


template<> inline
Number<__m256d> Number<__m256d>::operator*(const Number& other) const
{Number<__m256d> result;result.value=_mm256_mul_pd(value,other.value);return result;}


//------------------------------------//
//           SUBTRACTION              //
//------------------------------------//


template<> inline
Number<__m256d> Number<__m256d>::operator-(const Number& other) const
{Number<__m256d> result;result.value=_mm256_sub_pd(value,other.value);return result;}


//------------------------------------//
//            DIVISION                //
//------------------------------------//


template<> inline
Number<__m256d> Number<__m256d>::operator/(const Number& other) const
{Number<__m256d> result;result.value=_mm256_div_pd(value,other.value);return result;}


//==============================================================//
//                                                              //
//                      COMPARISON OPERATIONS                   //
//                                                              //
//==============================================================//

//------------------------------------//
//             LESS THAN              //
//------------------------------------//


template<> inline
Number<__m256d>::Mask Number<__m256d>::operator<(const Number& other) const
{Mask result;result.value=_mm256_cmp_pd(value,other.value,_CMP_LT_OS);return result;}


//------------------------------------//
//             GREATER THAN              //
//------------------------------------//


template<> inline
Number<__m256d>::Mask Number<__m256d>::operator>(const Number& other) const
{Mask result;result.value=_mm256_cmp_pd(value,other.value,_CMP_GT_OS);return result;}


//------------------------------------//
//             LESS EQUALS              //
//------------------------------------//


template<> inline
Number<__m256d>::Mask Number<__m256d>::operator<=(const Number& other) const
{Mask result;result.value=_mm256_cmp_pd(value,other.value,_CMP_LE_OS);return result;}


//------------------------------------//
//             GREATER EQUALS              //
//------------------------------------//



template<> inline
Number<__m256d>::Mask Number<__m256d>::operator>=(const Number& other) const
{Mask result;result.value=_mm256_cmp_pd(value,other.value,_CMP_GE_OS);return result;}



//------------------------------------//
//               EQUALS               //
//------------------------------------//


template<> inline
Number<__m256d>::Mask Number<__m256d>::operator==(const Number& other) const
{Mask result;result.value=_mm256_cmp_pd(value,other.value,_CMP_EQ_OS);return result;}


//==============================================================//
//                                                              //
//                      BITWISE  OPERATIONS                     //
//                                                              //
//==============================================================//

//------------------------------------//
//             BITWISE AND            //
//------------------------------------//



template<> inline
Number<__m256d> Number<__m256d>::operator&(const Number& other) const
{Number<__m256d> result;result.value=_mm256_and_pd(value,other.value);return result;}


//------------------------------------//
//             BITWISE OR             //
//------------------------------------//


template<> inline
Number<__m256d> Number<__m256d>::operator|(const Number& other) const
{Number<__m256d> result;result.value=_mm256_or_pd(value,other.value);return result;}


//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//



template<> inline
Number<__m256d> Number<__m256d>::operator^(const Number& other) const
{Number<__m256d> result;result.value=_mm256_xor_pd(value,other.value);return result;}


//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//


template<> inline
Number<__m256d> Number<__m256d>::operator~() const
{
    Number<__m256d> result;
    doubleConverter X;
    X.i = 0xFFFFFFFFFFFF;
    BUILD_CONSTANT_4_DOUBLE(__one, X.d);
    __m256d val2=_mm256_load_pd(__one);
    result.value=_mm256_andnot_pd(value,val2);
     return result;
}


//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//



template<> inline
Number<__m256d> Number<__m256d>::andnot(const Number& other) const
{
    Number<__m256d> result;
    result.value=_mm256_andnot_pd(value,other.value);
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
Number<__m256d> Number<__m256d>::sqrt() const
{Number<__m256d> result;result.value=_mm256_sqrt_pd(value);return result;}


//------------------------------------//
//       RECIPROCAL SQUARE ROOT       //
//------------------------------------//


template<> inline
Number<__m256d> Number<__m256d>::rsqrt() const
{Number<__m256d> result;result.value=_mm256_div_pd(_mm256_set1_pd(1.0),_mm256_sqrt_pd(value));;return result;}



//------------------------------------//
//                 LOG                //
//------------------------------------//



template<> inline
Number<__m256d> Number<__m256d>::log() const
{
    Number<__m256d> result;
#ifdef __INTEL_COMPILER
#ifndef FORCE_IDENTICAL_BEHAVIOR
    result.value=_mm256_log_pd(value);
#else
    __m128d val1=_mm256_extractf128_pd(value,0);
    __m128d val2=_mm256_extractf128_pd(value,1);
    val1 = _mm_log_pd(val1);
    val2 = _mm_log_pd(val2);
    result.value = _mm256_insertf128_pd(result.value,val1,0);
    result.value = _mm256_insertf128_pd(result.value,val2,1);
#endif
#else

#endif
return result;
}



//------------------------------------//
//                 EXP                //
//------------------------------------//



template<> inline
Number<__m256d> Number<__m256d>::exp() const
{
    Number<__m256d> result;
#ifdef __INTEL_COMPILER
#ifndef FORCE_IDENTICAL_BEHAVIOR
    result.value=_mm256_exp_pd(value);
#else
    __m128d val1=_mm256_extractf128_pd(value,0);
    __m128d val2=_mm256_extractf128_pd(value,1);
    val1 = _mm_exp_pd(val1);
    val2 = _mm_exp_pd(val2);
    result.value = _mm256_insertf128_pd(result.value,val1,0);
    result.value = _mm256_insertf128_pd(result.value,val2,1);
#endif
#else

#endif
return result;
}



//------------------------------------//
//        INVERSE (RECIPROCAL)        //
//------------------------------------//


template<> inline
Number<__m256d> Number<__m256d>::inverse() const
{
    // 22-bit precision -- See Intel IA32 Architecture Optimization Manual
   Number<__m256d> result;   
   __m256d val3=_mm256_div_pd(_mm256_set1_pd(1.0),value);
   __m256d val2=_mm256_add_pd(val3,val3);
   val3=_mm256_mul_pd(val3,val3);
   val3=_mm256_mul_pd(val3,value);  
   result.value=_mm256_sub_pd(val2,val3);
    
    return result;
}



//------------------------------------//
//        Absolulte Value             //
//------------------------------------//


template<> inline
Number<__m256d> Number<__m256d>::abs() const
{
    Number<__m256d> neg,result;
    neg.value=_mm256_sub_pd(neg.value,value);
    result.value=_mm256_max_pd(value,neg.value);
    return result;
}



//------------------------------------//
//               SIGN                 //
//------------------------------------//


template<> inline
Number<__m256d> Number<__m256d>::sign() const
{
    Number<__m256d> zero,negs,poss,result;
    BUILD_CONSTANT_4_DOUBLE(__one, 1.0);
    __m256d none=_mm256_load_pd(__one);
    none = _mm256_sub_pd(zero.value,none);
    __m256d one=_mm256_load_pd(__one);
    negs.value = _mm256_cmp_pd(value,zero.value,_CMP_GE_OS);
    result.value = _mm256_blendv_pd(none,one,negs.value);
    return result;
}



//------------------------------------//
//               MINIMUM              //
//------------------------------------//


template<> inline
Number<__m256d> min(const Number<__m256d>& A, const Number<__m256d>& B)
{Number<__m256d> result;result.value=_mm256_min_pd(A.value, B.value);return result;}


//------------------------------------//
//               MAXIMUM              //
//------------------------------------//


template<> inline
Number<__m256d> max(const Number<__m256d>& A, const Number<__m256d>& B)
{Number<__m256d> result;result.value=_mm256_max_pd(A.value, B.value);return result;}



//------------------------------------//
//      Masked assignment (blend)     //
//------------------------------------//


template<> inline
Number<__m256d> blend(const typename Number<__m256d>::Mask& mask, const Number<__m256d>& A, const Number<__m256d>& B)
{Number<__m256d> result;result.value=_mm256_blendv_pd(A.value,B.value,mask.value);return result;}

//------------------------------------//
//      Masked assignment (mask)     //
//------------------------------------//


template<> inline
Number<__m256d> Number<__m256d>::mask(const Mask& mask) const 
{Number<__m256d> result, temp;result.value=_mm256_blendv_pd(temp.value,value,mask.value);return result;}

//==============================================================//
//==============================================================//

#ifdef ENABLE_IO_SUPPORT
//==============================================================//
//                                                              //
//                  DEBUG PRINT SUPPORT                         //
//                                                              //
//==============================================================//


template<> inline
std:: ostream & operator<<( std:: ostream & os, const Number<__m256d>& number)
{
    double *fp = (double*)&number.value;
    os << "[ " << *(fp+3)
       << " " << *(fp+2)
       << " " << *(fp+1)
       << " " << *(fp) << " ]";
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
void Number<__m256d>::Load(const double* data)
{value=_mm256_loadu_pd(data);}




template<> inline
void Number<__m256d>::Load(const double& data)
{value=_mm256_loadu_pd(&data);}

//------------------------------------//
//              LOADS                 //
//------------------------------------//

#define UNTESTED

#if defined(UNTESTED)
#else
template<> inline
void Number<__m256d>::Gather(const double* data,const int* offsets)
{
    __m128d val1=_mm_load_ss(data+offsets[0]);
    __m128d val2=_mm_load_ss(data+offsets[1]);
    __m128d val3=_mm_load_ss(data+offsets[2]);
    __m128d val4=_mm_load_ss(data+offsets[3]);
    val1=_mm_insert_pd(val1,val2,0x10);
    val3=_mm_insert_pd(val3,val4,0x10);
    val1=_mm_movelh_pd(val1,val3);
    value=_mm256_insertf128_pd(value,val1,0);
    val1=_mm_load_ss(data+offsets[4]);
    val2=_mm_load_ss(data+offsets[5]);
    val3=_mm_load_ss(data+offsets[6]);
    val4=_mm_load_ss(data+offsets[7]);
    val1=_mm_insert_pd(val1,val2,0x10);
    val3=_mm_insert_pd(val3,val4,0x10);
    val1=_mm_movelh_pd(val1,val3);
    value=_mm256_insertf128_pd(value,val1,1);
}
#endif


#if defined(UNTESTED)
#else

template<> inline
void Number<__m256d>::Gather(const double* data,const int& offsets)
{
    __m128d val1=_mm_load_ss(data+offsets);
    __m128d val2=_mm_load_ss(data+(&offsets)[1]);
    __m128d val3=_mm_load_ss(data+(&offsets)[2]);
    __m128d val4=_mm_load_ss(data+(&offsets)[3]);
    val1=_mm_insert_pd(val1,val2,0x10);
    val3=_mm_insert_pd(val3,val4,0x10);
    val1=_mm_movelh_pd(val1,val3);
    value=_mm256_insertf128_pd(value,val1,0);
    val1=_mm_load_ss(data+(&offsets)[4]);
    val2=_mm_load_ss(data+(&offsets)[5]);
    val3=_mm_load_ss(data+(&offsets)[6]);
    val4=_mm_load_ss(data+(&offsets)[7]);
    val1=_mm_insert_pd(val1,val2,0x10);
    val3=_mm_insert_pd(val3,val4,0x10);
    val1=_mm_movelh_pd(val1,val3);
    value=_mm256_insertf128_pd(value,val1,1);
}
#endif


//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//


template<> inline
void Number<__m256d>::Load_Aligned(const double* data)
{value=_mm256_load_pd(data);}





template<> inline
void Number<__m256d>::Load_Aligned(const double& data)
{value=_mm256_load_pd(&data);}



//------------------------------------//
//             STORES                 //
//------------------------------------//


template<> inline
void Store(double* data,const Number<__m256d>& number)
{_mm256_store_pd(data,number.value);}



template<> inline
void Store(double& data,const Number<__m256d>& number)
{_mm256_store_pd(&data,number.value);}


//==============================================================//
//                                                              //
//                  HORIZONTAL OPERATIONS                       //
//                                                              //
//==============================================================//

template<> inline
Number<__m256d> Number<__m256d>::Spread(int i){
    Number<__m256d> result;
    switch(i){
    case 0:
        result.value = _mm256_permute2f128_pd( value, value, 0x00 );
        break;
    case 1:
        result.value = _mm256_permute2f128_pd( value, value, 0x11 );
        break;
    }
    return result;
};

#if defined(UNTESTED)
#else
template<> inline
Number<__m256d> Number<__m256d>::Distribute(int i){
    // 0 : AAAA BBBB
    // 2 : CCCC DDDD

    Number<__m256d> result, tmp;
    // __mmask16 mask1, mask2;

    switch(i){
    case 0:
        tmp.value = _mm256_permute_pd( value, 0x00 );
        result.value = _mm256_permute_pd( value, 0x55 );
        result.value = _mm256_blend_pd(tmp.value, result.value, 0xF0 );
        break;
    case 2:
        tmp.value = _mm256_permute_pd( value, 0xAA );
        result.value = _mm256_permute_pd( value, 0xFF );
        result.value = _mm256_blend_pd(tmp.value, result.value, 0xF0 );
        break;
    }
    return  result;
};
#endif

#if defined(UNTESTED)
#else
template<> inline
Number<__m256d> Number<__m256d>::SwizzleAdd(int i, const Number& other){
    // 2: NOP
    // 1: SWZ/ADD 
    Number<__m256d> result;
    result.value = value;
    if( i == 1 ){
        result.value = _mm256_permute2f128_pd( value, value, 0x01 );
        result.value = _mm256_add_pd( other.value, result.value );
    }
    return result;
};
#endif

template<> inline
Number<__m256d> Number<__m256d>::Horizontal_Quad_Add(){
    Number<__m256d> result;
    result.value = _mm256_hadd_pd( value, value );
    result.value = _mm256_hadd_pd( result.value, result.value );
    return result;
};

#if defined(UNTESTED)
#else
template<> inline
Number<__m256d> Number<__m256d>::Quad_Mask(int i){
    Number<__m256d> result;
    Number<__m256d> tmp;
    Number<__m256d> zero;
    switch( i ) {
        case 0: 
            result.value = _mm256_blend_pd( zero.value, value, 0x21); // [aaaa bbbb] --> [a000 0b00]
            //result.value = _mm256_permute2f128_pd(result.value, result.value, 0x00 ); // --> [0b00 a000]
            //result.value = _mm256_blend_pd( result.value, tmp.value, 0x12);
            break;
        case 1: 
            result.value = _mm256_blend_pd( zero.value, value, 0x84);
            // result.value = _mm256_permute2f128_pd(result.value, result.value, 0x00 ); // --> [0b00 a000]
            break;
    }
    return result;
};
#endif

#if defined(UNTESTED)
#else
template<> inline
void StoreQuadIn16(double& data, const Number<__m256d>& number, int quad){
    __m128d tmp;
    switch(quad){
    case 0:
    case 2:
        tmp = _mm256_extractf128_pd( number.value, 0x00 );
        break;
    case 1:
    case 3:
        tmp =  _mm256_extractf128_pd( number.value, 0x01 );
        break;
    }
    _mm_store_pd(((&data) + quad*4),tmp);
}
#endif

#if defined(UNTESTED)
#else
template<> inline
void StoreQuadIn16(double* data, const Number<__m256d>& number, int quad){
    __m128d tmp;
    switch(quad){
    case 0:
    case 2:
        tmp = _mm256_extractf128_pd( number.value, 0x00 );
        break;
    case 1:
    case 3:
        tmp =  _mm256_extractf128_pd( number.value, 0x01 );
        break;
    }
    _mm_store_pd(((data) + quad*4),tmp);
}
#endif
//==============================================================//
//==============================================================//



#endif
