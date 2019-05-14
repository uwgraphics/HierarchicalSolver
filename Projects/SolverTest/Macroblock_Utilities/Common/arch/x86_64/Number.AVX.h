//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_NUMBER_AVX_H__
#define __KERNEL_NUMBER_AVX_H__


//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//


template<> inline
Number<__m256>::Number()
{value=_mm256_xor_ps(value,value);}



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
Number<__m256> Number<__m256>::operator+(const Number& other) const
{Number<__m256> result;result.value=_mm256_add_ps(value,other.value);return result;}



//------------------------------------//
//           MULTIPLICATION           //
//------------------------------------//


template<> inline
Number<__m256> Number<__m256>::operator*(const Number& other) const
{Number<__m256> result;result.value=_mm256_mul_ps(value,other.value);return result;}


//------------------------------------//
//           SUBTRACTION              //
//------------------------------------//


template<> inline
Number<__m256> Number<__m256>::operator-(const Number& other) const
{Number<__m256> result;result.value=_mm256_sub_ps(value,other.value);return result;}


//------------------------------------//
//            DIVISION                //
//------------------------------------//


template<> inline
Number<__m256> Number<__m256>::operator/(const Number& other) const
{Number<__m256> result;result.value=_mm256_div_ps(value,other.value);return result;}


//==============================================================//
//                                                              //
//                      COMPARISON OPERATIONS                   //
//                                                              //
//==============================================================//

//------------------------------------//
//             LESS THAN              //
//------------------------------------//


template<> inline
Number<__m256>::Mask Number<__m256>::operator<(const Number& other) const
{Mask result;result.value=_mm256_cmp_ps(value,other.value,_CMP_LT_OS);return result;}


//------------------------------------//
//             GREATER THAN              //
//------------------------------------//


template<> inline
Number<__m256>::Mask Number<__m256>::operator>(const Number& other) const
{Mask result;result.value=_mm256_cmp_ps(value,other.value,_CMP_GT_OS);return result;}


//------------------------------------//
//             LESS EQUALS              //
//------------------------------------//


template<> inline
Number<__m256>::Mask Number<__m256>::operator<=(const Number& other) const
{Mask result;result.value=_mm256_cmp_ps(value,other.value,_CMP_LE_OS);return result;}


//------------------------------------//
//             GREATER EQUALS              //
//------------------------------------//



template<> inline
Number<__m256>::Mask Number<__m256>::operator>=(const Number& other) const
{Mask result;result.value=_mm256_cmp_ps(value,other.value,_CMP_GE_OS);return result;}



//------------------------------------//
//               EQUALS               //
//------------------------------------//


template<> inline
Number<__m256>::Mask Number<__m256>::operator==(const Number& other) const
{Mask result;result.value=_mm256_cmp_ps(value,other.value,_CMP_EQ_OS);return result;}


//==============================================================//
//                                                              //
//                      BITWISE  OPERATIONS                     //
//                                                              //
//==============================================================//

//------------------------------------//
//             BITWISE AND            //
//------------------------------------//



template<> inline
Number<__m256> Number<__m256>::operator&(const Number& other) const
{Number<__m256> result;result.value=_mm256_and_ps(value,other.value);return result;}


//------------------------------------//
//             BITWISE OR             //
//------------------------------------//


template<> inline
Number<__m256> Number<__m256>::operator|(const Number& other) const
{Number<__m256> result;result.value=_mm256_or_ps(value,other.value);return result;}


//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//



template<> inline
Number<__m256> Number<__m256>::operator^(const Number& other) const
{Number<__m256> result;result.value=_mm256_xor_ps(value,other.value);return result;}


//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//


template<> inline
Number<__m256> Number<__m256>::operator~() const
{
    Number<__m256> result;
    floatConverter X;
    X.i = 0xFFFFFFFF;
    BUILD_CONSTANT_8(__one, X.f);
    __m256 val2=_mm256_load_ps(__one);
    result.value=_mm256_andnot_ps(value,val2);
     return result;
}


//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//



template<> inline
Number<__m256> Number<__m256>::andnot(const Number& other) const
{
    Number<__m256> result;
    result.value=_mm256_andnot_ps(value,other.value);
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
Number<__m256> Number<__m256>::sqrt() const
{Number<__m256> result;result.value=_mm256_sqrt_ps(value);return result;}


//------------------------------------//
//       RECIPROCAL SQUARE ROOT       //
//------------------------------------//


template<> inline
Number<__m256> Number<__m256>::rsqrt() const
{Number<__m256> result;result.value=_mm256_rsqrt_ps(value);return result;}



//------------------------------------//
//                 LOG                //
//------------------------------------//



template<> inline
Number<__m256> Number<__m256>::log() const
{
    Number<__m256> result;
#ifdef __INTEL_COMPILER
#ifndef FORCE_IDENTICAL_BEHAVIOR
    result.value=_mm256_log_ps(value);
#else
    __m128 val1=_mm256_extractf128_ps(value,0);
    __m128 val2=_mm256_extractf128_ps(value,1);
    val1 = _mm_log_ps(val1);
    val2 = _mm_log_ps(val2);
    result.value = _mm256_insertf128_ps(result.value,val1,0);
    result.value = _mm256_insertf128_ps(result.value,val2,1);
#endif
#else

#endif
return result;
}



//------------------------------------//
//                 EXP                //
//------------------------------------//



template<> inline
Number<__m256> Number<__m256>::exp() const
{
    Number<__m256> result;
#ifdef __INTEL_COMPILER
#ifndef FORCE_IDENTICAL_BEHAVIOR
    result.value=_mm256_exp_ps(value);
#else
    __m128 val1=_mm256_extractf128_ps(value,0);
    __m128 val2=_mm256_extractf128_ps(value,1);
    val1 = _mm_exp_ps(val1);
    val2 = _mm_exp_ps(val2);
    result.value = _mm256_insertf128_ps(result.value,val1,0);
    result.value = _mm256_insertf128_ps(result.value,val2,1);
#endif
#else

#endif
return result;
}



//------------------------------------//
//        INVERSE (RECIPROCAL)        //
//------------------------------------//


template<> inline
Number<__m256> Number<__m256>::inverse() const
{
    // 22-bit precision -- See Intel IA32 Architecture Optimization Manual
   Number<__m256> result;
   __m256 val3=_mm256_rcp_ps(value);
   __m256 val2=_mm256_add_ps(val3,val3);
   val3=_mm256_mul_ps(val3,val3);
   val3=_mm256_mul_ps(val3,value);  
   result.value=_mm256_sub_ps(val2,val3);
    
    return result;
}



//------------------------------------//
//        Absolulte Value             //
//------------------------------------//


template<> inline
Number<__m256> Number<__m256>::abs() const
{
    Number<__m256> neg,result;
    neg.value=_mm256_sub_ps(neg.value,value);
    result.value=_mm256_max_ps(value,neg.value);
    return result;
}



//------------------------------------//
//               SIGN                 //
//------------------------------------//


template<> inline
Number<__m256> Number<__m256>::sign() const
{
    Number<__m256> zero,negs,poss,result;
    BUILD_CONSTANT_8(__one, 1.0f);
    __m256 none=_mm256_load_ps(__one);
    none = _mm256_sub_ps(zero.value,none);
    __m256 one=_mm256_load_ps(__one);
    negs.value = _mm256_cmp_ps(value,zero.value,_CMP_GE_OS);
    result.value = _mm256_blendv_ps(none,one,negs.value);
    return result;
}



//------------------------------------//
//               MINIMUM              //
//------------------------------------//


template<> inline
Number<__m256> min(const Number<__m256>& A, const Number<__m256>& B)
{Number<__m256> result;result.value=_mm256_min_ps(A.value, B.value);return result;}


//------------------------------------//
//               MAXIMUM              //
//------------------------------------//


template<> inline
Number<__m256> max(const Number<__m256>& A, const Number<__m256>& B)
{Number<__m256> result;result.value=_mm256_max_ps(A.value, B.value);return result;}



//------------------------------------//
//      Masked assignment (blend)     //
//------------------------------------//


template<> inline
Number<__m256> blend(const typename Number<__m256>::Mask& mask, const Number<__m256>& A, const Number<__m256>& B)
{Number<__m256> result;result.value=_mm256_blendv_ps(A.value,B.value,mask.value);return result;}

//------------------------------------//
//      Masked assignment (mask)     //
//------------------------------------//


template<> inline
Number<__m256> Number<__m256>::mask(const Mask& mask) const 
{Number<__m256> result, temp;result.value=_mm256_blendv_ps(temp.value,value,mask.value);return result;}

//==============================================================//
//==============================================================//

#ifdef ENABLE_IO_SUPPORT
//==============================================================//
//                                                              //
//                  DEBUG PRINT SUPPORT                         //
//                                                              //
//==============================================================//


template<> inline
std:: ostream & operator<<( std:: ostream & os, const Number<__m256>& number)
{
    float *fp = (float*)&number.value;
    os << "[ " << *(fp+7)
       << " " << *(fp+6)
       << " " << *(fp+5)
       << " " << *(fp+4)
       << " " << *(fp+3)
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
void Number<__m256>::Load(const float* data)
{value=_mm256_loadu_ps(data);}




template<> inline
void Number<__m256>::Load(const float& data)
{value=_mm256_loadu_ps(&data);}

//------------------------------------//
//              LOADS                 //
//------------------------------------//



template<> inline
void Number<__m256>::Gather(const float* data,const int* offsets)
{
    __m128 val1=_mm_load_ss(data+offsets[0]);
    __m128 val2=_mm_load_ss(data+offsets[1]);
    __m128 val3=_mm_load_ss(data+offsets[2]);
    __m128 val4=_mm_load_ss(data+offsets[3]);
    val1=_mm_insert_ps(val1,val2,0x10);
    val3=_mm_insert_ps(val3,val4,0x10);
    val1=_mm_movelh_ps(val1,val3);
    value=_mm256_insertf128_ps(value,val1,0);
    val1=_mm_load_ss(data+offsets[4]);
    val2=_mm_load_ss(data+offsets[5]);
    val3=_mm_load_ss(data+offsets[6]);
    val4=_mm_load_ss(data+offsets[7]);
    val1=_mm_insert_ps(val1,val2,0x10);
    val3=_mm_insert_ps(val3,val4,0x10);
    val1=_mm_movelh_ps(val1,val3);
    value=_mm256_insertf128_ps(value,val1,1);
}




template<> inline
void Number<__m256>::Gather(const float* data,const int& offsets)
{
    __m128 val1=_mm_load_ss(data+offsets);
    __m128 val2=_mm_load_ss(data+(&offsets)[1]);
    __m128 val3=_mm_load_ss(data+(&offsets)[2]);
    __m128 val4=_mm_load_ss(data+(&offsets)[3]);
    val1=_mm_insert_ps(val1,val2,0x10);
    val3=_mm_insert_ps(val3,val4,0x10);
    val1=_mm_movelh_ps(val1,val3);
    value=_mm256_insertf128_ps(value,val1,0);
    val1=_mm_load_ss(data+(&offsets)[4]);
    val2=_mm_load_ss(data+(&offsets)[5]);
    val3=_mm_load_ss(data+(&offsets)[6]);
    val4=_mm_load_ss(data+(&offsets)[7]);
    val1=_mm_insert_ps(val1,val2,0x10);
    val3=_mm_insert_ps(val3,val4,0x10);
    val1=_mm_movelh_ps(val1,val3);
    value=_mm256_insertf128_ps(value,val1,1);
}



//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//


template<> inline
void Number<__m256>::Load_Aligned(const float* data)
{value=_mm256_load_ps(data);}





template<> inline
void Number<__m256>::Load_Aligned(const float& data)
{value=_mm256_load_ps(&data);}



//------------------------------------//
//             STORES                 //
//------------------------------------//


template<> inline
void Store(float* data,const Number<__m256>& number)
{_mm256_store_ps(data,number.value);}



template<> inline
void Store(float& data,const Number<__m256>& number)
{_mm256_store_ps(&data,number.value);}


//==============================================================//
//                                                              //
//                  HORIZONTAL OPERATIONS                       //
//                                                              //
//==============================================================//

template<> inline
Number<__m256> Number<__m256>::Spread(int i){
    Number<__m256> result;
    switch(i){
    case 0:
        result.value = _mm256_permute2f128_ps( value, value, 0x00 );
        break;
    case 1:
        result.value = _mm256_permute2f128_ps( value, value, 0x11 );
        break;
    }
    return result;
};


template<> inline
Number<__m256> Number<__m256>::Distribute(int i){
    // 0 : AAAA BBBB
    // 2 : CCCC DDDD

    Number<__m256> result, tmp;
    // __mmask16 mask1, mask2;

    switch(i){
    case 0:
        tmp.value = _mm256_permute_ps( value, 0x00 );
        result.value = _mm256_permute_ps( value, 0x55 );
        result.value = _mm256_blend_ps(tmp.value, result.value, 0xF0 );
        break;
    case 2:
        tmp.value = _mm256_permute_ps( value, 0xAA );
        result.value = _mm256_permute_ps( value, 0xFF );
        result.value = _mm256_blend_ps(tmp.value, result.value, 0xF0 );
        break;
    }
    return  result;
};


template<> inline
Number<__m256> Number<__m256>::SwizzleAdd(int i, const Number& other){
    // 2: NOP
    // 1: SWZ/ADD 
    Number<__m256> result;
    result.value = value;
    if( i == 1 ){
        result.value = _mm256_permute2f128_ps( value, value, 0x01 );
        result.value = _mm256_add_ps( other.value, result.value );
    }
    return result;
};

template<> inline
Number<__m256> Number<__m256>::Horizontal_Quad_Add(){
    Number<__m256> result;
    result.value = _mm256_hadd_ps( value, value );
    result.value = _mm256_hadd_ps( result.value, result.value );
    return result;
};

template<> inline
Number<__m256> Number<__m256>::Quad_Mask(int i){
    Number<__m256> result;
    Number<__m256> tmp;
    Number<__m256> zero;
    switch( i ) {
        case 0: 
            result.value = _mm256_blend_ps( zero.value, value, 0x21); // [aaaa bbbb] --> [a000 0b00]
            //result.value = _mm256_permute2f128_ps(result.value, result.value, 0x00 ); // --> [0b00 a000]
            //result.value = _mm256_blend_ps( result.value, tmp.value, 0x12);
            break;
        case 1: 
            result.value = _mm256_blend_ps( zero.value, value, 0x84);
            // result.value = _mm256_permute2f128_ps(result.value, result.value, 0x00 ); // --> [0b00 a000]
            break;
    }
    return result;
};

template<> inline
void StoreQuadIn16(float& data, const Number<__m256>& number, int quad){
    __m128 tmp;
    switch(quad){
    case 0:
    case 2:
        tmp = _mm256_extractf128_ps( number.value, 0x00 );
        break;
    case 1:
    case 3:
        tmp =  _mm256_extractf128_ps( number.value, 0x01 );
        break;
    }
    _mm_store_ps(((&data) + quad*4),tmp);
}

template<> inline
void StoreQuadIn16(float* data, const Number<__m256>& number, int quad){
    __m128 tmp;
    switch(quad){
    case 0:
    case 2:
        tmp = _mm256_extractf128_ps( number.value, 0x00 );
        break;
    case 1:
    case 3:
        tmp =  _mm256_extractf128_ps( number.value, 0x01 );
        break;
    }
    _mm_store_ps(((data) + quad*4),tmp);
}
//==============================================================//
//==============================================================//



#endif
