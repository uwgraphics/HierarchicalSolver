//#####################################################################
//  Copyright (c) 2018 Haixiang Liu.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_NUMBER_AVX512_DOUBLE_H__
#define __KERNEL_NUMBER_AVX512_DOUBLE_H__


//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//


template<> inline
Number<__m512d>::Number()
{value=_mm512_castsi512_pd(_mm512_xor_epi32(_mm512_castpd_si512(value),_mm512_castpd_si512(value)));}



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
Number<__m512d> Number<__m512d>::operator+(const Number& other) const
{Number<__m512d> result;result.value=_mm512_add_pd(value,other.value);return result;}



//------------------------------------//
//           MULTIPLICATION           //
//------------------------------------//


template<> inline
Number<__m512d> Number<__m512d>::operator*(const Number& other) const
{Number<__m512d> result;result.value=_mm512_mul_pd(value,other.value);return result;}


//------------------------------------//
//           SUBTRACTION              //
//------------------------------------//


template<> inline
Number<__m512d> Number<__m512d>::operator-(const Number& other) const
{Number<__m512d> result;result.value=_mm512_sub_pd(value,other.value);return result;}


//------------------------------------//
//            DIVISION                //
//------------------------------------//


template<> inline
Number<__m512d> Number<__m512d>::operator/(const Number& other) const
{Number<__m512d> result;result.value=_mm512_div_pd(value,other.value);return result;}


//==============================================================//
//                                                              //
//                      COMPARISON OPERATIONS                   //
//                                                              //
//==============================================================//

//------------------------------------//
//             LESS THAN              //
//------------------------------------//


template<> inline
typename Number<__m512d>::Mask Number<__m512d>::operator<(const Number& other) const
{
    Mask result;
    result.value=_mm512_cmp_pd_mask(value,other.value,_MM_CMPINT_LT);
    return result;
}


//------------------------------------//
//             GREATER THAN              //
//------------------------------------//


template<> inline
typename Number<__m512d>::Mask Number<__m512d>::operator>(const Number& other) const
{
    Mask result;
    result.value=_mm512_cmp_pd_mask(value,other.value,_MM_CMPINT_GT);
    return result;
}


//------------------------------------//
//             LESS EQUALS              //
//------------------------------------//


template<> inline
typename Number<__m512d>::Mask Number<__m512d>::operator<=(const Number& other) const
{
    Mask result;
    result.value=_mm512_cmp_pd_mask(value,other.value,_MM_CMPINT_LE);
    return result;
}


//------------------------------------//
//             GREATER EQUALS              //
//------------------------------------//



template<> inline
typename Number<__m512d>::Mask Number<__m512d>::operator>=(const Number& other) const
{
    Mask result;
    result.value=_mm512_cmp_pd_mask(value,other.value,_MM_CMPINT_GE);
    return result;
}



//------------------------------------//
//               EQUALS               //
//------------------------------------//


template<> inline
typename Number<__m512d>::Mask Number<__m512d>::operator==(const Number& other) const
{
    Mask result;
    result.value=_mm512_cmp_pd_mask(value,other.value,_MM_CMPINT_EQ);
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
Number<__m512d> Number<__m512d>::operator&(const Number& other) const
{Number<__m512d> result;result.value=_mm512_castsi512_pd(_mm512_and_epi32(_mm512_castpd_si512(value),_mm512_castpd_si512(other.value)));return result;}


//------------------------------------//
//             BITWISE OR             //
//------------------------------------//


template<> inline
Number<__m512d> Number<__m512d>::operator|(const Number& other) const
{Number<__m512d> result;result.value=_mm512_castsi512_pd(_mm512_or_epi32(_mm512_castpd_si512(value),_mm512_castpd_si512(other.value)));return result;}


//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//



template<> inline
Number<__m512d> Number<__m512d>::operator^(const Number& other) const
{Number<__m512d> result;result.value=_mm512_castsi512_pd(_mm512_xor_epi32(_mm512_castpd_si512(value),_mm512_castpd_si512(other.value)));return result;}


//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//


template<> inline
Number<__m512d> Number<__m512d>::operator~() const
{
    Number<__m512d> result;
    doubleConverter X;
    X.i = 0xFFFFFFFFFFFFFFFF;
    BUILD_CONSTANT_8_DOUBLE(__one, X.d);
    __m512d val2=_mm512_load_pd(__one);
    result.value=_mm512_castsi512_pd(_mm512_andnot_epi32(_mm512_castpd_si512(value),_mm512_castpd_si512(val2)));
    return result;
}


//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//



template<> inline
Number<__m512d> Number<__m512d>::andnot(const Number& other) const
{
    Number<__m512d> result;
    __m512i A, B;
    A = _mm512_castpd_si512(value);
    B = _mm512_castpd_si512(other.value);
    result.value=_mm512_castsi512_pd(_mm512_and_epi32( B, _mm512_xor_epi32( B, A )));
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
Number<__m512d> Number<__m512d>::sqrt() const
{Number<__m512d> result;result.value=_mm512_sqrt_pd(value);return result;}


//------------------------------------//
//       RECIPROCAL SQUARE ROOT       //
//------------------------------------//


template<> inline
Number<__m512d> Number<__m512d>::rsqrt() const
{
    Number<__m512d> result;
    result.value=_mm512_rsqrt14_pd(value);
    return result;
}



//------------------------------------//
//                 LOG                //
//------------------------------------//


template<> inline
Number<__m512d> Number<__m512d>::log() const
{
    Number<__m512d> result;
    result.value = _mm512_log_pd(value);
return result;
}



//------------------------------------//
//                 EXP                //
//------------------------------------//


template<> inline
Number<__m512d> Number<__m512d>::exp() const
{
    Number<__m512d> result;
    result.value = _mm512_exp_pd(value);
    return result;
}



//------------------------------------//
//        INVERSE (RECIPROCAL)        //
//------------------------------------//


template<> inline
Number<__m512d> Number<__m512d>::inverse() const
{
   Number<__m512d> result;
   result.value = _mm512_rcp14_pd( value );
   return result;
}



//------------------------------------//
//        Absolulte Value             //
//------------------------------------//


template<> inline
Number<__m512d> Number<__m512d>::abs() const
{
    Number<__m512d> neg,result;
    neg.value=_mm512_sub_pd(neg.value,value);
    result.value=_mm512_max_pd(value,neg.value);
    return result;
}



//------------------------------------//
//               SIGN                 //
//------------------------------------//


template<> inline
Number<__m512d> Number<__m512d>::sign() const
{
    Number<__m512d> zero,poss,result;
    Mask negs;
    BUILD_CONSTANT_8_DOUBLE(__one,1.0);
    __m512d none=_mm512_load_pd(__one);
    none = _mm512_sub_pd(zero.value,none);
    __m512d one=_mm512_load_pd(__one);
    negs.value = _mm512_cmp_pd_mask(value,zero.value,_MM_CMPINT_GE);
    result.value = _mm512_mask_blend_pd(negs.value,none,one);
    
    return result;
}



//------------------------------------//
//               MINIMUM              //
//------------------------------------//


template<> inline
Number<__m512d> min(const Number<__m512d>& A, const Number<__m512d>& B)
{Number<__m512d> result;result.value=_mm512_min_pd(A.value, B.value);return result;}


//------------------------------------//
//               MAXIMUM              //
//------------------------------------//


template<> inline
Number<__m512d> max(const Number<__m512d>& A, const Number<__m512d>& B)
{Number<__m512d> result;result.value=_mm512_max_pd(A.value, B.value);return result;}



//------------------------------------//
//      Masked assignment (blend)     //
//------------------------------------//


template<> inline
Number<__m512d> blend(const typename Number<__m512d>::Mask& mask, const Number<__m512d>& A, const Number<__m512d>& B)
{Number<__m512d> result;result.value=_mm512_mask_blend_pd(mask.value,A.value,B.value);return result;}


//------------------------------------//
//      Masked assignment (mask)     //
//------------------------------------//


template<> inline
Number<__m512d> Number<__m512d>::mask(const Mask& mask) const
{
    Number<__m512d> result;
/*
    result.value=_mm512_castsi512_pd(_mm512_mask_xor_epi32(_mm512_castpd_si512(result.value),
                                                           mask.value,
                                                           _mm512_castpd_si512(value),
                                                           _mm512_castpd_si512(value)
                                                           )
                                                           );
*/
    result.value = _mm512_mask_mov_pd( result.value,  mask.value,  value );

    return result;
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
std:: ostream & operator<<( std:: ostream & os, const Number<__m512d>& number)
{
    double *fp = (double*)&number.value;
    os << "[ " 
       << *(fp+8)
       << " " << *(fp+7)
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
void Number<__m512d>::Load(const double* data)
{value=_mm512_load_pd((void*)data);}




template<> inline
void Number<__m512d>::Load(const double& data)
{value=_mm512_load_pd((void*)&data);}

//------------------------------------//
//              LOADS                 //
//------------------------------------//



template<> inline
void Number<__m512d>::Gather(const double* data,const int* offsets)
{
    __m256i index = _mm256_load_epi32((void*)offsets);
    value = _mm512_i32gather_pd(index, (void*)data, 8 );
}




template<> inline
void Number<__m512d>::Gather(const double* data,const int& offsets)
{
    __m256i index = _mm256_load_epi32((void*)&offsets);
    value = _mm512_i32gather_pd(index, (void*)data, 8 );
}



//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//


template<> inline
void Number<__m512d>::Load_Aligned(const double* data)
{value=_mm512_load_pd((void*)data);}





template<> inline
void Number<__m512d>::Load_Aligned(const double& data)
{value=_mm512_load_pd((void*)&data);}



//------------------------------------//
//             STORES                 //
//------------------------------------//


template<> inline
void Store(double* data,const Number<__m512d>& number)
{_mm512_store_pd(data,number.value);}



template<> inline
void Store(double& data,const Number<__m512d>& number)
{_mm512_store_pd(&data,number.value);}


//==============================================================//
//                                                              //
//                  HORIZONTAL OPERATIONS                       //
//                                                              //
//==============================================================//
#define UNTESTED

#if defined(UNTESTED)
#else
template<> inline
Number<__m512d> Number<__m512d>::Spread(int i){
    Number<__m512d> result;
    __mmask16 mask;
    mask = _mm512_int2mask( 0xFFFF );
    switch(i){
    case 0:
        result.value = _mm512_mask_permute4f128_pd( value, mask, value, (_MM_PERM_ENUM)0x0 );
        break;
    case 1:
        result.value = _mm512_mask_permute4f128_pd( value, mask, value, (_MM_PERM_ENUM)0x55 );
        break;
    case 2:
        result.value = _mm512_mask_permute4f128_pd( value, mask, value, (_MM_PERM_ENUM)0xAA );
        break;
    case 3:
        result.value = _mm512_mask_permute4f128_pd( value, mask, value, (_MM_PERM_ENUM)0xFF );
        break;
    }
    return result;
};
#endif

#if defined(UNTESTED)
#else
template<> inline
Number<__m512d> Number<__m512d>::Distribute(int i){
    Number<__m512d> result, tmp;
    tmp.value = _mm512_swizzle_pd(value, _MM_SWIZ_REG_BADC );
    __mmask16 mask1, mask2;
    mask1 = _mm512_int2mask( 0x33CC );
    result.value = _mm512_mask_blend_pd( mask1, value, tmp.value );
    tmp.value = _mm512_swizzle_pd(result.value, _MM_SWIZ_REG_CDAB );
    mask2 = _mm512_int2mask( 0x5A5A );
    result.value = _mm512_mask_blend_pd( mask2, result.value, tmp.value );
    return result;
};
#endif

#if defined(UNTESTED)
#else
template<> inline
Number<__m512d> Number<__m512d>::SwizzleAdd(int i, const Number& other){
    Number<__m512d> result;
    __mmask16 mask;
    mask = _mm512_int2mask( 0xFFFF );
    switch ( i ){
    case 1:        
        result.value = _mm512_mask_permute4f128_pd( value, mask, value, (_MM_PERM_ENUM)0xB1 ); // ???
        result.value = _mm512_add_pd( result.value, value );
        break;
    case 2:
        result.value = _mm512_mask_permute4f128_pd( value, mask, value, (_MM_PERM_ENUM)0x4E ); // ???
        result.value = _mm512_add_pd( result.value, value );
        break;       
    }; 
    return result;
};
#endif

#if defined(UNTESTED)
#else
template<> inline
Number<__m512d> Number<__m512d>::Horizontal_Quad_Add(){
    Number<__m512d> result;
    Number<__m512d> temp;
    temp.value = _mm512_swizzle_pd(value, _MM_SWIZ_REG_CDAB);
    result.value = _mm512_add_pd( value, temp.value);
    temp.value = _mm512_swizzle_pd(result.value, _MM_SWIZ_REG_BADC);
    result.value = _mm512_add_pd( temp.value, result.value);
    return result;
};
#endif

#if defined(UNTESTED)
#else
template<> inline
Number<__m512d> Number<__m512d>::Quad_Mask(int i){
    Number<__m512d> result;
    Number<__m512d> zero;
    __mmask16 mask;
    mask = _mm512_int2mask( 0x8421 );
    result.value = _mm512_mask_blend_pd( mask, zero.value, value );
    return result;
};
#endif

#if defined(UNTESTED)
#else
template<> inline
void StoreQuadIn16(double& data, const Number<__m512d>& number, int quad){
    __mmask16 mask;
    switch(quad){
    case 0:
        mask = _mm512_int2mask( 0x000F );
        break;
    case 1:
        mask = _mm512_int2mask( 0x00F0 );
        break;
    case 2:
        mask = _mm512_int2mask( 0x0F00 );
        break;
    case 3:
        mask = _mm512_int2mask( 0xF000 );
        break;
    }
    _mm512_mask_store_pd(&data,mask,number.value);
}
#endif

#if defined(UNTESTED)
#else
template<> inline
void StoreQuadIn16(double* data, const Number<__m512d>& number, int quad){
    __mmask16 mask;
    switch(quad){
    case 0:
        mask = _mm512_int2mask( 0x000F );
        break;
    case 1:
        mask = _mm512_int2mask( 0x00F0 );
        break;
    case 2:
        mask = _mm512_int2mask( 0x0F00 );
        break;
    case 3:
        mask = _mm512_int2mask( 0xF000 );
        break;
    }
    _mm512_mask_store_pd(data,mask,number.value);
}
#endif
//==============================================================//
//==============================================================//



#endif
