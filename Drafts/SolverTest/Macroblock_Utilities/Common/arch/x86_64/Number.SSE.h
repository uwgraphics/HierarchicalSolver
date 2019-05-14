//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_NUMBER_SSE_H__
#define __KERNEL_NUMBER_SSE_H__


//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//


template<> inline
Number<__m128>::Number()
{value=_mm_xor_ps(value,value);}


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
Number<__m128> Number<__m128>::operator+(const Number& other) const
{Number<__m128> result;result.value=_mm_add_ps(value,other.value);return result;}

//------------------------------------//
//           MULTIPLICATION           //
//------------------------------------//


template<> inline
Number<__m128> Number<__m128>::operator*(const Number& other) const
{Number<__m128> result;result.value=_mm_mul_ps(value,other.value);return result;}

//------------------------------------//
//           SUBTRACTION              //
//------------------------------------//


template<> inline
Number<__m128> Number<__m128>::operator-(const Number& other) const
{Number<__m128> result;result.value=_mm_sub_ps(value,other.value);return result;}

//------------------------------------//
//            DIVISION                //
//------------------------------------//

template<> inline
Number<__m128> Number<__m128>::operator/(const Number& other) const
{Number<__m128> result;result.value=_mm_div_ps(value,other.value);return result;}

//==============================================================//
//                                                              //
//                      COMPARISON OPERATIONS                   //
//                                                              //
//==============================================================//

//------------------------------------//
//             LESS THAN              //
//------------------------------------//

template<> inline
typename Number<__m128>::Mask Number<__m128>::operator<(const Number& other) const
{Mask result;result.value=_mm_cmplt_ps(value,other.value);return result;}

//------------------------------------//
//             GREATER THAN              //
//------------------------------------//

template<> inline
typename Number<__m128>::Mask Number<__m128>::operator>(const Number& other) const
{Mask result;result.value=_mm_cmpgt_ps(value,other.value);return result;}

//------------------------------------//
//             LESS EQUALS              //
//------------------------------------//

template<> inline
typename Number<__m128>::Mask Number<__m128>::operator<=(const Number& other) const
{Mask result;result.value=_mm_cmple_ps(value,other.value);return result;}

//------------------------------------//
//             GREATER EQUALS              //
//------------------------------------//

template<> inline
typename Number<__m128>::Mask Number<__m128>::operator>=(const Number& other) const
{Mask result;result.value=_mm_cmpge_ps(value,other.value);return result;}

//------------------------------------//
//               EQUALS               //
//------------------------------------//

template<> inline
typename Number<__m128>::Mask Number<__m128>::operator==(const Number& other) const
{Mask result;result.value=_mm_cmpeq_ps(value,other.value);return result;}

//==============================================================//
//                                                              //
//                      BITWISE  OPERATIONS                     //
//                                                              //
//==============================================================//

//------------------------------------//
//             BITWISE AND            //
//------------------------------------//

template<> inline
Number<__m128> Number<__m128>::operator&(const Number& other) const
{Number<__m128> result;result.value=_mm_and_ps(value,other.value);return result;}

//------------------------------------//
//             BITWISE OR             //
//------------------------------------//

template<> inline
Number<__m128> Number<__m128>::operator|(const Number& other) const
{Number<__m128> result;result.value=_mm_or_ps(value,other.value);return result;}

//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//

template<> inline
Number<__m128> Number<__m128>::operator^(const Number& other) const
{Number<__m128> result;result.value=_mm_xor_ps(value,other.value);return result;}

//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//

template<> inline
Number<__m128> Number<__m128>::operator~() const
{
    Number<__m128> result;
    floatConverter X;
    X.i = 0xFFFFFFFF;
    BUILD_CONSTANT_4(__one, X.f);
    __m128 val2=_mm_load_ps(__one);
    result.value=_mm_andnot_ps(value,val2);
    return result;
}

//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//

template<> inline
Number<__m128> Number<__m128>::andnot(const Number& other) const
{
    Number<__m128> result;
    result.value=_mm_andnot_ps(value,other.value);
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
Number<__m128> Number<__m128>::sqrt() const
{Number<__m128> result;result.value=_mm_sqrt_ps(value);return result;}

//------------------------------------//
//       RECIPROCAL SQUARE ROOT       //
//------------------------------------//

template<> inline
Number<__m128> Number<__m128>::rsqrt() const
{Number<__m128> result;result.value=_mm_rsqrt_ps(value);return result;}

//------------------------------------//
//                 LOG                //
//------------------------------------//


template<> inline
Number<__m128> Number<__m128>::log() const
{
    Number<__m128> result;
#ifdef __INTEL_COMPILER
    result.value=_mm_log_ps(value);    
#else

#endif
    return result;
}

//------------------------------------//
//                 EXP                //
//------------------------------------//


template<> inline
Number<__m128> Number<__m128>::exp() const
{
    Number<__m128> result;
#ifdef __INTEL_COMPILER
    result.value=_mm_exp_ps(value);
#else

#endif
    return result;
}


//------------------------------------//
//        INVERSE (RECIPROCAL)        //
//------------------------------------//

template<> inline
Number<__m128> Number<__m128>::inverse() const
{
    // 22-bit precision -- See Intel IA32 Architecture Optimization Manual
   Number<__m128> result;
   __m128 val3=_mm_rcp_ps(value);
   __m128 val2=_mm_add_ps(val3,val3);
   val3=_mm_mul_ps(val3,val3);
   val3=_mm_mul_ps(val3,value);  
   result.value=_mm_sub_ps(val2,val3);

   return result;
}

//------------------------------------//
//        Absolulte Value             //
//------------------------------------//

template<> inline
Number<__m128> Number<__m128>::abs() const
{
    Number<__m128> neg,result;
    neg.value=_mm_sub_ps(neg.value,value);
    result.value=_mm_max_ps(value,neg.value);
    return result;
}

//------------------------------------//
//               SIGN                 //
//------------------------------------//

template<> inline
Number<__m128> Number<__m128>::sign() const
{
    Number<__m128> zero,negs,poss,result;
    BUILD_CONSTANT_4(__one, 1.0f);
    __m128 none=_mm_load_ps(__one);
    none = _mm_sub_ps(zero.value,none);
    __m128 one=_mm_load_ps(__one);
    negs.value = _mm_cmp_ps(value,zero.value,_CMP_GE_OS);
    result.value = _mm_blendv_ps(none,one,negs.value);
    return result;
}

//------------------------------------//
//               MINIMUM              //
//------------------------------------//

template<> inline
Number<__m128> min(const Number<__m128>& A, const Number<__m128>& B)
{Number<__m128> result;result.value=_mm_min_ps(A.value, B.value);return result;}

//------------------------------------//
//               MAXIMUM              //
//------------------------------------//

template<> inline
Number<__m128> max(const Number<__m128>& A, const Number<__m128>& B)
{Number<__m128> result;result.value=_mm_max_ps(A.value, B.value);return result;}

//------------------------------------//
//      Masked assignment (blend)     //
//------------------------------------//

template<> inline
Number<__m128> blend(const typename Number<__m128>::Mask& mask, const Number<__m128>& A, const Number<__m128>& B)
{
    Number<__m128> result;
    Number<__m128> temp;
    result.value=_mm_blendv_ps(A.value,B.value,mask.value);
//    result.value = _mm_and_ps(mask.value, B.value);
//    temp.value = _mm_andnot_ps(mask.value, A.value);
//    result.value = _mm_or_ps(result.value, temp.value);
    return result;
}

//------------------------------------//
//      Masked assignment (mask)     //
//------------------------------------//

template<> inline
Number<__m128> Number<__m128>::mask(const Mask& mask) const
{
    Number<__m128> result;
    Number<__m128> temp;
    
    result.value=_mm_blendv_ps(temp.value,value,mask.value);
//    result.value = _mm_and_ps(mask.value, B.value);
//    temp.value = _mm_andnot_ps(mask.value, A.value);
//    result.value = _mm_or_ps(result.value, temp.value);
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
std:: ostream & operator<<( std:: ostream & os, const Number<__m128>& number)
{
    float *fp = (float*)&number.value;
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
void Number<__m128>::Load(const float* data)
{value=_mm_loadu_ps(data);}


template<> inline
void Number<__m128>::Load(const float& data)
{value=_mm_loadu_ps(&data);}

//------------------------------------//
//              LOADS                 //
//------------------------------------//

template<> inline
void Number<__m128>::Gather(const float* data,const int* offsets)
{
    __m128 val1=_mm_load_ss(data+offsets[0]);
    __m128 val2=_mm_load_ss(data+offsets[1]);
    __m128 val3=_mm_load_ss(data+offsets[2]);
    __m128 val4=_mm_load_ss(data+offsets[3]);
    val1=_mm_insert_ps(val1,val2,0x10);
    val3=_mm_insert_ps(val3,val4,0x10);
    value=_mm_movelh_ps(val1,val3);
}


template<> inline
void Number<__m128>::Gather(const float* data,const int& offsets)
{
    __m128 val1=_mm_load_ss(data+offsets);
    __m128 val2=_mm_load_ss(data+(&offsets)[1]);
    __m128 val3=_mm_load_ss(data+(&offsets)[2]);
    __m128 val4=_mm_load_ss(data+(&offsets)[3]);
    val1=_mm_insert_ps(val1,val2,0x10);
    val3=_mm_insert_ps(val3,val4,0x10);
    value=_mm_movelh_ps(val1,val3);
}

//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//

template<> inline
void Number<__m128>::Load_Aligned(const float* data)
{value=_mm_load_ps(data);}


template<> inline
void Number<__m128>::Load_Aligned(const float& data)
{value=_mm_load_ps(&data);}

//------------------------------------//
//             STORES                 //
//------------------------------------//

template<> inline
void Store(float* data,const Number<__m128>& number)
{_mm_store_ps(data,number.value);}


template<> inline
void Store(float& data,const Number<__m128>& number)
{_mm_store_ps(&data,number.value);}
//==============================================================//
//                                                              //
//                  HORIZONTAL OPERATIONS                       //
//                                                              //
//==============================================================//

template<> inline
Number<__m128> Number<__m128>::Spread(int i){
    return *this;
};


template<> inline
Number<__m128> Number<__m128>::Distribute(int i){
    Number<__m128> result;
    switch( i ){
    case 0:
        result.value = _mm_castsi128_ps( _mm_shuffle_epi32( _mm_castps_si128( value ), 0x00 ) );
        break;
    case 1:
        result.value = _mm_castsi128_ps( _mm_shuffle_epi32( _mm_castps_si128( value ), 0x55 ) );
        break;
    case 2:
        result.value = _mm_castsi128_ps( _mm_shuffle_epi32( _mm_castps_si128( value ), 0xAA ) );
        break;
    case 3:
        result.value = _mm_castsi128_ps( _mm_shuffle_epi32( _mm_castps_si128( value ), 0xFF ) );
        break;
    }
    return result;
};


template<> inline
Number<__m128> Number<__m128>::SwizzleAdd(int i, const Number& other){
    return *this;
};

template<> inline
Number<__m128> Number<__m128>::Horizontal_Quad_Add(){
    Number<__m128> result;
    result.value = _mm_hadd_ps( value, value );
    result.value = _mm_hadd_ps( result.value, result.value );
    return result;
};

template<> inline
Number<__m128> Number<__m128>::Quad_Mask(int i){
    Number<__m128> result;
    Number<__m128> zero;
    switch( i ) {
        case 0: 
            result.value = _mm_blend_ps( zero.value, value, 0x01);break;
        case 1: 
            result.value = _mm_blend_ps( zero.value, value, 0x02);break;
        case 2: 
            result.value = _mm_blend_ps( zero.value, value, 0x04);break;
        case 3: 
            result.value = _mm_blend_ps( zero.value, value, 0x08);break;
    }
    return result;
};



template<> inline
void StoreQuadIn16(float& data, const Number<__m128>& number, int quad){
    _mm_store_ps(((&data) + quad*4),number.value);
}

template<> inline
void StoreQuadIn16(float* data, const Number<__m128>& number, int quad){
    _mm_store_ps(((data) + quad*4),number.value);
}

//==============================================================//
//==============================================================//



#endif
