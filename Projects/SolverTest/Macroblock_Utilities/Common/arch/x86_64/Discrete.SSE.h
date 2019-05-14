//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_DISCRETE_SSE_H__
#define __KERNEL_DISCRETE_SSE_H__


//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//

template<> inline
Discrete<__m128i>::Discrete()
{value=_mm_xor_si128(value,value);}

template<> inline
Discrete<__m128i>::Discrete(const Mask& mask)
{
    Discrete A; A.value = _mm_set1_epi32( 0 );
    Discrete B; B.value = _mm_set1_epi32( 0xFFFFFFFF );
    value = _mm_castps_si128( _mm_blendv_ps( _mm_castsi128_ps(A.value), _mm_castsi128_ps(B.value), mask.value ));
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
Discrete<__m128i> Discrete<__m128i>::operator+(const Discrete& other) const
{Discrete<__m128i> result;result.value=_mm_add_epi32(value,other.value);return result;}

//------------------------------------//
//           MULTIPLICATION           //
//------------------------------------//

template<> inline
Discrete<__m128i> Discrete<__m128i>::operator*(const Discrete& other) const
{Discrete<__m128i> result;result.value=_mm_mul_epi32(value,other.value);return result;}



//------------------------------------//
//           SUBTRACTION              //
//------------------------------------//

template<> inline
Discrete<__m128i> Discrete<__m128i>::operator-(const Discrete& other) const
{Discrete<__m128i> result;result.value=_mm_sub_epi32(value,other.value);return result;}


//------------------------------------//
//            DIVISION                //
//------------------------------------//

template<> inline
Discrete<__m128i> Discrete<__m128i>::operator/(const Discrete& other) const
{
    Discrete<__m128i> result;
#ifdef __INTEL_COMPILER
    result.value=_mm_div_epi32(value,other.value);
#else
    __m128 this_ps = _mm_cvtepi32_ps(value);
    __m128 other_ps = _mm_cvtepi32_ps(other.value);
    __m128 result_ps = _mm_div_ps(this_ps,other_ps);
    result.value = _mm_cvttps_epi32(result_ps);    
#endif

    return result;
}



//==============================================================//
//                                                              //
//                      COMPARISON OPERATIONS                   //
//                                                              //
//==============================================================//

//------------------------------------//
//             LESS THAN              //
//------------------------------------//

template<> inline
Discrete<__m128i>::Mask Discrete<__m128i>::operator<(const Discrete& other) const
{
    Mask result;result.value=_mm_castsi128_ps(_mm_cmplt_epi32(value,other.value));return result;
}


//------------------------------------//
//             GREATER THAN              //
//------------------------------------//

template<> inline
Discrete<__m128i>::Mask Discrete<__m128i>::operator>(const Discrete& other) const
{
    Mask result;result.value=_mm_castsi128_ps(_mm_cmpgt_epi32(value,other.value));return result;
}

//------------------------------------//
//             LESS EQUALS              //
//------------------------------------//

template<> inline
Discrete<__m128i>::Mask Discrete<__m128i>::operator<=(const Discrete& other) const
{
    Mask result;
    result.value=_mm_castsi128_ps( _mm_or_si128(
                                                _mm_cmplt_epi32(value,other.value),
                                                _mm_cmpeq_epi32(value,other.value)
                                                )
                                   );
    return result;
}


//------------------------------------//
//             GREATER EQUALS              //
//------------------------------------//

template<> inline
Discrete<__m128i>::Mask Discrete<__m128i>::operator>=(const Discrete& other) const
{
    Mask result;
    result.value=_mm_castsi128_ps( _mm_or_si128(
                                                _mm_cmpgt_epi32(value,other.value),
                                                _mm_cmpeq_epi32(value,other.value)
                                                )
                                   );
    return result;
}


//------------------------------------//
//               EQUALS               //
//------------------------------------//

template<> inline
Discrete<__m128i>::Mask Discrete<__m128i>::operator==(const Discrete& other) const
{
    Mask result;result.value=_mm_castsi128_ps(_mm_cmpeq_epi32(value,other.value));return result;
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
Discrete<__m128i> Discrete<__m128i>::operator&(const Discrete& other) const
{
    Discrete result;
    result.value = _mm_and_si128( value, other.value );
    return result;
}


//------------------------------------//
//             BITWISE OR             //
//------------------------------------//

template<> inline
Discrete<__m128i> Discrete<__m128i>::operator|(const Discrete& other) const
{
    Discrete result;
    result.value = _mm_or_si128( value, other.value );
    return result;
}


//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//

template<> inline
Discrete<__m128i> Discrete<__m128i>::operator^(const Discrete& other) const
{
    Discrete result;
    result.value = _mm_xor_si128( value, other.value );
    return result;
}


//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//

template<> inline
Discrete<__m128i> Discrete<__m128i>::operator~() const
{
    Discrete result;
    BUILD_ICONSTANT_4(__one, 0xFFFFFFFF);
    __m128i val2=_mm_load_si128((const __m128i*)(__one));
    result.value=_mm_xor_si128(value,val2);
    return result;
}

//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//

template<> inline
Discrete<__m128i> Discrete<__m128i>::andnot(const Discrete& other) const
{
    Discrete result;
    result.value = _mm_andnot_si128( value, other.value );
    return result;
}


//==============================================================//
//==============================================================//

//------------------------------------//
//               MINIMUM              //
//------------------------------------//

template<> inline
Discrete<__m128i> min(const Discrete<__m128i>& A, const Discrete<__m128i>& B) 
{Discrete<__m128i> result;result.value=_mm_min_epi32(A.value, B.value);return result;}

//------------------------------------//
//               MAXIMUM              //
//------------------------------------//

template<> inline
Discrete<__m128i> max(const Discrete<__m128i>& A, const Discrete<__m128i>& B)
{Discrete<__m128i> result;result.value=_mm_min_epi32(A.value, B.value);return result;}


//------------------------------------//
//      Masked assignment (blend)     //
//------------------------------------//

template<> inline
Discrete<__m128i> blend(const NumberPolicy<Discrete<__m128i> >::MASK_TYPE& mask, const Discrete<__m128i>& A, const Discrete<__m128i>& B)
{
    Discrete<__m128i> result;
    result.value = _mm_castps_si128( _mm_blendv_ps( _mm_castsi128_ps(A.value), _mm_castsi128_ps(B.value), mask.value ));
    return result;
}

//------------------------------------//
//      Masked assignment (mask )     //
//------------------------------------//
template<> inline
Discrete<__m128i> Discrete<__m128i>::mask(const Mask& mask) const
{
    Discrete<__m128i> result;
    Discrete<__m128i> temp;
    result.value = _mm_castps_si128( _mm_blendv_ps( _mm_castsi128_ps(temp.value), _mm_castsi128_ps(value), mask.value ));
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
std:: ostream & operator<<( std:: ostream & os, const Discrete<__m128i>& discrete)
{
    int *fp = (int*)&discrete.value;
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
void Discrete<__m128i>::Load(const int* data)
{value=_mm_loadu_si128((const __m128i*)(data));}


template<> inline
void Discrete<__m128i>::Load(const int& data)
{value=_mm_loadu_si128((const __m128i*)(&data));}


//------------------------------------//
//              LOADS                 //
//------------------------------------//

template<> inline
void Discrete<__m128i>::Gather(const int* data,const int* offsets)
{
    // Cheap hacks here: treat int* as float* for loading and bit movement
    __m128 val1=_mm_load_ss(((float*)data)+offsets[0]);
    __m128 val2=_mm_load_ss(((float*)data)+offsets[1]);
    __m128 val3=_mm_load_ss(((float*)data)+offsets[2]);
    __m128 val4=_mm_load_ss(((float*)data)+offsets[3]);
    val1=_mm_insert_ps(val1,val2,0x10);
    val3=_mm_insert_ps(val3,val4,0x10);
    value=_mm_castps_si128( _mm_movelh_ps(val1,val3));
}

template<> inline
void Discrete<__m128i>::Gather(const int* data,const int& offsets)
{
    // Cheap hacks here: treat int* as float* for loading and bit movement
    __m128 val1=_mm_load_ss(((float*)data)+offsets);
    __m128 val2=_mm_load_ss(((float*)data)+(&offsets)[1]);
    __m128 val3=_mm_load_ss(((float*)data)+(&offsets)[2]);
    __m128 val4=_mm_load_ss(((float*)data)+(&offsets)[3]);
    val1=_mm_insert_ps(val1,val2,0x10);
    val3=_mm_insert_ps(val3,val4,0x10);
    value=_mm_castps_si128(_mm_movelh_ps(val1,val3));
}


//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//

template<> inline
void Discrete<__m128i>::Load_Aligned(const int* data)
{value=_mm_load_si128((const __m128i*)(data));}

template<> inline
void Discrete<__m128i>::Load_Aligned(const int& data)
{value=_mm_load_si128((const __m128i*)(&data));}

//------------------------------------//
//             STORES                 //
//------------------------------------//

template<> inline
void Store(int* data,const Discrete<__m128i>& discrete)
{_mm_store_si128((__m128i*)(data),discrete.value);}

template<> inline
void Store(int& data,const Discrete<__m128i>& discrete)
{_mm_store_si128((__m128i*)(&data),discrete.value);}


//==============================================================//
//==============================================================//
#if 0

#endif

#endif
