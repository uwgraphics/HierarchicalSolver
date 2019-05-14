//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_MASK_SSE_H__
#define __KERNEL_MASK_SSE_H__

//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//

template<> inline
Mask<__m128i>::Mask()
{value=_mm_xor_si128(value,value);}


template<> inline
Mask<__m128>::Mask()
{value=_mm_xor_ps(value,value);}


//==============================================================//
//                                                              //
//                      BITWISE  OPERATIONS                     //
//                                                              //
//==============================================================//

//------------------------------------//
//             BITWISE AND            //
//------------------------------------//

template<> inline
Mask<__m128i> Mask<__m128i>::operator&(const Mask& other) const
{Mask<__m128i> result;result.value=_mm_and_si128(value,other.value);return result;}

template<> inline
Mask<__m128> Mask<__m128>::operator&(const Mask& other) const
{Mask<__m128> result;result.value=_mm_and_ps(value,other.value);return result;}

//------------------------------------//
//             BITWISE OR             //
//------------------------------------//

template<> inline
Mask<__m128i> Mask<__m128i>::operator|(const Mask& other) const
{Mask<__m128i> result;result.value=_mm_or_si128(value,other.value);return result;}

template<> inline
Mask<__m128> Mask<__m128>::operator|(const Mask& other) const
{Mask<__m128> result;result.value=_mm_or_ps(value,other.value);return result;}

//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//

template<> inline
Mask<__m128i> Mask<__m128i>::operator^(const Mask& other) const
{Mask<__m128i> result;result.value=_mm_xor_si128(value,other.value);return result;}

template<> inline
Mask<__m128> Mask<__m128>::operator^(const Mask& other) const
{Mask<__m128> result;result.value=_mm_xor_ps(value,other.value);return result;}


//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//

template<> inline
Mask<__m128i> Mask<__m128i>::operator~() const
{
    Mask<__m128i> result;
    const int __one[4] = {0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF};
    __m128i val2=_mm_load_si128((__m128i*)&__one);
    result.value=_mm_andnot_si128(value,val2);
    return result;
}


template<> inline
Mask<__m128> Mask<__m128>::operator~() const
{
    Mask<__m128> result;
    MASK_HELPERS::floatConverter fc;
    fc.i = 0xFFFFFFFF;
    const float __one[4] = {fc.f,fc.f,fc.f,fc.f};
    __m128 val2=_mm_load_ps(__one);
    result.value=_mm_andnot_ps(value,val2);
    return result;
}


//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//

template<> inline
Mask<__m128i> Mask<__m128i>::andnot(const Mask& other) const
{
    Mask<__m128i> result;
    result.value=_mm_andnot_si128(value,other.value);
    return result;
}

template<> inline
Mask<__m128> Mask<__m128>::andnot(const Mask& other) const
{
    Mask<__m128> result;
    result.value=_mm_andnot_ps(value,other.value);
    return result;
}

//==============================================================//
//                                                              //
//                  DEBUG PRINT SUPPORT                         //
//                                                              //
//==============================================================//

template<> inline
std:: ostream & operator<<( std:: ostream & os, const Mask<__m128i>& number)
{
    int *fp = (int*)&number.value;
    os << "[ " << *(fp+3)
       << " " << *(fp+2)
       << " " << *(fp+1)
       << " " << *(fp) << " ]";
    return os;
}

template<> inline
std:: ostream & operator<<( std:: ostream & os, const Mask<__m128>& number)
{
    float *fp = (float*)&number.value;
    os << "[ " << *(fp+3)
       << " " << *(fp+2)
       << " " << *(fp+1)
       << " " << *(fp) << " ]";
    return os;
}

//==============================================================//
//                                                              //
//                  LOADS AND STORES                            //
//                                                              //
//==============================================================//

//------------------------------------//
//              LOADS                 //
//------------------------------------//

template<> inline
void Mask<__m128i>::Load(const Ext_Type* data)
{value=_mm_loadu_si128((__m128i*)data);}

template<> inline
void Mask<__m128i>::Load(const Ext_Type& data)
{value=_mm_loadu_si128((__m128i*)&data);}



template<> inline
void Mask<__m128>::Load(const Ext_Type* data)
{value=_mm_loadu_ps(data);}

template<> inline
void Mask<__m128>::Load(const Ext_Type& data)
{value=_mm_loadu_ps(&data);}

//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//

template<> inline
void Mask<__m128i>::Load_Aligned(const Ext_Type* data)
{value=_mm_load_si128((__m128i*)data);}

template<> inline
void Mask<__m128i>::Load_Aligned(const Ext_Type& data)
{value=_mm_load_si128((__m128i*)&data);}


template<> inline
void Mask<__m128>::Load_Aligned(const Ext_Type* data)
{value=_mm_load_ps(data);}

template<> inline
void Mask<__m128>::Load_Aligned(const Ext_Type& data)
{value=_mm_load_ps(&data);}


//------------------------------------//
//             STORES                 //
//------------------------------------//

template<> inline
void Store(typename Mask<__m128i>::Ext_Type* data,const Mask<__m128i>& number)
{_mm_store_si128((__m128i*)data,number.value);}

template<> inline
void Store(typename Mask<__m128i>::Ext_Type& data,const Mask<__m128i>& number)
{_mm_store_si128((__m128i*)&data,number.value);}



template<> inline
void Store(typename Mask<__m128>::Ext_Type* data,const Mask<__m128>& number)
{_mm_store_ps(data,number.value);}

template<> inline
void Store(typename Mask<__m128>::Ext_Type& data,const Mask<__m128>& number)
{_mm_store_ps(&data,number.value);}


#endif
