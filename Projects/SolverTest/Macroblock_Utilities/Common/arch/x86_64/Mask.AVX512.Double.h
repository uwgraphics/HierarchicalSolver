//#####################################################################
//  Copyright (c) 2018 Haixiang Liu.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_MASK_AVX512_DOUBLE_H__
#define __KERNEL_MASK_AVX512_DOUBLE_H__

#include <cmath>
#include <iostream>

//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//

template<> inline
Mask<__mmask8>::Mask()
{value=(__mmask8)_mm512_kxor((__mmask16)value, (__mmask16)value);}

template<> inline
Mask<__mmask8> Mask<__mmask8>::True()
{Mask result; result.value = (__mmask8)_mm512_knot((__mmask16)result.value); return result;}

template<> inline
Mask<__mmask8> Mask<__mmask8>::False()
{Mask result; return result;}

//==============================================================//
//                                                              //
//                      BITWISE  OPERATIONS                     //
//                                                              //
//==============================================================//

//------------------------------------//
//             BITWISE AND            //
//------------------------------------//

template<> inline
Mask<__mmask8> Mask<__mmask8>::operator&(const Mask& other) const
{
    Mask<__mmask8> result;
    result.value = (__mmask8)_mm512_kand((__mmask16)value,other.value);
    return result;
}

//------------------------------------//
//             BITWISE OR             //
//------------------------------------//

template<> inline
Mask<__mmask8> Mask<__mmask8>::operator|(const Mask& other) const
{
    Mask<__mmask8> result;
    result.value = (__mmask8)_mm512_kor((__mmask16)value, (__mmask16)other.value);
    return result;
}


//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//

template<> inline
Mask<__mmask8> Mask<__mmask8>::operator^(const Mask& other) const
{
    Mask<__mmask8> result;
    result.value  = (__mmask8)_mm512_kxor((__mmask16)value, (__mmask16)other.value);
    return result;
}


//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//

template<> inline
Mask<__mmask8> Mask<__mmask8>::operator~() const
{
    Mask<__mmask8> result;
    result.value = (__mmask8)_mm512_knot((__mmask16)value);
    return result;
}


//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//

template<> inline
Mask<__mmask8> Mask<__mmask8>::andnot(const Mask& other) const
{
    Mask<__mmask8> result;
    result.value = (__mmask8)_mm512_kandn((__mmask16)value, (__mmask16)other.value);
    return result;
}

//==============================================================//
//                                                              //
//                  DEBUG PRINT SUPPORT                         //
//                                                              //
//==============================================================//

template<> inline
std:: ostream & operator<<( std:: ostream & os, const Mask<__mmask8>& number)
{
    int mask = _mm512_mask2int( (__mmask16)number.value );
    os << "[ ";
    os << (int)(mask & 0x1) << " ";
    os << (int)(mask & 0x2) << " ";
    os << (int)(mask & 0x4) << " ";
    os << (int)(mask & 0x8) << " ";
    os << (int)(mask & 0x10) << " ";
    os << (int)(mask & 0x20) << " ";
    os << (int)(mask & 0x40) << " ";
    os << (int)(mask & 0x80) << " ";
    os << " ]";
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
void Mask<__mmask8>::Load(const Ext_Type* data)
{value=(__mmask8)_mm512_int2mask(*data);}

template<> inline
void Mask<__mmask8>::Load(const Ext_Type& data)
{value=(__mmask8)_mm512_int2mask(data);}

//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//

template<> inline
void Mask<__mmask8>::Load_Aligned(const Ext_Type* data)
{value=(__mmask8)_mm512_int2mask(*data);}

template<> inline
void Mask<__mmask8>::Load_Aligned(const Ext_Type& data)
{value=(__mmask8)_mm512_int2mask(data);}

//------------------------------------//
//             STORES                 //
//------------------------------------//

template<> inline
void Store(typename Mask<__mmask8>::Ext_Type* data,const Mask<__mmask8>& number)
{*data=_mm512_mask2int((__mmask16)number.value);}

template<> inline
void Store(typename Mask<__mmask8>::Ext_Type& data,const Mask<__mmask8>& number)
{data=_mm512_mask2int((__mmask16)number.value);}

#endif
