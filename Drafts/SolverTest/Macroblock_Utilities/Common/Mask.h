//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_MASK_H__
#define __KERNEL_MASK_H__

#include <cmath>
#include <iostream>

template<class Tw> class Mask;
template<class Tw> void Store(typename MaskPolicy<Mask<Tw> >::MASK_EXTERNAL_TYPE* data,const Mask<Tw>& number);
template<class Tw> void Store(typename MaskPolicy<Mask<Tw> >::MASK_EXTERNAL_TYPE& data,const Mask<Tw>& number);
template<class Tw> std::ostream& operator<<( std::ostream& os, const Mask<Tw>& number);

namespace MASK_HELPERS{
    typedef union {
        int i;
        float f;
    } floatConverter;
}

template<class Tw>
class Mask
{
public:
    Tw value;
    typedef typename MaskPolicy<Mask<Tw> >::MASK_EXTERNAL_TYPE Ext_Type;

    Mask();

    static Mask True();
    static Mask False();

    Mask operator&(const Mask& other) const;
    Mask operator|(const Mask& other) const;
    Mask operator^(const Mask& other) const;
    Mask andnot(const Mask& other) const;
    Mask operator~() const;

    void Load_Aligned(const Ext_Type* data);
    void Load_Aligned(const Ext_Type& data);

    void Load(const Ext_Type* data);
    void Load(const Ext_Type& data);

    friend void Store<>(Ext_Type* data,const Mask& number);
    friend void Store<>(Ext_Type& data,const Mask& number);

    friend std::ostream& operator<< <>( std:: ostream & os, const Mask& number);
};

//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//

template<> inline
Mask<bool>::Mask()
{value=false;}

template<> inline
Mask<bool> Mask<bool>::True()
{Mask result; result.value = true; return result;}

template<> inline
Mask<bool> Mask<bool>::False()
{Mask result; result.value = false; return result;}

//==============================================================//
//                                                              //
//                      BITWISE  OPERATIONS                     //
//                                                              //
//==============================================================//

//------------------------------------//
//             BITWISE AND            //
//------------------------------------//

template<> inline
Mask<bool> Mask<bool>::operator&(const Mask& other) const
{
    Mask<bool> result;
    result.value = value && other.value;
    return result;
}

//------------------------------------//
//             BITWISE OR             //
//------------------------------------//

template<> inline
Mask<bool> Mask<bool>::operator|(const Mask& other) const
{
    Mask<bool> result;
    result.value = value || other.value;
    return result;
}


//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//

template<> inline
Mask<bool> Mask<bool>::operator^(const Mask& other) const
{
    Mask<bool> result;
    result.value  = (bool)(((int)value) ^ ((int)other.value));
    return result;
}


//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//

template<> inline
Mask<bool> Mask<bool>::operator~() const
{
    Mask<bool> result;
    result.value = !value;
    return result;
}


//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//

template<> inline
Mask<bool> Mask<bool>::andnot(const Mask& other) const
{
    Mask<bool> result;
    result.value = value && (!other.value);
    return result;
}

//==============================================================//
//                                                              //
//                  DEBUG PRINT SUPPORT                         //
//                                                              //
//==============================================================//

template<> inline
std:: ostream & operator<<( std:: ostream & os, const Mask<bool>& number)
{
    os << "[ " << number.value << " ]";
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
void Mask<bool>::Load(const Ext_Type* data)
{value=*data;}

template<> inline
void Mask<bool>::Load(const Ext_Type& data)
{value=data;}

//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//

template<> inline
void Mask<bool>::Load_Aligned(const Ext_Type* data)
{value=*data;}

template<> inline
void Mask<bool>::Load_Aligned(const Ext_Type& data)
{value=data;}

//------------------------------------//
//             STORES                 //
//------------------------------------//

template<> inline
void Store(Mask<bool>::Ext_Type* data,const Mask<bool>& number)
{
    MASK_HELPERS::floatConverter ALL_ONES;
    ALL_ONES.i = (0xFFFFFFFF);
    if(number.value)
        *data=ALL_ONES.f;
    else
        *data=0.0f;
}

template<> inline
void Store(Mask<bool>::Ext_Type& data,const Mask<bool>& number)
{
    MASK_HELPERS::floatConverter ALL_ONES;
    ALL_ONES.i = (0xFFFFFFFF);
    if(number.value)
        data=ALL_ONES.f;
    else
        data=0.0f;
}

#endif
