//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_NUMBER_DOUBLE_H__
#define __KERNEL_NUMBER_DOUBLE_H__

//#include <math.h>

#ifdef ENABLE_IO_SUPPORT
#include <iostream>
#endif

// Don't know if needed
namespace {
    typedef union {
        unsigned long long i;
        double d;
    } doubleConverter;
}

//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//

template<> inline
Number<double>::Number()
{value=0.;}

// don't know what to do with this
/* template<> inline */
/* Number<double>::Number(const Mask& mask) */
/* { */
/*     if(mask.value){ */
/*         floatConverter ALL_ONES; */
/*         ALL_ONES.i = (0xFFFFFFFFFFFFFFFF); */
/*         value = ALL_ONES.d; */
/*     } */
/*     else */
/*         value = 0; */
/* } */

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
Number<double> Number<double>::operator+(const Number& other) const
{

    Number result;
    result.value=value+other.value;
    return result;
}

//------------------------------------//
//           MULTIPLICATION           //
//------------------------------------//

template<> inline
Number<double> Number<double>::operator*(const Number& other) const
{Number result;result.value=value*other.value;return result;}


//------------------------------------//
//           SUBTRACTION              //
//------------------------------------//

template<> inline
Number<double> Number<double>::operator-(const Number& other) const
{Number result;result.value=value-other.value;return result;}

//------------------------------------//
//            DIVISION                //
//------------------------------------//

template<> inline
Number<double> Number<double>::operator/(const Number& other) const
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
Number<double>::Mask Number<double>::operator<(const Number& other) const
{
    Mask result;
    result.value = value < other.value ? true : false;
    return result;
}


//------------------------------------//
//             GREATER THAN              //
//------------------------------------//

template<> inline
Number<double>::Mask Number<double>::operator>(const Number& other) const
{
    Mask result;
    result.value = value > other.value ? true : false;
    return result;
}

//------------------------------------//
//             LESS EQUALS              //
//------------------------------------//

template<> inline
Number<double>::Mask Number<double>::operator<=(const Number& other) const
{
    Mask result;
    result.value = value <= other.value ? true : false;
    return result;
}


//------------------------------------//
//             GREATER EQUALS              //
//------------------------------------//

template<> inline
Number<double>::Mask Number<double>::operator>=(const Number& other) const
{
    Mask result;
    result.value = value >= other.value ? true : false;
    return result;
}


//------------------------------------//
//               EQUALS               //
//------------------------------------//

template<> inline
Number<double>::Mask Number<double>::operator==(const Number& other) const
{
    Mask result;
    result.value = value == other.value ? true : false;
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
Number<double> Number<double>::operator&(const Number& other) const
{
    Number result;
    doubleConverter VAL1;
    doubleConverter VAL2;
    doubleConverter RESULT;
    VAL1.d = value;
    VAL2.d = other.value;
    RESULT.i = VAL1.i & VAL2.i;
    result.value = RESULT.d;
    return result;
}


//------------------------------------//
//             BITWISE OR             //
//------------------------------------//

template<> inline
Number<double> Number<double>::operator|(const Number& other) const
{
    Number result;
    doubleConverter VAL1;
    doubleConverter VAL2;
    doubleConverter RESULT;
    VAL1.d = value;
    VAL2.d = other.value;
    RESULT.i = VAL1.i | VAL2.i;
    result.value = RESULT.d;
    return result;
}


//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//

template<> inline
Number<double> Number<double>::operator^(const Number& other) const
{
    Number result;
    doubleConverter VAL1;
    doubleConverter VAL2;
    doubleConverter RESULT;
    VAL1.d = value;
    VAL2.d = other.value;
    RESULT.i = VAL1.i ^ VAL2.i;
    result.value = RESULT.d;
    return result;
}


//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//

template<> inline
Number<double> Number<double>::operator~() const
{
    Number result;
    doubleConverter ALL_ONES;
    ALL_ONES.i = (0xFFFFFFFFFFFFFFFF);
    doubleConverter VAL1;
    doubleConverter RESULT;
    VAL1.d = value;
    RESULT.i = ~VAL1.i & ALL_ONES.i;
    result.value = RESULT.d;
    return result;
}

//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//

template<> inline
Number<double> Number<double>::andnot(const Number& other) const
{
    Number result;
    doubleConverter VAL1;
    doubleConverter VAL2;
    doubleConverter RESULT;
    VAL1.d = value;
    VAL2.d = other.value;
    RESULT.i = ~VAL1.i & VAL2.i;
    result.value = RESULT.d;
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
Number<double> Number<double>::sqrt() const
{Number result;result.value=::sqrt(value);return result;}


//------------------------------------//
//       RECIPROCAL SQUARE ROOT       //
//------------------------------------//

template<> inline
Number<double> Number<double>::rsqrt() const
{
Number result;
 result.value=1.0/::sqrt(value);
 return result;
}

//------------------------------------//
//                 LOG                //
//------------------------------------//

template<> inline
Number<double> Number<double>::log() const
{
    Number result;
    result.value=::log(value);
    return result;
}



//------------------------------------//
//                 EXP                //
//------------------------------------//

template<> inline
Number<double> Number<double>::exp() const
{
    Number result;
    result.value=::exp(value);
    return result;
}


//------------------------------------//
//        INVERSE (RECIPROCAL)        //
//------------------------------------//

template<> inline
Number<double> Number<double>::inverse() const
{
    Number result;
    result.value=1./value;
    return result;
}



//------------------------------------//
//        Absolulte Value             //
//------------------------------------//

template<> inline
Number<double> Number<double>::abs() const
{
    Number result;
    result.value=(value>=0.)?value:-value;
    return result;
}

//------------------------------------//
//               SIGN                 //
//------------------------------------//

template<> inline
Number<double> Number<double>::sign() const
{
    Number result;
    result.value=(value<0.)?-1.0f:1.0f;
    return result;
}

//------------------------------------//
//               MINIMUM              //
//------------------------------------//

template<> inline
Number<double> min(const Number<double>& A, const Number<double>& B)
{Number<double> result;result.value=A.value < B.value ? A.value : B.value;return result;}

//------------------------------------//
//               MAXIMUM              //
//------------------------------------//

template<> inline
Number<double> max(const Number<double>& A, const Number<double>& B)
{Number<double> result;result.value=A.value > B.value ? A.value : B.value;return result;}


//------------------------------------//
//      Masked assignment (blend)     //
//------------------------------------//

template<> inline
Number<double> blend(const NumberPolicy<Number<double> >::MASK_TYPE& mask, const Number<double>& A, const Number<double>& B)
{
    Number<double> result;
    result.value = mask.value ? B.value : A.value;
    return result;
}

//------------------------------------//
//      Masked assignment (mask )     //
//------------------------------------//
template<> inline
Number<double> Number<double>::mask(const Mask& mask) const
{
    Number<double> result;
    result.value = mask.value ? value : (double)(0.0);
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
std:: ostream & operator<<( std:: ostream & os, const Number<double>& number)
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
void Number<double>::Load(const double* data)
{value=*data;}


template<> inline
void Number<double>::Load(const double& data)
{value=data;}


//------------------------------------//
//              LOADS                 //
//------------------------------------//

template<> inline
void Number<double>::Gather(const double* data,const int* offsets)
{value=data[*offsets];}

template<> inline
void Number<double>::Gather(const double* data,const int& offsets)
{value=data[offsets];}


//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//

template<> inline
void Number<double>::Load_Aligned(const double* data)
{value=*data;}

template<> inline
void Number<double>::Load_Aligned(const double& data)
{value=data;}

//------------------------------------//
//             STORES                 //
//------------------------------------//

template<> inline
void Store(double* data,const Number<double>& number)
{*data=number.value;}

template<> inline
void Store(double& data,const Number<double>& number)
{data=number.value;}
//==============================================================//
//                                                              //
//                  HORIZONTAL OPERATIONS                       //
//                                                              //
//==============================================================//

template<> inline
Number<double> Number<double>::Spread(int i){
    return *this;
};


template<> inline
Number<double> Number<double>::Distribute(int i){
    return *this;
};


template<> inline
Number<double> Number<double>::SwizzleAdd(int i, const Number& other){
    return *this;
};

template<> inline
Number<double> Number<double>::Horizontal_Quad_Add(){
    return *this;
};

template<> inline
Number<double> Number<double>::Quad_Mask(int i){
    return *this;
};

// don't understand these two
template<> inline
void StoreQuadIn16(double& data, const Number<double>& number, int quad){
    data=number.value;
}

template<> inline
void StoreQuadIn16(double* data, const Number<double>& number, int quad){
    *data=number.value;
}

//==============================================================//
//==============================================================//

#endif
