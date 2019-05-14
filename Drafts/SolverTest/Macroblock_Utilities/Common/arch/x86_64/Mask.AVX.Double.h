//#####################################################################
//  Copyright (c) 2018 Haixiang Liu.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_MASK_AVX_DOUBLE_H__
#define __KERNEL_MASK_AVX_DOUBLE_H__

//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//

template<> inline
Mask<__m256d>::Mask()
{value=_mm256_xor_pd(value,value);}

//==============================================================//
//                                                              //
//                      BITWISE  OPERATIONS                     //
//                                                              //
//==============================================================//

//------------------------------------//
//             BITWISE AND            //
//------------------------------------//

template<> inline
Mask<__m256d> Mask<__m256d>::operator&(const Mask& other) const
{Mask<__m256d> result;result.value=_mm256_and_pd(value,other.value);return result;}




//------------------------------------//
//             BITWISE OR             //
//------------------------------------//

template<> inline
Mask<__m256d> Mask<__m256d>::operator|(const Mask& other) const
{Mask<__m256d> result;result.value=_mm256_or_pd(value,other.value);return result;}




//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//


template<> inline
Mask<__m256d> Mask<__m256d>::operator^(const Mask& other) const
{Mask<__m256d> result;result.value=_mm256_xor_pd(value,other.value);return result;}



//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//


template<> inline
Mask<__m256d> Mask<__m256d>::operator~() const
{
    Mask<__m256d> result;
    const double __one[4] = {1.0,1.0,1.0,1.0};
    __m256d val2=_mm256_load_pd(__one);
    result.value=_mm256_andnot_pd(value,val2);
    return result;
}


//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//


template<> inline
Mask<__m256d> Mask<__m256d>::andnot(const Mask& other) const
{
    Mask<__m256d> result;
    result.value=_mm256_andnot_pd(value,other.value);
    return result;
}



//==============================================================//
//                                                              //
//                  DEBUG PRINT SUPPORT                         //
//                                                              //
//==============================================================//

template<> inline
std:: ostream & operator<<( std:: ostream & os, const Mask<__m256d>& number)
{
    double *fp = (double*)&number.value;
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
void Mask<__m256d>::Load(const double* data)
{value=_mm256_loadu_pd(data);}

template<> inline
void Mask<__m256d>::Load(const double& data)
{value=_mm256_loadu_pd(&data);}

//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//

template<> inline
void Mask<__m256d>::Load_Aligned(const double* data)
{value=_mm256_load_pd(data);}

template<> inline
void Mask<__m256d>::Load_Aligned(const double& data)
{value=_mm256_load_pd(&data);}


//------------------------------------//
//             STORES                 //
//------------------------------------//
template<> inline
void Store(double* data,const Mask<__m256d>& number)
{_mm256_store_pd(data,number.value);}

template<> inline
void Store(double& data,const Mask<__m256d>& number)
{_mm256_store_pd(&data,number.value);}




#endif
