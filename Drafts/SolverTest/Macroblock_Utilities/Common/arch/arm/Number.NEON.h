//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_NUMBER_NEON_H__
#define __KERNEL_NUMBER_NEON_H__


//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//


template<> inline
Number<float32x4_t>::Number()
{value=vreinterpretq_f32_u32(veorq_u32(vreinterpretq_u32_f32(value),vreinterpretq_u32_f32(value)));}


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
Number<float32x4_t> Number<float32x4_t>::operator+(const Number& other) const
{Number<float32x4_t> result;result.value=vaddq_f32(value,other.value);return result;}

//------------------------------------//
//           MULTIPLICATION           //
//------------------------------------//


template<> inline
Number<float32x4_t> Number<float32x4_t>::operator*(const Number& other) const
{Number<float32x4_t> result;result.value=vmulq_f32(value,other.value);return result;}

//------------------------------------//
//           SUBTRACTION              //
//------------------------------------//


template<> inline
Number<float32x4_t> Number<float32x4_t>::operator-(const Number& other) const
{Number<float32x4_t> result;result.value=vsubq_f32(value,other.value);return result;}

//------------------------------------//
//            DIVISION                //
//------------------------------------//

template<> inline
Number<float32x4_t> Number<float32x4_t>::operator/(const Number& other) const
{
   Number<float32x4_t> result;
   float32x4_t val3=vrecpeq_f32(other.value);
   float32x4_t val2=vaddq_f32(val3,val3);
   val3=vmulq_f32(val3,val3);
   val3=vmulq_f32(val3,other.value);  
   result.value=vsubq_f32(val2,val3);
   result.value=vmulq_f32(value,result.value);
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
Number<float32x4_t> Number<float32x4_t>::operator<(const Number& other) const
{Number<float32x4_t> result;result.value=vreinterpretq_f32_u32(vcltq_f32(value,other.value));return result;}

//------------------------------------//
//             GREATER THAN              //
//------------------------------------//

template<> inline
Number<float32x4_t> Number<float32x4_t>::operator>(const Number& other) const
{Number<float32x4_t> result;result.value=vreinterpretq_f32_u32(vcgtq_f32(value,other.value));return result;}

//------------------------------------//
//             LESS EQUALS              //
//------------------------------------//

template<> inline
Number<float32x4_t> Number<float32x4_t>::operator<=(const Number& other) const
{Number<float32x4_t> result;result.value=vreinterpretq_f32_u32(vcleq_f32(value,other.value));return result;}

//------------------------------------//
//             GREATER EQUALS              //
//------------------------------------//

template<> inline
Number<float32x4_t> Number<float32x4_t>::operator>=(const Number& other) const
{Number<float32x4_t> result;result.value=vreinterpretq_f32_u32(vcgeq_f32(value,other.value));return result;}

//------------------------------------//
//               EQUALS               //
//------------------------------------//

template<> inline
Number<float32x4_t> Number<float32x4_t>::operator==(const Number& other) const
{Number<float32x4_t> result;result.value=vreinterpretq_f32_u32(vceqq_f32(value,other.value));return result;}

//==============================================================//
//                                                              //
//                      BITWISE  OPERATIONS                     //
//                                                              //
//==============================================================//

//------------------------------------//
//             BITWISE AND            //
//------------------------------------//

template<> inline
Number<float32x4_t> Number<float32x4_t>::operator&(const Number& other) const
{Number<float32x4_t> result;result.value=vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(value),vreinterpretq_u32_f32(other.value)));return result;}

//------------------------------------//
//             BITWISE OR             //
//------------------------------------//

template<> inline
Number<float32x4_t> Number<float32x4_t>::operator|(const Number& other) const
{Number<float32x4_t> result;result.value=vreinterpretq_f32_u32(vorrq_u32(vreinterpretq_u32_f32(value),vreinterpretq_u32_f32(other.value)));return result;}

//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//

template<> inline
Number<float32x4_t> Number<float32x4_t>::operator^(const Number& other) const
{Number<float32x4_t> result;result.value=vreinterpretq_f32_u32(veorq_u32(vreinterpretq_u32_f32(value),vreinterpretq_u32_f32(other.value)));return result;}

//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//

template<> inline
Number<float32x4_t> Number<float32x4_t>::operator~() const
{
    Number<float32x4_t> result;
    result.value=vreinterpretq_f32_u32(vmvnq_u32(vreinterpretq_u32_f32(value)));
    return result;
}

//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//

template<> inline
Number<float32x4_t> Number<float32x4_t>::andnot(const Number& other) const
{
    Number<float32x4_t> result;
    result.value=vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(value),vmvnq_u32(vreinterpretq_u32_f32(other.value))));
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
Number<float32x4_t> Number<float32x4_t>::sqrt() const
{
   Number<float32x4_t> result;
   float32x4_t val3=vrsqrteq_f32(value);
   float32x4_t val2=vaddq_f32(val3,val3);
   val3=vmulq_f32(val3,val3);
   val3=vmulq_f32(val3,value);  
   result.value=vsubq_f32(val2,val3);
   result.value=vmulq_f32(result.value,value);
   return result;
}

//------------------------------------//
//       RECIPROCAL SQUARE ROOT       //
//------------------------------------//

template<> inline
Number<float32x4_t> Number<float32x4_t>::rsqrt() const
{
   Number<float32x4_t> result;
   float32x4_t val3=vrsqrteq_f32(value);
   float32x4_t val2=vaddq_f32(val3,val3);
   val3=vmulq_f32(val3,val3);
   val3=vmulq_f32(val3,value);  
   result.value=vsubq_f32(val2,val3);

   return result;
}


//------------------------------------//
//                 LOG                //
//------------------------------------//

template<> inline
Number<float32x4_t> Number<float32x4_t>::log() const
{Number<float32x4_t> result;/*result.value=_mm_log_ps(value);*/return result;}

//------------------------------------//
//                 EXP                //
//------------------------------------//

template<> inline
Number<float32x4_t> Number<float32x4_t>::exp() const
{Number<float32x4_t> result;/*result.value=_mm_exp_ps(value);*/return result;}

//------------------------------------//
//        INVERSE (RECIPROCAL)        //
//------------------------------------//

template<> inline
Number<float32x4_t> Number<float32x4_t>::inverse() const
{
    // 22-bit precision -- See Intel IA32 Architecture Optimization Manual
   Number<float32x4_t> result;
   float32x4_t val3=vrecpeq_f32(value);
   float32x4_t val2=vaddq_f32(val3,val3);
   val3=vmulq_f32(val3,val3);
   val3=vmulq_f32(val3,value);  
   result.value=vsubq_f32(val2,val3);

   return result;
}

//------------------------------------//
//        Absolulte Value             //
//------------------------------------//

template<> inline
Number<float32x4_t> Number<float32x4_t>::abs() const
{
    Number<float32x4_t> result;
    result.value=vabsq_f32(value);
    return result;
}

//------------------------------------//
//               SIGN                 //
//------------------------------------//

template<> inline
Number<float32x4_t> Number<float32x4_t>::sign() const
{
    Number<float32x4_t> zero,negs,poss,result,temp;
    const float __one[4] = {1.0f,1.0f,1.0f,1.0f};
    float32x4_t none=vld1q_f32(__one);
    none = vsubq_f32(zero.value,none);
    float32x4_t one=vld1q_f32(__one);
    negs.value = vreinterpretq_f32_u32(vcgeq_f32(value,zero.value));
    result.value = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(negs.value),vreinterpretq_u32_f32(one)));
    temp.value = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(negs.value),vmvnq_u32(vreinterpretq_u32_f32(none))));
    result.value = vreinterpretq_f32_u32(vorrq_u32(vreinterpretq_u32_f32(result.value),vreinterpretq_u32_f32(temp.value)));
    return result;
}

//------------------------------------//
//               MINIMUM              //
//------------------------------------//

template<> inline
Number<float32x4_t> min(const Number<float32x4_t>& A, const Number<float32x4_t>& B)
{Number<float32x4_t> result;result.value=vminq_f32(A.value, B.value);return result;}

//------------------------------------//
//               MAXIMUM              //
//------------------------------------//

template<> inline
Number<float32x4_t> max(const Number<float32x4_t>& A, const Number<float32x4_t>& B)
{Number<float32x4_t> result;result.value=vmaxq_f32(A.value, B.value);return result;}

//------------------------------------//
//      Masked assignment (blend)     //
//------------------------------------//

template<> inline
Number<float32x4_t> blend(const Number<float32x4_t>& mask, const Number<float32x4_t>& A, const Number<float32x4_t>& B)
{
    Number<float32x4_t> result;
    Number<float32x4_t> temp;
    result.value = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(mask.value),vreinterpretq_u32_f32(B.value)));
    temp.value = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(mask.value),vmvnq_u32(vreinterpretq_u32_f32(A.value))));
    result.value = vreinterpretq_f32_u32(vorrq_u32(vreinterpretq_u32_f32(result.value),vreinterpretq_u32_f32(temp.value)));
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
std:: ostream & operator<<( std:: ostream & os, const Number<float32x4_t>& number)
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
void Number<float32x4_t>::Load(const float* data)
{value=vld1q_f32(data);}


template<> inline
void Number<float32x4_t>::Load(const float& data)
{value=vld1q_f32(&data);}

//------------------------------------//
//              LOADS                 //
//------------------------------------//

template<> inline
void Number<float32x4_t>::Gather(const float* data,const int* offsets)
{
    /* 
    float32x4_t val1=vld1q_dup_f32(data+offsets[0]);
    float32x4_t val2=vld1q_dup_f32(data+offsets[1]);
    float32x4_t val3=vld1q_dup_f32(data+offsets[2]);
    float32x4_t val4=vld1q_dup_f32(data+offsets[3]);



    float32x4x2_t val12=vtrnq_f32(val1, val2);
    float32x4x2_t val34=vtrnq_f32(val3, val4);
    val1=_mm_insert_ps(val1,val2,0x10);
    val3=_mm_insert_ps(val3,val4,0x10);
    value=_mm_movelh_ps(val1,val3);
    */
}


template<> inline
void Number<float32x4_t>::Gather(const float* data,const int& offsets)
{
    /*
    float32x4_t val1=_mm_load_ss(data+offsets);
    float32x4_t val2=_mm_load_ss(data+(&offsets)[1]);
    float32x4_t val3=_mm_load_ss(data+(&offsets)[2]);
    float32x4_t val4=_mm_load_ss(data+(&offsets)[3]);
    val1=_mm_insert_ps(val1,val2,0x10);
    val3=_mm_insert_ps(val3,val4,0x10);
    value=_mm_movelh_ps(val1,val3);
    */
}

//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//

template<> inline
void Number<float32x4_t>::Load_Aligned(const float* data)
{value=vld1q_f32(data);}


template<> inline
void Number<float32x4_t>::Load_Aligned(const float& data)
{value=vld1q_f32(&data);}

//------------------------------------//
//             STORES                 //
//------------------------------------//

template<> inline
void Store(float* data,const Number<float32x4_t>& number)
{vst1q_f32(data,number.value);}


template<> inline
void Store(float& data,const Number<float32x4_t>& number)
{vst1q_f32(&data,number.value);}

//==============================================================//
//==============================================================//


#endif
