//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_CONSTANTS_H__
#define __KERNEL_CONSTANTS_H__

// Or any other 16-wide arch
#if defined(ENABLE_MIC_INSTRUCTION_SET) || defined(ENABLE_AVX512_INSTRUCTION_SET)
#define BUILD_CONSTANT(name,value) const float name[16] __attribute__((aligned(64))) = {value,value,value,value, \
                                                                                      value,value,value,value, \
                                                                                      value,value,value,value, \
                                                                                      value,value,value,value};
#define BUILD_ICONSTANT(name,value) const int name[16] __attribute__((aligned(64))) = {value,value,value,value, \
                                                                                      value,value,value,value, \
                                                                                      value,value,value,value, \
                                                                                      value,value,value,value};
#define BUILD_CONSTANT_DOUBLE(name,value) const double name[8] __attribute__((aligned(64))) = {value,value,value,value, \
                                                                                               value,value,value,value};
#define BUILD_ICONSTANT_DOUBLE(name,value) const int64_t name[8] __attribute__((aligned(64))) = {value,value,value,value, \
                                                                                                 value,value,value,value};
#define BUILD_TDATA(name,size) T_DATA name size __attribute__((aligned(64)));
#define BUILD_IDATA(name,size) I_DATA name size __attribute__((aligned(64)));
#else

// Or any other 8-wide arch
#if defined(ENABLE_AVX_INSTRUCTION_SET)  
#define BUILD_CONSTANT(name,value) const float name[8] __attribute__((aligned(32))) = {value,value,value,value, \
                                                                                     value,value,value,value};
#define BUILD_ICONSTANT(name,value) const int name[8] __attribute__((aligned(32))) = {value,value,value,value, \
                                                                                     value,value,value,value};
#define BUILD_CONSTANT_DOUBLE(name,value) const double name[8] __attribute__((aligned(32))) = {value,value,value,value};
#define BUILD_ICONSTANT_DOUBLE(name,value) const int64_t name[8] __attribute__((aligned(32))) = {value,value,value,value};
#define BUILD_TDATA(name,size) T_DATA name size __attribute__((aligned(32)));
#define BUILD_IDATA(name,size) I_DATA name size __attribute__((aligned(32)));
#else

// Or any other 4-wide arch
#if defined(ENABLE_NEON_INSTRUCTION_SET) || defined(ENABLE_SSE_INSTRUCTION_SET) 
#define BUILD_CONSTANT(name,value) const float name[4] __attribute__((aligned(16))) = {value,value,value,value};
#define BUILD_ICONSTANT(name,value) const int name[4] __attribute__((aligned(16))) = {value,value,value,value};
#define BUILD_CONSTANT_DOUBLE(name,value) const double name[2] __attribute__((aligned(16))) = {value,value};
#define BUILD_ICONSTANT_DOUBLE(name,value) const int64_t name[2] __attribute__((aligned(16))) = {value,value};
#define BUILD_TDATA(name,size) T_DATA name size __attribute__((aligned(16)));
#define BUILD_IDATA(name,size) I_DATA name size __attribute__((aligned(16)));
#else

// Fallback if no vector arch's are enabled.
#define BUILD_CONSTANT(name,value) const float name[1] __attribute__((aligned(4))) = {value};
#define BUILD_ICONSTANT(name,value) const int name[1] __attribute__((aligned(4))) = {value};
#define BUILD_CONSTANT_DOUBLE(name,value) const double name[1] __attribute__((aligned(4))) = {value};
#define BUILD_ICONSTANT_DOUBLE(name,value) const int64_t name[1] __attribute__((aligned(4))) = {value};
#define BUILD_TDATA(name,size) T_DATA name size __attribute__((aligned(4)));
#define BUILD_IDATA(name,size) I_DATA name size __attribute__((aligned(4)));

#endif // End Size 4 group
#endif // End Size 8 group
#endif // End Size 16 group

#define BUILD_CONSTANT_16(name,value) const float name[16] __attribute__((aligned(64))) = { \
        value,value,value,value,                                        \
        value,value,value,value,                                        \
        value,value,value,value,                                        \
        value,value,value,value};

#define BUILD_CONSTANT_8(name,value) const float name[8] __attribute__((aligned(32))) = { \
        value,value,value,value,                                        \
        value,value,value,value};

#define BUILD_CONSTANT_4(name,value) const float name[4] __attribute__((aligned(16))) = { \
        value,value,value,value};

#define BUILD_CONSTANT_1(name,value) const float name[1] __attribute__((aligned(4))) = { \
        value};

#define BUILD_ICONSTANT_16(name,value) const int name[16] __attribute__((aligned(64))) = { \
        value,value,value,value,                                        \
        value,value,value,value,                                        \
        value,value,value,value,                                        \
        value,value,value,value};

#define BUILD_ICONSTANT_8(name,value) const int name[8] __attribute__((aligned(32))) = { \
        value,value,value,value,                                        \
        value,value,value,value};

#define BUILD_ICONSTANT_4(name,value) const int name[4] __attribute__((aligned(16))) = { \
        value,value,value,value};

#define BUILD_ICONSTANT_1(name,value) const int name[1] __attribute__((aligned(4))) = { \
        value};

#define BUILD_CONSTANT_16_DOUBLE(name,value) const double name[16] __attribute__((aligned(64))) = { \
        value,value,value,value,                                        \
        value,value,value,value,                                        \
        value,value,value,value,                                        \
        value,value,value,value};

#define BUILD_CONSTANT_8_DOUBLE(name,value) const double name[8] __attribute__((aligned(32))) = { \
        value,value,value,value,                                        \
        value,value,value,value};

#define BUILD_CONSTANT_4_DOUBLE(name,value) const double name[4] __attribute__((aligned(16))) = { \
        value,value,value,value};

#define BUILD_CONSTANT_1_DOUBLE(name,value) const double name[1] __attribute__((aligned(4))) = { \
        value};

#define BUILD_ICONSTANT_16_DOUBLE(name,value) const int64_t name[16] __attribute__((aligned(64))) = { \
        value,value,value,value,                                        \
        value,value,value,value,                                        \
        value,value,value,value,                                        \
        value,value,value,value};

#define BUILD_ICONSTANT_8_DOUBLE(name,value) const int64_t name[8] __attribute__((aligned(32))) = { \
        value,value,value,value,                                        \
        value,value,value,value};

#define BUILD_ICONSTANT_4_DOUBLE(name,value) const int64_t name[4] __attribute__((aligned(16))) = { \
        value,value,value,value};

#define BUILD_ICONSTANT_1_DOUBLE(name,value) const int64_t name[1] __attribute__((aligned(4))) = { \
        value};



#endif
