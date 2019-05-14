//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Sparse_Blocked_CSR_Multiply_Metaprogram_h__
#define __Sparse_Blocked_CSR_Multiply_Metaprogram_h__

#include "BlockedCSRMatrix3.h"
#include <type_traits>

template<class Tw,class T_DATA,class I_DATA>
void Sparse_Blocked_CSR_Multiply_Metaprogram(const BlockedCSRMatrix3<T_DATA> &U,
                                             const typename std::remove_extent<T_DATA>::type *b,
                                             typename std::remove_extent<T_DATA>::type *x,
                                             bool transpose);


template<class Tw,class T_DATA,class I_DATA>
void Sparse_Blocked_CSR_Masked_Multiply_Metaprogram(const BlockedCSRMatrix3<T_DATA> &U,
                                                    const typename std::remove_extent<T_DATA>::type *b,
                                                    typename std::remove_extent<T_DATA>::type *x,
                                                    const T_DATA (&mask),
                                                    bool transpose);

#endif
