//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Sparse_Blocked_Forwards_Substitution_h__
#define __Sparse_Blocked_Forwards_Substitution_h__

#include "BlockedCSCSymmetricMatrix3.h"
#include "BlockedCSRMatrix3.h"
#include <type_traits>

template<class Tw,class T_DATA,class I_DATA>
void Sparse_Blocked_Forwards_Substitution(const BlockedCSCSymmetricMatrix3<T_DATA> &U,
                                          typename std::remove_extent<T_DATA>::type *b,
                                          bool use_diagonal);

template<class Tw,class T_DATA,class I_DATA>
void Sparse_Blocked_Forwards_Substitution_3x3(const BlockedCSCSymmetricMatrix3<T_DATA> &U,
                                              typename std::remove_extent<T_DATA>::type *B,
                                              bool use_diagonal);

template<class Tw,class T_DATA,class I_DATA>
    void Sparse_Blocked_Forwards_Substitution_Matrix(const BlockedCSCSymmetricMatrix3<T_DATA> &U,
                                                     const BlockedCSRMatrix3<T_DATA> &A,
                                                     BlockedCSRMatrix3<T_DATA> &B);

#endif
