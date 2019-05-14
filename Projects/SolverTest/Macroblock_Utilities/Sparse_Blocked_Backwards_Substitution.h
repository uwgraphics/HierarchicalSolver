//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Sparse_Blocked_Backwards_Substitution_h__
#define __Sparse_Blocked_Backwards_Substitution_h__

#include "BlockedCSCSymmetricMatrix3.h"
#include <type_traits>

template<class Tw,class T_DATA,class I_DATA>
void Sparse_Blocked_Backwards_Substitution(const BlockedCSCSymmetricMatrix3<T_DATA> &U,
                                           typename std::remove_extent<T_DATA>::type *b,
                                           bool use_diagonal);

#endif
