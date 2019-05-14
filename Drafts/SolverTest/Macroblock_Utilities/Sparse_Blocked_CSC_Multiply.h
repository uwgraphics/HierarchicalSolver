//#####################################################################
// Copyright 2018, Qisi Wang, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Sparse_Blocked_CSC_Multiply_h__
#define __Sparse_Blocked_CSC_Multiply_h__

#include "BlockedCSCSymmetricMatrix3.h"
#include <type_traits>

template<class Tw,class T_DATA,class I_DATA>
    void Sparse_Blocked_CSC_Multiply(BlockedCSCSymmetricMatrix3<T_DATA> &A, int end);
#endif
