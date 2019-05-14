//#####################################################################
// Copyright 2012-2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Sparse_Cholesky_h__
#define __Sparse_Cholesky_h__

#include "BlockedCSCSymmetricMatrix3.h"

template<class T_RAW,class T_DATA,class I_DATA>
    void Sparse_Cholesky(BlockedCSCSymmetricMatrix3<T_DATA> &A);

template<class T_RAW,class T_DATA,class I_DATA>
    void Sparse_Cholesky_Partial(BlockedCSCSymmetricMatrix3<T_DATA> &A, const int stop, const int start=0);

#endif
