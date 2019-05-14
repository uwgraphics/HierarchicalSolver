//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __BlockedMatrixNXN_Cholesky_h__
#define __BlockedMatrixNXN_Cholesky_h__

#include "BlockedMatrixNXN.h"

template<class T>
void Block_Cholesky(BlockedMatrixNXN<T> &A);

template<class T>
void Block_Cholesky_From_Scalar(BlockedMatrixNXN<T> &A);

template<class Tw,class T>
void Block_Forward_Substitution(const BlockedMatrixNXN<T> &A,T* x,const bool use_diagonal);

template<class Tw,class T>
void Block_Backward_Substitution(const BlockedMatrixNXN<T> &A,T* x,const bool use_diagonal);

#endif
