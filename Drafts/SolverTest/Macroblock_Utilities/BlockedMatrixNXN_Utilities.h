//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef  __BlockedMatrixNXN_Utilities_h__
#define  __BlockedMatrixNXN_Utilities_h__

#include <ostream>

namespace PhysBAM{
template<class T> class MATRIX_MXN;
}

using namespace PhysBAM;

template<class T>
void Print_Matrix(const MATRIX_MXN<T>& A,std::ostream& out);

namespace BlockedMatrixNXN_Utilities
{
template<class T>
void Print_Matrix(const BlockedMatrixNXN<T>& A,std::ostream& out);

template<class T>
void Allocate(const int block_dim,void *&D_raw,void *&L_raw);

template<class T>
void LoadFromPhysbam(BlockedMatrixNXN<T>& A,const MATRIX_MXN<T>& Aref);
}

#endif
