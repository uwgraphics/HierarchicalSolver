//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef  __BlockedCSCSymmetricMatrix3_Utilities_h__
#define  __BlockedCSCSymmetricMatrix3_Utilities_h__

#include <ostream>

template<class T>
void Print_Matrix(const BlockedCSCSymmetricMatrix3_Reference<T>& A,std::ostream& out);

template<class T>
void Print_Matrix(const BlockedCSCSymmetricMatrix3<T>& A,std::ostream& out);

template<class T_DATA>
void Dense_Cholesky_Allocate(const int block_dim,void *&D_raw,void *&L_raw,int *&offsets,int *&row);

template<class T_DATA>
void Sparse_Cholesky_Initialize(BlockedCSCSymmetricMatrix3<T_DATA> &A,const int seed);

template<class T, class T_DATA>
void LoadFromPhysbam( BlockedCSCSymmetricMatrix3<T_DATA> &A, const BlockedCSCSymmetricMatrix3_Reference<T> &S );

template<class T, class T_DATA>
void StoreToPhysbam( BlockedCSCSymmetricMatrix3_Reference<T> &S, const BlockedCSCSymmetricMatrix3<T_DATA> &A );

template<class T, class T_DATA>
void Sparse_Cholesky_Allocate_From_Reference(void *&D_raw,void *&L_raw,int *&offsets,int *&row,
                                             const BlockedCSCSymmetricMatrix3_Reference<T> &A);
template <class T>
void InitializeBlockedLapacian( BlockedCSCSymmetricMatrix3_Reference<T> &A );

template <class T>
void InitializeBlockedIdentity( BlockedCSCSymmetricMatrix3_Reference<T> &A );

template<class T>
void InitializeBlockedDense( BlockedCSCSymmetricMatrix3_Reference<T> &A, int seed, int size );

template<class T, class T_DATA>
void StoreToPhysbam( SYMMETRIC_MATRIX_NXN<T>& S, const BlockedCSCSymmetricMatrix3<T_DATA> &A );

template<class T>
void StoreToPhysbam( MATRIX_MXN<T>& S, const BlockedCSCSymmetricMatrix3_Reference<T> &A );

template<class T, class T_DATA>
void LoadFromPhysbam( BlockedCSCSymmetricMatrix3<T_DATA> &A, const SYMMETRIC_MATRIX_NXN<T>& S );

#endif
