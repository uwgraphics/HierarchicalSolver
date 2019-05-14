//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef  __BlockedMatrixNXN_h__
#define  __BlockedMatrixNXN_h__

#include <type_traits>

//################################################
// Class BlockedMatrixNXN
//################################################
template<class T>
struct BlockedMatrixNXN
// Dense symmetric matrix in block column major order, with 4x4 matrix blocks
{
    typedef T (*Block_ptr)[16];
    typedef T (&Block_type)[16];
    typedef const T (&Block_const_type)[16];

    int n;       // Matrix dimension (in 4-blocks)
    void *D_raw; // Pointer to (n) diagonal blocks - 16x elements each (of type T)
    void *L_raw; // Pointer to (n*(n-1)/2) strictly lower diagonal blocks - 16x elements each (of type T)

    Block_type D(const int col) 
    {return static_cast<Block_ptr>(D_raw)[col];}

    Block_const_type D(const int col) const
    {return static_cast<const Block_ptr>(D_raw)[col];}

    Block_type L(const int index) 
    {return static_cast<Block_ptr>(L_raw)[index];}

    Block_const_type L(const int index) const
    {return static_cast<const Block_ptr>(L_raw)[index];}
};
//#####################################################################
#endif
