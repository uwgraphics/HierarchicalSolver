//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef  __BlockedCSRMatrix3_h__
#define  __BlockedCSRMatrix3_h__

#include <type_traits>

//################################################
// Class BlockedCSRMatrix3
//################################################
template<class T_DATA>
struct BlockedCSRMatrix3
// Sparse matrix in blocked Compressed Sparse Row format, with Matrix3x3 blocks
{
    int rows;        // Matrix row dimension in blocks
    int columns;     // Matrix column dimension in blocks
    void *E_raw;  // Pointer to blocks - 9x T_DATA elements each
    int *offsets; // First index in E for each of n rows (plus a sentinel value)
    int *column;     // Column index of every block

    typedef typename std::remove_extent<T_DATA>::type T;
    static_assert(std::is_floating_point<T>::value,"Must templatize over floating point type, or array");

    typedef T_DATA (*Matrix3_ptr)[9];
    typedef T_DATA (&Matrix3_type)[9];
    typedef const T_DATA (&Matrix3_const_type)[9];

    BlockedCSRMatrix3 Offset(const int offset) const
    {return BlockedCSRMatrix3{rows,columns,&static_cast<T*>(E_raw)[offset],offsets,column};}

    Matrix3_const_type E(const int index) const
    {return static_cast<const Matrix3_ptr>(E_raw)[index];}

    Matrix3_type E(const int index)
    {return static_cast<const Matrix3_ptr>(E_raw)[index];}

};
//#####################################################################
#endif
