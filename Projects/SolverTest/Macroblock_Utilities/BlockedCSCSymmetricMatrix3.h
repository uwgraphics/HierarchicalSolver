//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef  __BlockedCSCSymmetricMatrix3_h__
#define  __BlockedCSCSymmetricMatrix3_h__

#include <type_traits>

//################################################
// Class BlockedCSCSymmetricMatrix3
//################################################
template<class T_DATA>
struct BlockedCSCSymmetricMatrix3
// Sparse symmetric matrix in blocked Compressed Sparse Column format, with Matrix3x3 blocks
{
    int n;        // Matrix dimension in blocks
    void *D_raw;  // Pointer to (n) diagonal blocks - 6x T_DATA elements each
    void *L_raw;  // Pointer to strictly lower diagonal blocks - 9x T_DATA elements each
    int *offsets; // First index in L for each of n columns (plus a sentinel value)
    int *row;     // Row index of every lower triangular block

    typedef typename std::remove_extent<T_DATA>::type T;
    static_assert(std::is_floating_point<T>::value,"Must templatize over floating point type, or array");
    const static int T_SIZE = std::rank<T_DATA>::value == 0 ? 1 : std::extent<T_DATA>::value;

    typedef T_DATA (*Lower_Triangular_Matrix3_ptr)[6];
    typedef T_DATA (&Lower_Triangular_Matrix3_type)[6];
    typedef const T_DATA (&Lower_Triangular_Matrix3_const_type)[6];

    typedef T_DATA (*Matrix3_ptr)[9];
    typedef T_DATA (&Matrix3_type)[9];
    typedef const T_DATA (&Matrix3_const_type)[9];

    BlockedCSCSymmetricMatrix3 Offset(const int offset) const
    // For a T_DATA != float, moves the matrix down one channel of T_DATA
    {return BlockedCSCSymmetricMatrix3{n,&static_cast<T*>(D_raw)[offset],&static_cast<T*>(L_raw)[offset],offsets,row};}

    BlockedCSCSymmetricMatrix3 SubMatrix(const int col ) const
    // Return a new sub-matrix which starts (col)s down from the orginal
    {return BlockedCSCSymmetricMatrix3{n-col,D_raw,L_raw,offsets,row};}

    int Row(const int i) const
    {return row[i];}

    int Offsets(const int i) const
    {return offsets[i];}

    Lower_Triangular_Matrix3_const_type D(const int col) const
    {return static_cast<const Lower_Triangular_Matrix3_ptr>(D_raw)[col];}

    Lower_Triangular_Matrix3_type D(const int col)
    {return static_cast<const Lower_Triangular_Matrix3_ptr>(D_raw)[col];}

    Matrix3_const_type L(const int index) const
    {return static_cast<const Matrix3_ptr>(L_raw)[index];}

    Matrix3_type L(const int index)
    {return static_cast<const Matrix3_ptr>(L_raw)[index];}
};
//#####################################################################
#endif
