#ifndef __KERNEL_DIAGONALMATRIX3_H__
#define __KERNEL_DIAGONALMATRIX3_H__


#include "stdlib.h"
#include <type_traits>


template<class Tn> struct Matrix3;


template<class Tn> struct DiagonalMatrix3 {
    typedef Tn TnType;
    typedef typename std::decay<Tn>::type TnBase;
    typedef Matrix3< typename std::add_const< typename std::remove_reference< Tn >::type >::type > ConstNonref_Matrix3;
    typedef Vector3< typename std::add_const< typename std::remove_reference< Tn >::type >::type > ConstNonref_Vector3;
    typedef DiagonalMatrix3< typename std::add_const< typename std::remove_reference< Tn >::type >::type > ConstNonref_DiagMatrix3;
    static const bool IS_REF =  std::is_reference< Tn >::value;
    static const bool IS_CONST = std::is_const< typename std::remove_reference< Tn >::type >::value;

    Tn x11, x22, x33;

    DiagonalMatrix3() :
    x11(), x22(), x33()
    {}

#define CTYPE DiagonalMatrix3
#define INITIALIZER_LIST x11(v.x11), x22(v.x22), x33(v.x33)
    CTYPE( CTYPE< TnBase >& v ):
        INITIALIZER_LIST
    {}
    CTYPE( CTYPE< const TnBase >& v ):
        INITIALIZER_LIST
    {}
    CTYPE( CTYPE< TnBase& >& v ):
        INITIALIZER_LIST
    {}
    CTYPE( CTYPE< const TnBase& >& v ):
        INITIALIZER_LIST
    {}
    CTYPE( const CTYPE< TnBase >& v ):
        INITIALIZER_LIST
    {}
    CTYPE( const CTYPE< const TnBase >& v ):
        INITIALIZER_LIST
    {}
    CTYPE( const CTYPE< TnBase& >& v ):
        INITIALIZER_LIST
    {}
    CTYPE( const CTYPE< const TnBase& >& v ):
        INITIALIZER_LIST
    {}
    template<class MatrixClass>
    explicit CTYPE(MatrixClass& v):
        INITIALIZER_LIST
    {}
#undef INITIALIZER_LIST
#undef CTYPE

    explicit DiagonalMatrix3(const Tn& x11_in, const Tn& x22_in, const Tn& x33_in)
        :x11(x11_in),x22(x22_in),x33(x33_in) {}


    template< class U>
    DiagonalMatrix3& operator=(const U& v)
    {
        static_assert( !IS_CONST,
                       "Implicit assignment operator cannot be used for const-types." );
        x11 = v.x11;
        x22 = v.x22;
        x33 = v.x33;
        return *this;
    }

    template< class U>
    DiagonalMatrix3& operator=(U& v)
    {
        static_assert( !IS_CONST,
                       "Implicit assignment operator cannot be used for const-types." );
        x11 = v.x11;
        x22 = v.x22;
        x33 = v.x33;
        return *this;
    }

    template<class T_DATA>
    void Load_Aligned(const T_DATA (&A)[3])
    {
        static_assert( !IS_CONST,
                       "Load Aligned cannot be used for const-types." );
        x11.Load_Aligned(A[0]);
        x22.Load_Aligned(A[1]);
        x33.Load_Aligned(A[2]);
    }

    template<class T_DATA>
    void Load_Aligned(const T_DATA (&A)[6])
    {
        static_assert( !IS_CONST,
                       "Load Aligned cannot be used for const-types." );
        x11.Load_Aligned(A[0]);
        x22.Load_Aligned(A[3]);
        x33.Load_Aligned(A[5]);
    }

    template<class T_DATA>
    void Load_Aligned(const T_DATA (&A)[9])
    {
        static_assert( !IS_CONST,
                       "Load Aligned cannot be used for const-types." );
        x11.Load_Aligned(A[0]);
        x22.Load_Aligned(A[4]);
        x33.Load_Aligned(A[8]);
    }

    template<class T_DATA>
    void Store(T_DATA (&A)[3]) const
    {
        ::Store(A[0],x11);
        ::Store(A[1],x22);
        ::Store(A[2],x33);
    }

    template<class T_DATA>
    void Store(T_DATA (&A)[6]) const
    {
        ::Store(A[0],x11);
        ::Store(A[3],x22);
        ::Store(A[5],x33);
    }

    template<class T_DATA>
    void Store(T_DATA (&A)[9]) const
    {
        ::Store(A[0],x11);
        ::Store(A[4],x22);
        ::Store(A[8],x33);
    }

    // Accessors

    Tn& operator()(int i, int j){
        if( i<1 || i > 3 || j<1 || j>3 || i!=j)
            exit(1);
        if( i==1 && j==1 )
            return x11;
        if( i==2 && j==2 )
            return x22;
        if( i==3 && j==3 )
            return x33;
    }

    const Tn& operator()(int i, int j) const {
        if( i<1 || i > 3 || j<1 || j>3 || i!=j)
            exit(1);
        if( i==1 && j==1 )
            return x11;
        if( i==2 && j==2 )
            return x22;
        if( i==3 && j==3 )
            return x33;
    }

    // Scalar Math

    DiagonalMatrix3& operator+=( const Tn& n ){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Scalar plus-equals operator (DiagonalMatrix3) cannot be used for const-types." );
        x11 = x11+n; x22 = x22+n; x33 = x33+n;
        return *this;
    }

    DiagonalMatrix3& operator*=( const Tn& n ){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Scalar plus-equals operator (DiagonalMatrix3) cannot be used for const-types." );
        x11 = x11*n; x22 = x22*n; x33 = x33*n;
        return *this;
    }

    DiagonalMatrix3& operator-=( const Tn& n ){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Scalar plus-equals operator (DiagonalMatrix3) cannot be used for const-types." );
        x11 = x11-n; x22 = x22-n; x33 = x33-n;
        return *this;
    }

    DiagonalMatrix3& operator/=( const Tn& n ){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Scalar plus-equals operator (DiagonalMatrix3) cannot be used for const-types." );
        x11 = x11/n; x22 = x22/n; x33 = x33/n;
        return *this;
    }


    // Matrix Math

    template<class U>
    DiagonalMatrix3& operator+=( const DiagonalMatrix3<U>& v ) {
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Matrix plus-equals operator (DiagonalMatrix3) cannot be used for const-types." );
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "Matrix plus-equals operator (DiagonalMatrix3) must operate on same base type." );
        x11 = x11+v.x11; x22 = x22+v.x22; x33 = x33+v.x33;
        return *this;
    }

    template<class U>
    DiagonalMatrix3& operator-=( const DiagonalMatrix3<U>& v ) {
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Matrix minus-equals operator (DiagonalMatrix3) cannot be used for const-types." );
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "Matrix minus-equals operator (DiagonalMatrix3) must operate on same base type." );
        x11 = x11-v.x11; x22 = x22-v.x22; x33 = x33-v.x33;
        return *this;
    }

    DiagonalMatrix3& Inverse() {
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Matrix Inverse (DiagonalMatrix3) cannot be used for const-types." );

        x11 = x11.inverse(); x22 = x22.inverse(); x33 = x33.inverse();
        return *this;
    }

    template<class U>
    void In_Place_Column_Scale( Matrix3<U>& m ) {
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "In_Place_Column_Scale (DiagonalMatrix3) cannot be used for const-types." );
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "In_Place_Column_Scale (DiagonalMatrix3) must operate on same base type." );
        m.x11 = m.x11 * x11;
        m.x21 = m.x21 * x11;
        m.x31 = m.x31 * x11;
        m.x12 = m.x12 * x22;
        m.x22 = m.x22 * x22;
        m.x32 = m.x32 * x22;
        m.x13 = m.x13 * x33;
        m.x23 = m.x23 * x33;
        m.x33 = m.x33 * x33;
    }

    template<class U>
    void In_Place_Row_Scale( Matrix3<U>& m ) {
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "In_Place_Row_Scale (DiagonalMatrix3) cannot be used for const-types." );
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "In_Place_Row_Scale (DiagonalMatrix3) must operate on same base type." );
        m.x11 = m.x11 * x11;
        m.x12 = m.x12 * x11;
        m.x13 = m.x13 * x11;
        m.x21 = m.x21 * x22;
        m.x22 = m.x22 * x22;
        m.x23 = m.x23 * x22;
        m.x31 = m.x31 * x33;
        m.x32 = m.x32 * x33;
        m.x33 = m.x33 * x33;
    }


    template<class U>
    ConstNonref_Vector3 operator*( const Vector3<U>& v ) {
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "Vector multiply (DiagonalMatrix3) must operate on same base type." );
        return ConstNonref_Vector3(x11*v.x, x22*v.y, x33*v.z);
    }


};


#endif
