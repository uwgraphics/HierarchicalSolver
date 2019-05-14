#ifndef __KERNEL_MATRIX3_H__
#define __KERNEL_MATRIX3_H__

#include "stdlib.h"
#include <type_traits>

#include "Common/KernelCommon.h"
#include "Common/Vector3.h"

template<class Tn> struct DiagonalMatrix3;
template<class Tn> struct SymmetricMatrix3;
template<class Tn> struct UnitriangularMatrix3;
//template<class Tn> struct Vector3;

template<class Tn> struct Matrix3 {
    typedef Tn TnType;
    typedef typename std::decay<Tn>::type TnBase;
    typedef Matrix3< typename std::add_const< typename std::remove_reference< Tn >::type >::type > ConstNonref_Matrix3;
    typedef Vector3< typename std::add_const< typename std::remove_reference< Tn >::type >::type > ConstNonref_Vector3;
    typedef SymmetricMatrix3< typename std::add_const< typename std::remove_reference< Tn >::type >::type > ConstNonref_SymMatrix3;
    typedef UnitriangularMatrix3< typename std::add_const< typename std::remove_reference< Tn >::type >::type > ConstNonref_UniMatrix3;
    typedef DiagonalMatrix3< typename std::add_const< typename std::remove_reference< Tn >::type >::type > ConstNonref_DiagMatrix3;
    static const bool IS_REF =  std::is_reference< Tn >::value;
    static const bool IS_CONST = std::is_const< typename std::remove_reference< Tn >::type >::value;


    Tn x11, x21, x31, x12, x22, x32, x13, x23, x33;

    Matrix3() :
        x11(), x21(), x31(),
        x12(), x22(), x32(),
        x13(), x23(), x33()
    {}

#define CTYPE Matrix3
#define INITIALIZER_LIST x11(v.x11), x21(v.x21), x31(v.x31), \
                         x12(v.x12), x22(v.x22), x32(v.x32), \
                         x13(v.x13), x23(v.x23), x33(v.x33)
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

    explicit Matrix3(const Tn& x11_in,const Tn& x21_in,const Tn& x31_in,
                     const Tn& x12_in,const Tn& x22_in,const Tn& x32_in,
                     const Tn& x13_in,const Tn& x23_in,const Tn& x33_in )
        :x11(x11_in),x21(x21_in),x31(x31_in),
         x12(x12_in),x22(x22_in),x32(x32_in),
         x13(x13_in),x23(x23_in),x33(x33_in)
    {}


    template< class U>
    Matrix3& operator=( U& v)
    {
        static_assert( !IS_CONST,
                       "Implicit assignment operator cannot be used for const-types." );

        x11 = v.x11;
        x21 = v.x21;
        x31 = v.x31;
        x12 = v.x12;
        x22 = v.x22;
        x32 = v.x32;
        x13 = v.x13;
        x23 = v.x23;
        x33 = v.x33;
        return *this;
    }

    template< class U>
    Matrix3& operator=(const U& v)
    {
        static_assert( !IS_CONST,
                       "Implicit assignment operator cannot be used for const-types." );

        x11 = v.x11;
        x21 = v.x21;
        x31 = v.x31;
        x12 = v.x12;
        x22 = v.x22;
        x32 = v.x32;
        x13 = v.x13;
        x23 = v.x23;
        x33 = v.x33;
        return *this;
    }

    template<class T_DATA>
    void Load_Aligned(const T_DATA (&A)[9])
    {
        static_assert( !IS_CONST,
                       "Load Aligned cannot be used for const-types." );

        x11.Load_Aligned(A[0]);
        x21.Load_Aligned(A[1]);
        x31.Load_Aligned(A[2]);
        x12.Load_Aligned(A[3]);
        x22.Load_Aligned(A[4]);
        x32.Load_Aligned(A[5]);
        x13.Load_Aligned(A[6]);
        x23.Load_Aligned(A[7]);
        x33.Load_Aligned(A[8]);
    }

    template<class T_DATA>
    void Load_Aligned_Prefetch(const T_DATA (&A)[9], int offset)
    {
        static_assert( !IS_CONST,
                       "Load Aligned cannot be used for const-types." );


        enum {prefetch_offset=9};


        x11.Load_Aligned(A[0]); _mm_prefetch((char*)&A[0+prefetch_offset], _MM_HINT_NTA);
        x21.Load_Aligned(A[1]); _mm_prefetch((char*)&A[1+prefetch_offset], _MM_HINT_NTA);
        x31.Load_Aligned(A[2]); _mm_prefetch((char*)&A[2+prefetch_offset], _MM_HINT_NTA);
        x12.Load_Aligned(A[3]); _mm_prefetch((char*)&A[3+prefetch_offset], _MM_HINT_NTA);
        x22.Load_Aligned(A[4]); _mm_prefetch((char*)&A[4+prefetch_offset], _MM_HINT_NTA);
        x32.Load_Aligned(A[5]); _mm_prefetch((char*)&A[5+prefetch_offset], _MM_HINT_NTA);
        x13.Load_Aligned(A[6]); _mm_prefetch((char*)&A[6+prefetch_offset], _MM_HINT_NTA);
        x23.Load_Aligned(A[7]); _mm_prefetch((char*)&A[7+prefetch_offset], _MM_HINT_NTA);
        x33.Load_Aligned(A[8]); _mm_prefetch((char*)&A[8+prefetch_offset], _MM_HINT_NTA);


    }


    template<class T_DATA>
    void Store(T_DATA (&A)[9]) const
    {
        ::Store(A[0],x11);
        ::Store(A[1],x21);
        ::Store(A[2],x31);
        ::Store(A[3],x12);
        ::Store(A[4],x22);
        ::Store(A[5],x32);
        ::Store(A[6],x13);
        ::Store(A[7],x23);
        ::Store(A[8],x33);
    }

    // Accessors

    Tn& operator()(int i, int j){
        if( i<1 || i > 3 || j<1 || j>3)
            exit(1);
        if( i==1 && j==1 )
            return x11;
        if( i==1 && j==2 )
            return x12;
        if( i==1 && j==3 )
            return x13;
        if( i==2 && j==1 )
            return x21;
        if( i==2 && j==2 )
            return x22;
        if( i==2 && j==3 )
            return x23;
        if( i==3 && j==1 )
            return x31;
        if( i==3 && j==2 )
            return x32;
        if( i==3 && j==3 )
            return x33;
    }

    const Tn& operator()(int i, int j) const {
        if( i<1 || i > 3 || j<1 || j>3)
            exit(1);
        if( i==1 && j==1 )
            return x11;
        if( i==1 && j==2 )
            return x12;
        if( i==1 && j==3 )
            return x13;
        if( i==2 && j==1 )
            return x21;
        if( i==2 && j==2 )
            return x22;
        if( i==2 && j==3 )
            return x23;
        if( i==3 && j==1 )
            return x31;
        if( i==3 && j==2 )
            return x32;
        if( i==3 && j==3 )
            return x33;
    }

    // Scalar Math

    Matrix3& operator+=( const Tn& n ){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Scalar plus-equals operator (Matrix3) cannot be used for const-types." );
        x11 = x11+n; x21 = x21+n; x31 = x31+n;
        x12 = x12+n; x22 = x22+n; x32 = x32+n;
        x13 = x13+n; x23 = x23+n; x33 = x33+n;
        return *this;
    }

    Matrix3& operator*=( const Tn& n ){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Scalar times-equals operator (Matrix3) cannot be used for const-types." );
        x11 = x11*n; x21 = x21*n; x31 = x31*n;
        x12 = x12*n; x22 = x22*n; x32 = x32*n;
        x13 = x13*n; x23 = x23*n; x33 = x33*n;
        return *this;
    }

    Matrix3& operator-=( const Tn& n ){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Scalar minus-equals operator (Matrix3) cannot be used for const-types." );
        x11 = x11-n; x21 = x21-n; x31 = x31-n;
        x12 = x12-n; x22 = x22-n; x32 = x32-n;
        x13 = x13-n; x23 = x23-n; x33 = x33-n;
        return *this;
    }

    Matrix3& operator/=( const Tn& n ){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Scalar divide-equals operator (Matrix3) cannot be used for const-types." );
        x11 = x11/n; x21 = x21/n; x31 = x31/n;
        x12 = x12/n; x22 = x22/n; x32 = x32/n;
        x13 = x13/n; x23 = x23/n; x33 = x33/n;
        return *this;
    }


    // Matrix Math

    template<class U>
    Matrix3& operator+=( const Matrix3<U>& v ) {
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Matrix plus-equals operator (Matrix3) cannot be used for const-types." );
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "Matrix plus-equals operator (Matrix3) must operate on same base type." );
        x11 = x11+v.x11; x21 = x21+v.x21; x31 = x31+v.x31;
        x12 = x12+v.x12; x22 = x22+v.x22; x32 = x32+v.x32;
        x13 = x13+v.x13; x23 = x23+v.x23; x33 = x33+v.x33;
        return *this;
    }

    template<class U>
    Matrix3& operator-=( const Matrix3<U>& v ) {
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Matrix minus-equals operator (Matrix3) cannot be used for const-types." );
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "Matrix minus-equals operator (Matrix3) must operate on same base type." );
        x11 = x11-v.x11; x21 = x21-v.x21; x31 = x31-v.x31;
        x12 = x12-v.x12; x22 = x22-v.x22; x32 = x32-v.x32;
        x13 = x13-v.x13; x23 = x23-v.x23; x33 = x33-v.x33;
        return *this;
    }

    template<class U>
    Matrix3& operator*=( const Matrix3<U>& v ) {
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Matrix times-equals operator (Matrix3) cannot be used for const-types." );
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "Matrix times-equals operator (Matrix3) must operate on same base type." );
        TnBase n11 = x11*v.x11 + x12*v.x21 + x13*v.x31;
        TnBase n21 = x21*v.x11 + x22*v.x21 + x23*v.x31;
        TnBase n31 = x31*v.x11 + x32*v.x21 + x33*v.x31;

        TnBase n12 = x11*v.x12 + x12*v.x22 + x13*v.x32;
        TnBase n22 = x21*v.x12 + x22*v.x22 + x23*v.x32;
        TnBase n32 = x31*v.x12 + x32*v.x22 + x33*v.x32;

        TnBase n13 = x11*v.x13 + x12*v.x23 + x13*v.x33;
        TnBase n23 = x21*v.x13 + x22*v.x23 + x23*v.x33;
        TnBase n33 = x31*v.x13 + x32*v.x23 + x33*v.x33;
        x11=n11;
        x21=n21;
        x31=n31;
        x12=n12;
        x22=n22;
        x32=n32;
        x13=n13;
        x23=n23;
        x33=n33;
        return *this;
    }

    template<class U>
    ConstNonref_Matrix3 operator*( const Matrix3<U>& v ) {
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "Matrix times operator (Matrix3) must operate on same base type." );
        return ConstNonref_Matrix3( x11*v.x11 + x12*v.x21 + x13*v.x31,
                                    x21*v.x11 + x22*v.x21 + x23*v.x31,
                                    x31*v.x11 + x32*v.x21 + x33*v.x31,
                                    x11*v.x12 + x12*v.x22 + x13*v.x32,
                                    x21*v.x12 + x22*v.x22 + x23*v.x32,
                                    x31*v.x12 + x32*v.x22 + x33*v.x32,
                                    x11*v.x13 + x12*v.x23 + x13*v.x33,
                                    x21*v.x13 + x22*v.x23 + x23*v.x33,
                                    x31*v.x13 + x32*v.x23 + x33*v.x33);
    }

    template<class U> inline
    void In_Place_Times_Transpose( const Matrix3<U>& v ) {
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "In Place Times Transpose (Matrix3) must operate on same base type." );

        Vector3<TnBase> C1;
        Vector3<TnBase> C2;
        Vector3<TnBase> C3;
        {
            Vector3<TnBase&> A( x11, x21, x31 );
            C1 = A * v.x11;
            C2 = A * v.x21;
            C3 = A * v.x31;
        }
        {
            Vector3<TnBase&> A( x12, x22, x32 );
            C1 = C1 + A * v.x12;
            C2 = C2 + A * v.x22;
            C3 = C3 + A * v.x32;
        }
        {
            Vector3<TnBase&> A( x13, x23, x33 );
            C1 = C1 + A * v.x13;
            C2 = C2 + A * v.x23;
            C3 = C3 + A * v.x33;
        }
        x11 = C1.x;
        x12 = C2.x;
        x13 = C3.x;
        x21 = C1.y;
        x22 = C2.y;
        x23 = C3.y;
        x31 = C1.z;
        x32 = C2.z;
        x33 = C3.z;
    }

    // Vector3<TnBase&> somehow doesn't compile here
    template<class U> inline
    void In_Place_Transpose_Times( const Matrix3<U>& v ) {
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "In Place Times Transpose (Matrix3) must operate on same base type." );

        Vector3<TnBase> C1;
        Vector3<TnBase> C2;
        Vector3<TnBase> C3;
        {
            Vector3<TnBase> A( x11, x12, x13 );
            C1 = A * v.x11;
            C2 = A * v.x12;
            C3 = A * v.x13;
        }
        {
            Vector3<TnBase> A( x21, x22, x23 );
            C1 = C1 + A * v.x21;
            C2 = C2 + A * v.x22;
            C3 = C3 + A * v.x23;
        }
        {
            Vector3<TnBase> A( x31, x32, x33 );
            C1 = C1 + A * v.x31;
            C2 = C2 + A * v.x32;
            C3 = C3 + A * v.x33;
        }
        x11 = C1.x;
        x12 = C2.x;
        x13 = C3.x;
        x21 = C1.y;
        x22 = C2.y;
        x23 = C3.y;
        x31 = C1.z;
        x32 = C2.z;
        x33 = C3.z;
    }

    template<class U>
    ConstNonref_Matrix3 Times_Transpose( const Matrix3<U>& v ) {
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "Times Transpose (Matrix3) must operate on same base type." );
        return ConstNonref_Matrix3(x11*v.x11 + x12*v.x12 + x13*v.x13,
                            x21*v.x11 + x22*v.x12 + x23*v.x13,
                            x31*v.x11 + x32*v.x12 + x33*v.x13,
                            x11*v.x21 + x12*v.x22 + x13*v.x23,
                            x21*v.x21 + x22*v.x22 + x23*v.x23,
                            x31*v.x21 + x32*v.x22 + x33*v.x23,
                            x11*v.x31 + x12*v.x32 + x13*v.x33,
                            x21*v.x31 + x22*v.x32 + x23*v.x33,
                            x31*v.x31 + x32*v.x32 + x33*v.x33);
    }

    template<class U>
    ConstNonref_Matrix3 Transpose_Times( const Matrix3<U>& v ) {
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "Times Transpose (Matrix3) must operate on same base type." );
        return ConstNonref_Matrix3(x11*v.x11 + x21*v.x21 + x31*v.x31,
                                   x12*v.x11 + x22*v.x21 + x32*v.x31,
                                   x13*v.x11 + x32*v.x21 + x33*v.x31,
                                   x11*v.x12 + x21*v.x22 + x31*v.x32,
                                   x12*v.x12 + x22*v.x22 + x32*v.x32,
                                   x13*v.x12 + x23*v.x22 + x33*v.x32,
                                   x11*v.x13 + x21*v.x23 + x31*v.x33,
                                   x12*v.x13 + x22*v.x23 + x32*v.x33,
                                   x13*v.x13 + x23*v.x23 + x33*v.x33);
    }

    // Compute the symmetric matrix A*D*At (A=*this, D=v)
    template<class U>
    ConstNonref_SymMatrix3 Conjugate( const DiagonalMatrix3<U>& v ) {
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "Conjugate (Matrix3) must operate on same base type." );
        return ConstNonref_SymMatrix3(v.x11*x11*x11+v.x22*x12*x12+v.x33*x13*x13,
                                      v.x11*x21*x11+v.x22*x22*x12+v.x33*x23*x13,
                                      v.x11*x31*x11+v.x22*x32*x12+v.x33*x33*x13,
                                      v.x11*x21*x21+v.x22*x22*x22+v.x33*x23*x23,
                                      v.x11*x31*x21+v.x22*x32*x22+v.x33*x33*x23,
                                      v.x11*x31*x31+v.x22*x32*x32+v.x33*x33*x33);
    }

    // Compute the symmetric matrix At*D*A (A=*this, D=v)
    template<class U>
    ConstNonref_SymMatrix3 Transpose_Conjugate( const DiagonalMatrix3<U>& v ) {
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "Conjugate (Matrix3) must operate on same base type." );
        return ConstNonref_SymMatrix3(v.x11*x11*x11+v.x22*x21*x21+v.x33*x31*x31,
                                      v.x11*x12*x11+v.x22*x22*x21+v.x33*x32*x31,
                                      v.x11*x13*x11+v.x22*x23*x21+v.x33*x33*x31,
                                      v.x11*x12*x12+v.x22*x22*x22+v.x33*x32*x32,
                                      v.x11*x13*x12+v.x22*x23*x22+v.x33*x33*x32,
                                      v.x11*x13*x13+v.x22*x23*x23+v.x33*x33*x33);
    }

    template<class U>
    ConstNonref_Vector3 operator*( const Vector3<U>& v ) {
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "Vector multiply (Matrix3) must operate on same base type." );
        return ConstNonref_Vector3(x11*v.x + x12*v.y + x13*v.z,
                                   x21*v.x + x22*v.y + x23*v.z,
                                   x31*v.x + x32*v.y + x33*v.z);
    }


    template<class U>
    ConstNonref_Vector3 Times_Transpose( const Vector3<U>& v ) {
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "Vector multiply (Matrix3) must operate on same base type." );
        return ConstNonref_Vector3(x11*v.x + x21*v.y + x31*v.z,
                                   x12*v.x + x22*v.y + x32*v.z,
                                   x13*v.x + x23*v.y + x33*v.z);
    }


    ConstNonref_DiagMatrix3 Extract_Diagonal_Component(){
        return ConstNonref_DiagMatrix3( x11, x22, x33 );
    }

    ConstNonref_UniMatrix3 Extract_Lower_Unitriangular_Component(){
        return ConstNonref_UniMatrix3( x21, x31, x32 );
    }

    ConstNonref_UniMatrix3 Extract_Upper_Unitriangular_Component(){
        return ConstNonref_UniMatrix3( x12, x13, x23 );
    }


};


#endif
