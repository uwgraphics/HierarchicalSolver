#ifndef __KERNEL_UNITRIANGULARMATRIX3_H__
#define __KERNEL_UNITRIANGULARMATRIX3_H__

#include "stdlib.h"
#include <type_traits>


template<class Tn> struct Vector3;
template<class Tn> struct Matrix3;

template<class Tn> struct UnitriangularMatrix3 {
    typedef Tn TnType;
    typedef typename std::decay<Tn>::type TnBase;
    typedef Matrix3< typename std::add_const< typename std::remove_reference< Tn >::type >::type > ConstNonref_Matrix3;
    typedef Vector3< typename std::add_const< typename std::remove_reference< Tn >::type >::type > ConstNonref_Vector3;
    typedef UnitriangularMatrix3< typename std::add_const< typename std::remove_reference< Tn >::type >::type > ConstNonref_UniMatrix3;
    static const bool IS_REF =  std::is_reference< Tn >::value;
    static const bool IS_CONST = std::is_const< typename std::remove_reference< Tn >::type >::value;

    Tn x21, x31, x32;

    UnitriangularMatrix3():
    x21(), x31(), x32() 
    {}

#define CTYPE UnitriangularMatrix3
#define INITIALIZER_LIST x21(v.x21), x31(v.x31), x32(v.x32)
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

    explicit UnitriangularMatrix3(const Tn& x21_in,const Tn& x31_in,const Tn& x32_in)
        :x21(x21_in),x31(x31_in),x32(x32_in) {}
    
    template< class U>
    UnitriangularMatrix3& operator=(const U& v)
    {
        static_assert( !IS_CONST,
                       "Implicit assignment operator cannot be used for const-types." );
        x21 = v.x21;
        x31 = v.x31;
        x32 = v.x32;
        return *this;
    }

    template< class U>
    UnitriangularMatrix3& operator=(U& v)
    {
        static_assert( !IS_CONST,
                       "Implicit assignment operator cannot be used for const-types." );
        x21 = v.x21;
        x31 = v.x31;
        x32 = v.x32;
        return *this;
    }

    template<class T_DATA>
    void Load_Aligned(const T_DATA (&A)[3])
    {
        static_assert( !IS_CONST,
                       "Load Aligned cannot be used for const-types." );
        x21.Load_Aligned(A[0]);
        x31.Load_Aligned(A[1]);
        x32.Load_Aligned(A[2]);
    }

    template<class T_DATA>
    void Load_Aligned(const T_DATA (&A)[6])
    {
        static_assert( !IS_CONST,
                       "Load Aligned cannot be used for const-types." );
        x21.Load_Aligned(A[1]);
        x31.Load_Aligned(A[2]);
        x32.Load_Aligned(A[4]);
    }

    template<class T_DATA>
    void Load_Aligned(const T_DATA (&A)[9])
    {
        static_assert( !IS_CONST,
                       "Load Aligned cannot be used for const-types." );
        x21.Load_Aligned(A[1]);
        x31.Load_Aligned(A[2]);
        x32.Load_Aligned(A[5]);
    }

    template<class T_DATA>
    void Store(T_DATA (&A)[3]) const
    {
        ::Store(A[0],x21);
        ::Store(A[1],x31);
        ::Store(A[2],x32);   
    }

    template<class T_DATA>
    void Store(T_DATA (&A)[6]) const
    {
        ::Store(A[1],x21);
        ::Store(A[2],x31);
        ::Store(A[4],x32);   
    }

    template<class T_DATA>
    void Store(T_DATA (&A)[9]) const
    {
        ::Store(A[1],x21);
        ::Store(A[2],x31);
        ::Store(A[5],x32);   
    }


    // Accessors

    Tn& operator()(int i, int j){
        if( i<1 || i > 3 || j<1 || j>3 || i==j)
            exit(1);
        if( i==1 && j==2 )
            return x21;
        if( i==1 && j==3 )
            return x31;
        if( i==2 && j==1 )
            return x21;
        if( i==2 && j==3 )
            return x32;
        if( i==3 && j==1 )
            return x31;
        if( i==3 && j==2 )
            return x32;
    }

    const Tn& operator()(int i, int j) const {
        if( i<1 || i > 3 || j<1 || j>3 || i==j)
            exit(1);
        if( i==1 && j==2 )
            return x21;
        if( i==1 && j==3 )
            return x31;
        if( i==2 && j==1 )
            return x21;
        if( i==2 && j==3 )
            return x32;
        if( i==3 && j==1 )
            return x31;
        if( i==3 && j==2 )
            return x32;
    }

    // Methods
    
    
    // Scalar Math

    UnitriangularMatrix3& operator+=( const Tn& n ){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Scalar plus-equals operator (UnitriangularMatrix3) cannot be used for const-types." );
        x21 = x21+n; x31 = x31+n; x32 = x32+n;
        return *this;
    }

    UnitriangularMatrix3& operator*=( const Tn& n ){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Scalar times-equals operator (UnitriangularMatrix3) cannot be used for const-types." );
        x21 = x21*n; x31 = x31*n; x32 = x32*n;
        return *this;
    }

    UnitriangularMatrix3& operator-=( const Tn& n ){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Scalar minus-equals operator (UnitriangularMatrix3) cannot be used for const-types." );
        x21 = x21-n; x31 = x31-n; x32 = x32-n;
        return *this;
    }

    UnitriangularMatrix3& operator/=( const Tn& n ){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Scalar divide-equals operator (UnitriangularMatrix3) cannot be used for const-types." );
        x21 = x21/n; x31 = x31/n; x32 = x32/n;
        return *this;
    }
    

    // Matrix Math

    template<class U>
    UnitriangularMatrix3& operator+=( const UnitriangularMatrix3<U>& v ) {
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Matrix plus-equals operator (UnitriangularMatrix3) cannot be used for const-types." );
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "Matrix plus-equals operator (UnitriangularMatrix3) must operate on same base type." );
        x21 = x21+v.x21; x31 = x31+v.x31; x32 = x32+v.x32; 
        return *this;
    }

    template<class U>
    UnitriangularMatrix3& operator-=( const UnitriangularMatrix3<U>& v ) {
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Matrix minus-equals operator (UnitriangularMatrix3) cannot be used for const-types." );
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "Matrix minus-equals operator (UnitriangularMatrix3) must operate on same base type." );
        x21 = x21-v.x21; x31 = x31-v.x31; x32 = x32-v.x32; 
        return *this;
    }

    template<class U>
    void In_Place_Forward_Substitution( Vector3<U>& v ) {
        static_assert( !std::is_const<std::remove_reference<U> >::value,
                       "In_Place_Forward_Substitution (UnitriangularMatrix3) cannot be used for const arg types." );
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "In_Place_Forward_Substitution (UnitriangularMatrix3) must operate on same base type." );
        v.x = v.x;
        v.y = (v.y - (x21*v.x));
        v.z = (v.z - ((x32*v.y)+(x31*v.x)));        
    }

    template<class U>
    void In_Place_Backward_Substitution( Vector3<U>& v ) {
        static_assert( !std::is_const<std::remove_reference<U> >::value,
                       "In_Place_Backward_Substitution (UnitriangularMatrix3) cannot be used for const arg types." );
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "In_Place_Backward_Substitution (UnitriangularMatrix3) must operate on same base type." );
        v.z = v.z;
        v.y = (v.y - (x32*v.z));
        v.x = (v.x - ((x21*v.y)+(x31*v.z)));        
    }

    template<class U>
    void In_Place_Forward_Substitution( Matrix3<U>& m ) {
        static_assert( !std::is_const<std::remove_reference<U> >::value,
                       "In_Place_Forward_Substitution (UnitriangularMatrix3) cannot be used for const arg types." );
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "In_Place_Forward_Substitution (UnitriangularMatrix3) must operate on same base type." );
        {   Vector3<TnBase&> column( m.x11, m.x21, m.x31 );
            In_Place_Forward_Substitution( column );  }
        {   Vector3<TnBase&> column( m.x12, m.x22, m.x32 );
            In_Place_Forward_Substitution( column );  }
        {   Vector3<TnBase&> column( m.x13, m.x23, m.x33 );
            In_Place_Forward_Substitution( column );  }
    }

    template<class U>
    void In_Place_Forward_Substitution_On_Transpose( Matrix3<U>& m ) {
        static_assert( !std::is_const<std::remove_reference<U> >::value,
                       "In_Place_Forward_Substitution_On_Transpose (UnitriangularMatrix3) cannot be used for const arg types." );
        static_assert( std::is_same<TnBase, typename std::decay<U>::type >::value,
                       "In_Place_Forward_Substitution_On_Transpose (UnitriangularMatrix3) must operate on same base type." );
        {   Vector3<TnBase&> row( m.x11, m.x12, m.x13 );
            In_Place_Forward_Substitution( row );  }
        {   Vector3<TnBase&> row( m.x21, m.x22, m.x23 );
            In_Place_Forward_Substitution( row );  }
        {   Vector3<TnBase&> row( m.x31, m.x32, m.x33 );
            In_Place_Forward_Substitution( row );  }
    }

};


#endif
