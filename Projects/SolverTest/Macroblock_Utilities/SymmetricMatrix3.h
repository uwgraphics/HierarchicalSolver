#ifndef __KERNEL_SYMMETRICMATRIX3_H__
#define __KERNEL_SYMMETRICMATRIX3_H__

#include "stdlib.h"
#include <type_traits>

template<class Tn> struct DiagonalMatrix3;
template<class Tn> struct UnitriangularMatrix3;

template<class Tn> struct SymmetricMatrix3 {
    typedef Tn TnType;
    typedef typename std::decay<Tn>::type TnBase;
    typedef SymmetricMatrix3< typename std::add_const< typename std::remove_reference< Tn >::type >::type > ConstNonref_SymMatrix3;
    typedef UnitriangularMatrix3< typename std::add_const< typename std::remove_reference< Tn >::type >::type > ConstNonref_UniMatrix3;
    typedef DiagonalMatrix3< typename std::add_const< typename std::remove_reference< Tn >::type >::type > ConstNonref_DiagMatrix3;
    static const bool IS_REF =  std::is_reference< Tn >::value;
    static const bool IS_CONST = std::is_const< typename std::remove_reference< Tn >::type >::value;

    Tn x11, x21, x31, x22, x32, x33;

    SymmetricMatrix3() :
    x11(), x21(), x31(), x22(), x32(), x33() 
    {}

#define CTYPE SymmetricMatrix3
#define INITIALIZER_LIST x11(v.x11), x21(v.x21), x31(v.x31), \
                         x22(v.x22), x32(v.x32), x33(v.x33)
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


    explicit SymmetricMatrix3(const Tn& x11_in,const Tn& x21_in,const Tn& x31_in,
                              const Tn& x22_in,const Tn& x32_in,const Tn& x33_in)
        :x11(x11_in),x21(x21_in),x31(x31_in),x22(x22_in),x32(x32_in),x33(x33_in) {}
    
    template< class U>
    SymmetricMatrix3& operator=(const U& v)
    {
        static_assert( !IS_CONST,
                       "Implicit assignment operator cannot be used for const-types." );
        x11 = v.x11; 
        x21 = v.x21;
        x31 = v.x31;
        x22 = v.x22;
        x32 = v.x32;
        x33 = v.x33;
        return *this;
    }

    template< class U>
    SymmetricMatrix3& operator=(U& v)
    {
        static_assert( !IS_CONST,
                       "Implicit assignment operator cannot be used for const-types." );
        x11 = v.x11; 
        x21 = v.x21;
        x31 = v.x31;
        x22 = v.x22;
        x32 = v.x32;
        x33 = v.x33;
        return *this;
    }

    template<class T_DATA>
    void Load_Aligned(const T_DATA (&A)[6])
    {
        static_assert( !IS_CONST,
                       "Load Aligned cannot be used for const-types." );
        x11.Load_Aligned(A[0]);
        x21.Load_Aligned(A[1]);
        x31.Load_Aligned(A[2]);
        x22.Load_Aligned(A[3]);
        x32.Load_Aligned(A[4]);
        x33.Load_Aligned(A[5]);   
    }

    template<class T_DATA>
    void Load_Aligned(const T_DATA (&A)[9])
    {
        static_assert( !IS_CONST,
                       "Load Aligned cannot be used for const-types." );
        x11.Load_Aligned(A[0]);
        x21.Load_Aligned(A[1]);
        x31.Load_Aligned(A[2]);
        x22.Load_Aligned(A[4]);
        x32.Load_Aligned(A[5]);
        x33.Load_Aligned(A[8]);   
    }

    template<class T_DATA>
    void Store(T_DATA (&A)[6]) const
    {
        ::Store(A[0],x11);
        ::Store(A[1],x21);
        ::Store(A[2],x31);
        ::Store(A[3],x22);
        ::Store(A[4],x32);
        ::Store(A[5],x33);   
    }

    template<class T_DATA>
    void Store(T_DATA (&A)[9]) const
    {
        ::Store(A[0],x11);
        ::Store(A[1],x21);
        ::Store(A[2],x31);
        ::Store(A[4],x22);
        ::Store(A[5],x32);
        ::Store(A[8],x33);   
    }

    template<class T_DATA>
    void StoreFull(T_DATA (&A)[9]) const
    {
        ::Store(A[0],x11);
        ::Store(A[1],x21);
        ::Store(A[2],x31);
        ::Store(A[3],x21);
        ::Store(A[4],x22);
        ::Store(A[5],x32);
        ::Store(A[6],x31);
        ::Store(A[7],x32);
        ::Store(A[8],x33);   
    }

    // Accessors

    Tn& operator()(int i, int j){
        if( i<1 || i > 3 || j<1 || j>3)
            exit(1);
        if( i==1 && j==1 )
            return x11;
        if( i==1 && j==2 )
            return x21;
        if( i==1 && j==3 )
            return x31;
        if( i==2 && j==1 )
            return x21;
        if( i==2 && j==2 )
            return x22;
        if( i==2 && j==3 )
            return x32;
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
            return x21;
        if( i==1 && j==3 )
            return x31;
        if( i==2 && j==1 )
            return x21;
        if( i==2 && j==2 )
            return x22;
        if( i==2 && j==3 )
            return x32;
        if( i==3 && j==1 )
            return x31;
        if( i==3 && j==2 )
            return x32;
        if( i==3 && j==3 )
            return x33;    
    }

    // Scalar Math

    SymmetricMatrix3& operator+=( const Tn& n ){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Scalar plus-equals operator (SymmetricMatrix3) cannot be used for const-types." );
        x11 = x11+n; x21 = x21+n; x31 = x31+n; x22 = x22+n; x32 = x32+n; x33 = x33+n
        return *this;
    }

    SymmetricMatrix3& operator*=( const Tn& n ){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Scalar times-equals operator (SymmetricMatrix3) cannot be used for const-types." );
        x11 = x11*n; x21 = x21*n; x31 = x31*n; x22 = x22*n; x32 = x32*n; x33 = x33*n;
        return *this;
    }

    SymmetricMatrix3& operator-=( const Tn& n ){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Scalar minus-equals operator (SymmetricMatrix3) cannot be used for const-types." );
        x11 = x11-n; x21 = x21-n; x31 = x31-n; x22 = x22-n; x32 = x32-n; x33 = x33-n;
        return *this;
    }

    SymmetricMatrix3& operator/=( const Tn& n ){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Scalar divide-equals operator (SymmetricMatrix3) cannot be used for const-types." );
        x11 = x11/n; x21 = x21/n; x31 = x31/n; x22 = x22/n; x32 = x32/n; x33 = x33/n;
        return *this;
    }
    

    // Matrix Math

    template<class U>
    SymmetricMatrix3& operator+=( const SymmetricMatrix3<U>& v ) {
        static_assert( std::is_same<Tn, typename std::decay<U>::type >::value,
                       "Matrix plus-equals operator (SymmetricMatrix3) must operate on same base type." );
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Matrix plus-equals operator (SymmetricMatrix3) cannot be used for const-types." );
        x11 = x11+v.x11; x21 = x21+v.x21; x31 = x31+v.x31; x22 = x22+v.x22; x32 = x32+v.x32; x33 = x33+v.x33;
        return *this;
    }

    template<class U>
    SymmetricMatrix3& operator-=( const SymmetricMatrix3<U>& v ) {
        static_assert( std::is_same<Tn, typename std::decay<U>::type >::value,
                       "Matrix minus-equals operator (SymmetricMatrix3) must operate on same base type." );
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "Matrix minus-equals operator (SymmetricMatrix3) cannot be used for const-types." );
        x11 = x11-v.x11; x21 = x21-v.x21; x31 = x31-v.x31; x22 = x22-v.x22; x32 = x32-v.x32; x33 = x33-v.x33;
        return *this;
    }


    // Factorization
    void Decomposition_LDL(){
        static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                       "LDLt Factorization (SymmetricMatrix3) cannot be used for const-types." );
#if 0
        SymmetricMatrix3& A = (*this);

        for(int j=1;j<=3;j++){
            A(j,j)=A(j,j).inverse();
            for(int l=j+1;l<=3;l++){
                TnBase temp=A(l,j)*A(j,j);
                for(int k=l;k<=3;k++)
                    A(k,l)=A(k,l)-A(k,j)*temp;
                A(l,j)=temp;}}
#else
        TnBase temp;
        x11=x11.inverse(); // i=1
        temp=x21*x11;      //     l=2
        x22=x22-x21*temp;  //         k=2
        x32=x32-x31*temp;  //         k=3
        x21=temp;
        temp=x31*x11;      //     l=3
        x33=x33-x31*temp;  //         k=3
        x31=temp;
        x22=x22.inverse(); // i=2
        temp=x32*x22;      //     l=3
        x33=x33-x32*temp;  //         k=3
        x32=temp;
        x33=x33.inverse(); // i=3
#endif

    }

    template<class U>
    void LDLt(UnitriangularMatrix3<U>& L,DiagonalMatrix3<U>& D) const {
        static_assert( std::is_same<Tn, typename std::decay<U>::type >::value,
                       "LDLt (SymmetricMatrix3) must operate on same base type." );
        D.x11=x11;
        L.x21=x21;
        L.x31=x31;
        D.x22=x22;
        L.x32=x32;
        D.x33=x33;

        TnBase temp;
        D.x11=D.x11.inverse(); // i=1
        temp=L.x21*D.x11;      //     l=2
        D.x22=D.x22-L.x21*temp;  //         k=2
        L.x32=L.x32-L.x31*temp;  //         k=3
        L.x21=temp;
        temp=L.x31*D.x11;      //     l=3
        D.x33=D.x33-L.x31*temp;  //         k=3
        L.x31=temp;
        D.x22=D.x22.inverse(); // i=2
        temp=L.x32*D.x22;      //     l=3
        D.x33=D.x33-L.x32*temp;  //         k=3
        L.x32=temp;
        D.x33=D.x33.inverse(); // i=3

    }

    ConstNonref_DiagMatrix3 Extract_Diagonal_Component(){
        return ConstNonref_DiagMatrix3( x11, x22, x33 );
    }

    ConstNonref_UniMatrix3 Extract_Unitriangular_Component(){
        return ConstNonref_UniMatrix3( x21, x31, x32 );
    }

};


#endif
