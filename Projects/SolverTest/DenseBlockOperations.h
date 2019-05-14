//#####################################################################
// Copyright 2018, Eftychios Sifakis, Qisi Wang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef  __DenseBlockOperations_h__
#define  __DenseBlockOperations_h__

#include <cassert>
#include "SymmetricMatrix.h"

template<class Arch_Type>
struct Arch_Policy;

template<>
struct Arch_Policy<double> {
    typedef double Data_Type;
};

template<>
struct Arch_Policy<float> {
    typedef float Data_Type;
};


//################################################
// Class DenseBlockOperations <double>
//################################################
template<int Block_Size, class Arch_Type>
    struct DenseBlockOperations;

template<int Block_Size>
struct DenseBlockOperations<Block_Size, double> {
    using T = typename  Arch_Policy<double>::Data_Type;

    static void Cholesky(int n, T* data, int k){
        assert(k<=n);
        SymmetricMatrix<double> mat{n, data};
        for (int j=0; j<k; j+=Block_Size)
        {
            T* head_ptr = &mat(j,j);
            int size = n-j;
            int remaining = k-j>Block_Size?Block_Size:k-j;
            LDLT_Decomposition(size, head_ptr, remaining);
            for (int i=Block_Size; i<size; i+=Block_Size)
                Forwards_Substitution_On_Transpose(size, head_ptr, remaining, i);
            for (int i=Block_Size; i<size; i+=Block_Size) {
                Subtract_Conjugate_By_Diagonal(size, head_ptr, remaining, i);
                Scale_Column_By_Diagonal(size, head_ptr, remaining, i);
                for (int l=i+Block_Size; l<size; l+=Block_Size)
                    Subtract_Times_Transpose(size, head_ptr, remaining, l, i);
            }
        }
    }

    static void LDLT_Decomposition(int n, T* data, int k) {
        assert(k<=Block_Size && k<=n);
        SymmetricMatrix<double> mat{n, data};

        for (int j=0; j<k; j++) {
            mat(j,j) = 1/mat(j,j);
            for (int i=j+1; i<Block_Size&&i<n; i++){
                mat(i,i) -= mat(i,j)*mat(i,j)*mat(j,j);     // A_ii -= A_ij*A_ij*(1/D_jj)
                mat(i,j) *= mat(j,j);
                for (int l=i+1; l<Block_Size&&l<n; l++)
                    mat(l,i) -= mat(l,j)*mat(i,j);          // A_li -= A_lj*L_ij
            }
        }
    }

    static void Forwards_Substitution_On_Transpose(int n, T* data, int k, int p) {
        // p is the starting row index
        assert(k<=Block_Size && k<=n && p>=Block_Size && p<n);
        SymmetricMatrix<double> mat{n, data};

        for (int j=0; j<k; j++)
            for (int i=j+1; i<Block_Size && i<n; i++)
                for (int l=p; l<p+Block_Size&&l<n; l++)
                    mat(l,i) -= mat(l,j)*mat(i,j);
    }

    static void Scale_Column_By_Diagonal(int n, T* data, int k, int p) {
        // m is the starting row index
        assert(k<=Block_Size && k<=n && p>=Block_Size && p<n);
        SymmetricMatrix<double> mat{n, data};

        for (int j=0; j<k; j++)
            for (int i=p; i<p+Block_Size&&i<n; i++)
                mat(i,j) *= mat(j,j);
    }

    static void Subtract_Times_Transpose(int n, T* data, int k, int p, int q) {
        assert(k<=Block_Size && k<=n && p>=Block_Size && p<n && q>=Block_Size && q<n && p>=q+Block_Size);
        SymmetricMatrix<double> mat{n, data};
        for (int l=0; l<k; l++)
            for (int j=q; j<q+Block_Size; j++)
                for (int i=p; i<p+Block_Size&&i<n; i++)
                    mat(i,j) -= mat(i,l)*mat(j,l);
    }

    static void Subtract_Conjugate_By_Diagonal(int n, T* data, int k, int p) {
        assert(k<=Block_Size && k<=n && p>=Block_Size && p<n);
        SymmetricMatrix<double> mat{n, data};
        for (int l=0; l<k; l++)
            for (int j=p; j<p+Block_Size&&j<n; j++)
                for (int i=j; i<p+Block_Size&&i<n; i++)
                    mat(i,j) -= mat(i,l)*mat(j,l)*mat(l,l);
    }
};

//#####################################################################

//################################################
// Class DenseBlockOperations <float>
//################################################
template<int Block_Size, class Arch_Type>
    struct DenseBlockOperations;

template<int Block_Size>
struct DenseBlockOperations<Block_Size, float> {
    using T = typename  Arch_Policy<float>::Data_Type;

    static void Cholesky(int n, T* data, int k){
        assert(k<=n);
        SymmetricMatrix<float> mat{n, data};
        for (int j=0; j<k; j+=Block_Size)
        {
            T* head_ptr = &mat(j,j);
            int size = n-j;
            int remaining = k-j>Block_Size?Block_Size:k-j;
            LDLT_Decomposition(size, head_ptr, remaining);
            for (int i=Block_Size; i<size; i+=Block_Size)
                Forwards_Substitution_On_Transpose(size, head_ptr, remaining, i);
            for (int i=Block_Size; i<size; i+=Block_Size) {
                Subtract_Conjugate_By_Diagonal(size, head_ptr, remaining, i);
                Scale_Column_By_Diagonal(size, head_ptr, remaining, i);
                for (int l=i+Block_Size; l<size; l+=Block_Size)
                    Subtract_Times_Transpose(size, head_ptr, remaining, l, i);
            }
        }
    }

    static void LDLT_Decomposition(int n, T* data, int k) {
        assert(k<=Block_Size && k<=n);
        SymmetricMatrix<float> mat{n, data};

        for (int j=0; j<k; j++) {
            mat(j,j) = 1/mat(j,j);
            for (int i=j+1; i<Block_Size&&i<n; i++){
                mat(i,i) -= mat(i,j)*mat(i,j)*mat(j,j);     // A_ii -= A_ij*A_ij*(1/D_jj)
                mat(i,j) *= mat(j,j);
                for (int l=i+1; l<Block_Size&&l<n; l++)
                    mat(l,i) -= mat(l,j)*mat(i,j);          // A_li -= A_lj*L_ij
            }
        }
    }

    static void Forwards_Substitution_On_Transpose(int n, T* data, int k, int p) {
        // p is the starting row index
        assert(k<=Block_Size && k<=n && p>=Block_Size && p<n);
        SymmetricMatrix<float> mat{n, data};

        for (int j=0; j<k; j++)
            for (int i=j+1; i<Block_Size && i<n; i++)
                for (int l=p; l<p+Block_Size&&l<n; l++)
                    mat(l,i) -= mat(l,j)*mat(i,j);
    }

    static void Scale_Column_By_Diagonal(int n, T* data, int k, int p) {
        // m is the starting row index
        assert(k<=Block_Size && k<=n && p>=Block_Size && p<n);
        SymmetricMatrix<float> mat{n, data};

        for (int j=0; j<k; j++)
            for (int i=p; i<p+Block_Size&&i<n; i++)
                mat(i,j) *= mat(j,j);
    }

    static void Subtract_Times_Transpose(int n, T* data, int k, int p, int q) {
        assert(k<=Block_Size && k<=n && p>=Block_Size && p<n && q>=Block_Size && q<n && p>=q+Block_Size);
        SymmetricMatrix<float> mat{n, data};
        for (int l=0; l<k; l++)
            for (int j=q; j<q+Block_Size; j++)
                for (int i=p; i<p+Block_Size&&i<n; i++)
                    mat(i,j) -= mat(i,l)*mat(j,l);
    }

    static void Subtract_Conjugate_By_Diagonal(int n, T* data, int k, int p) {
        assert(k<=Block_Size && k<=n && p>=Block_Size && p<n);
        SymmetricMatrix<float> mat{n, data};
        for (int l=0; l<k; l++)
            for (int j=p; j<p+Block_Size&&j<n; j++)
                for (int i=j; i<p+Block_Size&&i<n; i++)
                    mat(i,j) -= mat(i,l)*mat(j,l)*mat(l,l);
    }
};

//#####################################################################
#endif //__DenseBlockOperations_h__
