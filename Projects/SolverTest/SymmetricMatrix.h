//#####################################################################
// Copyright 2018, Eftychios Sifakis, Qisi Wang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef  __SymmetricMatrix_h__
#define  __SymmetricMatrix_h__

#include <cassert>
#include <iostream>
#include <iomanip>

//################################################
// Class SymmetricMatrix
//################################################
template<class T>
struct SymmetricMatrix
// Symmetric matrix in packed form in column major
{
    size_t n;   // Matrix dimension
    T *data; // Pointer to lower triangular data

    inline const T& operator()(size_t i, size_t j) const {
        assert((0<=j) && (j<=i) && (i<=n));
        return data[j*n+i-(j*(j+1))/2];
    }

    inline T& operator()(size_t i, size_t j) {
        assert((0<=j) && (j<=i) && (i<=n));
        return data[j*n+i-(j*(j+1))/2];
    }

    void Conjugate() {
        // assume the diagonal terms are inverted
        for (int j=n-1; j>=0; j--)
        {
            (*this)(j,j) = 1/(*this)(j,j);
            for (int l=j+1; l<n; l++)
                for (int k=l; k<n; k++)
                    (*this)(k,l) += (*this)(k,j)*(*this)(l,j)*(*this)(j,j);
            for (int i=j+1; i<n; i++)
                (*this)(i,j) *= (*this)(j,j);
            //(*this)(j,j) *= (*this)(j,j);

        }
    }

    void Cholesky(int k) {
        // assume the diagonal terms are inverted
        for (int j=0; j<k; j++)
        {
             (*this)(j,j) = 1/(*this)(j,j);
             for (int i=j+1; i<n; i++){
                 (*this)(i,i) -= (*this)(i,j)*(*this)(i,j)*(*this)(j,j); // A_ii -= A_ij*A_ij*(1/D_jj)
                 (*this)(i,j) *= (*this)(j,j);
                 for (int l=i+1; l<n; l++)
                     (*this)(l,i) -= (*this)(l,j)*(*this)(i,j);          // A_li -= A_lj*L_ij
             }
        }
    }

    void print(int offset = 4) const
        {
            for (int i=0; i<n; i++){
                for (int j=0; j<=i; j++)
                    std::cout<<std::setw(offset)<<(*this)(i,j);
                std::cout<<std::endl;
            }
        }
};
//#####################################################################
#endif
