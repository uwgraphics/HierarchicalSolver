//#####################################################################
// Copyright 2002-2009, Eftychios Sifakis, Michael Doescher, Qisi Wang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_CSC_MATRIX
//#####################################################################
#ifndef __ARRAY_CSC_MATRIX__
#define __ARRAY_CSC_MATRIX__


#include <utility>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include "MATRIX_CONVERSIONS.h"

#include <fstream>

namespace PhysBAM{

template <class T>
class ARRAY_CSC_MATRIX {

public:
    ARRAY<ARRAY<int> > row_data;
    ARRAY<ARRAY<T> > data;

    // constructor
ARRAY_CSC_MATRIX(int n) :
    row_data(n), data(n)
    {
    }

    void Resize(int n)
    {
        row_data.Resize(n);
        data.Resize(n);
    }

    // count nonzero entries
    int nnz() {
        const int n = row_data.Size();
        int A_nnz = 0;
        for(int j=1;j<=row_data.Size();j++)
                A_nnz+=row_data(j).Size();
        return A_nnz;
    }

    int nnz(std::pair<int, int> range)
    {

        const int start_index = range.first;
        const int end_index = range.second;
        //LOG::cout<<start_index<<", "<<end_index<<std::endl;
        PHYSBAM_ASSERT(start_index>=1);
        PHYSBAM_ASSERT(end_index<=row_data.Size());
        int count = 0;
        for (int j=start_index; j<=end_index; j++)
            for (int i=1; i<=row_data(j).Size(); i++)
                if (row_data(j)(i)<=end_index)
                    count++;
        return count;
    }



    // print the matrix pattern
    void print_nz_patterns() {
        const int n = row_data.Size();
        ARRAY<std::string> nz_patterns(n);
        for(int j=1;j<=n;j++){
            int ii=1;
            for(int i=j;i<=n;i++){
                if(ii<=row_data(j).Size() && row_data(j)(ii)==i) { nz_patterns(i)+="* "; ii++; }
                else nz_patterns(i)+="  ";}
            LOG::cout<<nz_patterns(j)<<std::endl;}
    }

    void Multiply(const T* x, T* b) const
       {
        const int n = row_data.Size();
        for (int j=0;j<n;j++)
        {
            b[j] = (T)0;
        }

        for (int j=1; j<=n; j++){
            b[j-1] += data(j)(1)*x[j-1];
            for (int ii=2; ii<=row_data(j).Size();ii++)
            {
                const int i = row_data(j)(ii);
                b[i-1] += data(j)(ii)*x[j-1];
                b[j-1] += data(j)(ii)*x[i-1];
            }
        }
    }

    void Submatrix(std::pair<int, int> range, ARRAY_CSC_MATRIX<T>& CSC)
    {
        ARRAY<ARRAY<int> >& sub_Mat=CSC.row_data;

        const int start_index = range.first;
        const int end_index = range.second;

        PHYSBAM_ASSERT(sub_Mat.Size() == (end_index-start_index+1));
        for (int j = start_index; j<=end_index; j++)
            for (int i=1; i<=row_data(j).Size(); i++)
                if(row_data(j)(i)<=end_index)
                    sub_Mat(j-start_index+1).Append(row_data(j)(i)-start_index+1);
    }

    // symbolic_cholesky decomposition
    void symbolic_cholesky(ARRAY_CSC_MATRIX<T>& L_data) {
        ARRAY<ARRAY<int> >& L_CSC=L_data.row_data;
        PHYSBAM_ASSERT(row_data.Size() == L_CSC.Size());
        const int n = row_data.Size();
        ARRAY<int> parent(n);
        for(int i=1;i<=n;i++){
            HASHTABLE<int> hashtable;
            for(int ii=1;ii<=row_data(i).Size();ii++)
                hashtable.Set(row_data(i)(ii));
            for(int j=1;j<i;j++)
                if(parent(j)==i)
                    for(int ii=1;ii<=L_CSC(j).Size();ii++)
                        if(L_CSC(j)(ii)!=j) hashtable.Set(L_CSC(j)(ii));
            for(HASHTABLE_ITERATOR<int> iterator(hashtable);iterator.Valid();iterator.Next())
                L_CSC(i).Append(iterator.Key());
            Sort(L_CSC(i));
            if(L_CSC(i).Size()>1) parent(i)=L_CSC(i)(2);}
    }

    // fill with random symmetric positive definite entries
    void Generate_Random_Positive_Definite_Matrix(const int seed=1) {
        RANDOM_NUMBERS<T> random;
        random.Set_Seed(seed);
        const int n=row_data.m;
        data.Resize(n);
        for (int i=1;i<=n;i++) {
            data(i).Resize(row_data(i).m);
            for (int j=1;j<=row_data(i).m;j++) {
                data(i)(j)=random.Get_Number();}
            data(i)(1)+=(n+1);
        }

    }

    //      ARRAY_CSC_MATRIX_To_Hashtable
    void To_Hashtable(HASHTABLE<VECTOR<int,2>,T>& hashtable_matrix) {
        for (int i=1;i<=row_data.m;i++) {
            for (int j=1;j<=row_data(i).m;j++) {
                hashtable_matrix.Set(VECTOR<int,2>(i,row_data(i)(j)),data(i)(j));
                hashtable_matrix.Set(VECTOR<int,2>(row_data(i)(j),i),data(i)(j));}}
    }

    // PARDISO_Export
    void PARDISO_Export() {
        MATRIX_CONVERSIONS<T> mc;
        HASHTABLE<VECTOR<int,2>,T> hashtable_matrix;
        To_Hashtable(hashtable_matrix);
        mc.Hashtable_To_PARDISO(hashtable_matrix);
    }

    // void Subdomain_Matricies(ARRAY<HASHTABLE<VECTOR<int,2>,T> >& data) {}

    // void Extract_Subdomain_Matrix(const int subdomain_id, Type& data
    void Dump_Sparsity(std::string base = "z")
    {
        std::ofstream iout(base+"_i.out",std::ios::binary|std::ios::out|std::ios::trunc);
        std::ofstream jout(base+"_j.out",std::ios::binary|std::ios::out|std::ios::trunc);
        LOG::cout<<"opened file"<<std::endl;
        //ARRAY<ARRAY<int> > row_data;
        for (int j=1; j<=row_data.Size(); j++)
            for (int i=1; i<=row_data(j).Size(); i++)
            {
                iout.write((char*)&row_data(j)(i), sizeof(int));
                jout.write((char*)&j, sizeof(int));
            }
        iout.close();
        jout.close();

    }


};

}
#endif
