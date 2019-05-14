//#####################################################################
// Copyright 2009, Michael Doescher, Qisi Wang, Eftychios Sifakis
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_CONVERSIONS
//#####################################################################
#ifndef __MATRIX_CONVERSIONS__
#define __MATRIX_CONVERSIONS__

#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_HASHTABLE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>

#include "./ARRAY_CSC_MATRIX.h"

#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include <fstream>

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl.h"
#include "omp.h"

namespace PhysBAM{

template<class T>
bool Sort_Function(PAIR<VECTOR<int,2>,T> a, PAIR<VECTOR<int,2>,T> b) {
        if (a.x.x<b.x.x) return true;
        // for sorting the array version of the matrix
        if (a.x.x==b.x.x && a.x.y<b.x.y) return true;
        return false;
    }

template<class T>
struct MATRIX_CONVERSIONS
{
public:

    // #################################################################################################3
    //      Read_Hashtable_Matrix
    // #################################################################################################3
    void Read_Hashtable_Matrix(STREAM_TYPE& stream_type, HASHTABLE<VECTOR<int,2>,T>& hashtable_matrix, const std::string& filename) {
        FILE_UTILITIES::Read_From_File(stream_type,filename,hashtable_matrix);
    }

    // #################################################################################################3
    //      Write_Hashtable_Matrix
    // #################################################################################################3
    void Write_Hashtable_Matrix(STREAM_TYPE& stream_type, HASHTABLE<VECTOR<int,2>,T>& hashtable_matrix, const std::string& filename) {
        FILE_UTILITIES::Write_To_File(stream_type,filename,hashtable_matrix);
    }

    // #################################################################################################3
    //      Hashtable_To_Matlab
    // #################################################################################################3
    void Hashtable_To_Matlab(HASHTABLE<VECTOR<int,2>,T>& hashtable_matrix, const std::string& filename) {
        std::ofstream out(filename);
        for (HASHTABLE_ITERATOR<VECTOR<int,2>,T> iterator(hashtable_matrix);iterator.Valid();iterator.Next()) {
            out<<iterator.Key().x<<" "<<iterator.Key().y<<" "<<iterator.Data()<<std::endl;}
        out.close();

        // ** load into matlab with **
        // load matlab_formatted_matrix
        // S=spconvert(matlab_formatted_matrix);
    }

    // #################################################################################################3
    //      Matlab_To_Hashtable
    // #################################################################################################3
    void Matlab_To_Hashtable(HASHTABLE<VECTOR<int,2>,T>& hashtable_matrix, const std::string& filename) {
        int row,column;
        T data;
        std::ifstream in(filename,std::fstream::binary);
        while (in>>row  && in>>column && in>>data) {
            hashtable_matrix.Set(VECTOR<int,2>(row,column),data);}
    }

    // #################################################################################################3
    //      Hashtable_To_PARDISO
    // #################################################################################################3
    void Hashtable_To_PARDISO(HASHTABLE<VECTOR<int,2>,T>& hashtable_matrix, std::string suffix="") {
        // put the data into an array
        ARRAY<PAIR<VECTOR<int,2>,T> > array_matrix;
        for (HASHTABLE_ITERATOR<VECTOR<int,2>,T> iterator(hashtable_matrix);iterator.Valid();iterator.Next())
            if (iterator.Key().x<=iterator.Key().y) // use upper triangle part for PARDISO
                array_matrix.Append(PAIR<VECTOR<int,2>,T>(iterator.Key(),iterator.Data()));
        // sort the array
        Sort(array_matrix,Sort_Function<T>);

        // assemble the data into upper triangle compressed sparse row CSR format as required by PARDISO for symmetric matrices
        const int number_non_zeros=array_matrix.m;
        const int number_rows=array_matrix(array_matrix.m).x.x;

        // Matrix data
        T* a=new T[number_non_zeros];
        MKL_INT n = number_rows;
        MKL_INT* ia=new MKL_INT[number_rows+1];
        MKL_INT* ja=new MKL_INT[number_non_zeros];

        // fill in the matrix data
        int current_row=0;
        for (int i=1;i<=array_matrix.m;i++) {
            a[i-1]=array_matrix(i).y;
            ja[i-1]=array_matrix(i).x.y;
            while (array_matrix(i).x.x>current_row) {
                current_row++;
                ia[current_row-1]=i;}}
        ia[number_rows]=number_non_zeros+1;

        // write files
        {
        std::string ndat="n"+suffix+".dat";
        std::ofstream out(ndat,std::fstream::binary);
        out.write((const char*)(&n),sizeof(MKL_INT));
        out.close();
        }

        {
        std::string iadat="ia"+suffix+".dat";
        std::ofstream out(iadat,std::fstream::binary);
        out.write((const char*)ia,(n+1)*sizeof(MKL_INT));
        out.close();
        }

        {
        std::string jadat="ja"+suffix+".dat";
        std::ofstream out(jadat,std::fstream::binary);
        out.write((const char*)ja,(ia[n]-1)*sizeof(MKL_INT));
        out.close();
        }

        {
        std::string adat="a"+suffix+".dat";
        std::ofstream out(adat,std::fstream::binary);
        out.write((const char*)a,(ia[n]-1)*sizeof(T));
        out.close();
        }

        // epilogue
        delete a;
        delete ia;
        delete ja;
    }


};
}
#endif
