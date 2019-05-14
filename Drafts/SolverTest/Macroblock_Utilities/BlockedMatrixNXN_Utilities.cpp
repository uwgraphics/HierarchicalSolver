//#####################################################################
// Copyright 2015, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "BlockedMatrixNXN.h"
#include "BlockedMatrixNXN_Utilities.h"
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <stdlib.h>
using namespace PhysBAM;

//################################################
// Function Print_Matrix (for MATRIX_MXN)
//################################################
template<class T>
void Print_Matrix(const MATRIX_MXN<T>& A,std::ostream& out)
{
    PHYSBAM_ASSERT(A.m==A.n);
    int n=A.n;
    for(int i=1;i<=A.m;i++){
        for(int j=1;j<=A.n;j++){ out.width(12);out<<A(i,j);if(j<A.n) out<<" ";}
        out<<std::endl;}
}
//################################################
// Function Print_Matrix (for BlockedMatrixNXN)
//################################################
template<class T>
void BlockedMatrixNXN_Utilities::
Print_Matrix(const BlockedMatrixNXN<T>& A,std::ostream& out)
{
    ARRAY<std::string> rows(4*A.n);
    int index=0;
    for(int j=0;j<A.n;j++){
        for(int i=j+1;i<A.n;i++){
            const auto L=A.L(index++);
            for(int k=0;k<16;k++)
                rows(4*i+k%4+1)+=STRING_UTILITIES::string_sprintf("%12.6g ",L[k]);}
        const auto D=A.D(j);
        for(int k=0;k<16;k++)
            rows(4*j+k%4+1)+=STRING_UTILITIES::string_sprintf("%12.6g ",D[k]);}
    for(int i=1;i<=4*A.n;i++)
        out<<rows(i)<<std::endl;
}
//################################################
// Function Allocate
//################################################
template<class T>
void BlockedMatrixNXN_Utilities::
Allocate(const int block_dim,void *&D_raw,void *&L_raw)
{
    PHYSBAM_ASSERT(block_dim>0);
    posix_memalign(&D_raw, 64, block_dim*16*sizeof(T));
    posix_memalign(&L_raw, 64, block_dim*(block_dim-1)/2*16*sizeof(T));
}
//################################################
// Function LoadFromPhysbam
//################################################
template<class T>
void BlockedMatrixNXN_Utilities::
LoadFromPhysbam(BlockedMatrixNXN<T>& A,const MATRIX_MXN<T>& Aref)
{
    PHYSBAM_ASSERT(Aref.m==Aref.n && Aref.m==4*A.n);
    int index=0;
    for(int j=0;j<A.n;j++){
        for(int i=j+1;i<A.n;i++){
            auto L=A.L(index++);
            for(int k=0;k<16;k++)
                L[k]=Aref(4*i+k%4+1,4*j+k/4+1);}
        auto D=A.D(j);
        for(int k=0;k<16;k++)
            D[k]=Aref(4*j+k%4+1,4*j+k/4+1);}
}
//################################################

template void Print_Matrix(const MATRIX_MXN<float>& A,std::ostream& out);
namespace BlockedMatrixNXN_Utilities{
template void Print_Matrix(const BlockedMatrixNXN<float>& A,std::ostream& out);
template void Allocate<float>(const int block_dim,void *&D_raw,void *&L_raw);
template void LoadFromPhysbam(BlockedMatrixNXN<float>& A,const MATRIX_MXN<float>& Aref);
}
