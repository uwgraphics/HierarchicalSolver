//#####################################################################
// Copyright 2015, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BLOCK_ELASTICITY_OPERATIONS
//#####################################################################
#ifndef __BLOCK_ELASTICITY_OPERATIONS_h__
#define __BLOCK_ELASTICITY_OPERATIONS_h__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Common_Tools/Grids_Uniform_PDE_Linear/STENCIL.h>

using namespace PhysBAM;

template<int d>
struct BLOCK_ELASTICITY_OPERATIONS
{
private:
    typedef VECTOR<int,d> T_INDEX;
public:
    template<class T> using SYSTEM_MATRIX = ARRAY<STENCIL<T,d>,T_INDEX>;

    enum SIMD_CONV_MODE { COPY_UNIQUE, COPY_UNIQUE_CLEAR, COPY_DUPLICATE, ACCUM_ADD, ACCUM_SUB };
    enum SIMD_OP_MODE { PLUS_EQUAL, SUB_EQUAL };
    enum SIMD_UNARY_OP_MODE { NEGATE, ACCUMULATE_AND_DISTRIBUTE, COLLAPSE};

//#####################################################################
    static void Build_Association_Matrix(ARRAY<ARRAY<T_INDEX>,T_INDEX>& association_matrix,const RANGE<T_INDEX>& node_range);
    static void Initialize_Subdomains(ARRAY<RANGE<T_INDEX> >& subdomains,const bool verbose=false);
    static void Verify_Subdomains(const ARRAY<RANGE<T_INDEX> >& subdomains,
                                  ARRAY<ARRAY<T_INDEX>,T_INDEX> global_association_matrix,ARRAY<ARRAY<T_INDEX>,T_INDEX> subdomain_association_matrix);
    static void Identify_Duplicated_Variables(const ARRAY<RANGE<T_INDEX> >& subdomains,HASHTABLE<T_INDEX, ARRAY<PAIR<T_INDEX,int> > >& real_index_to_canonical_s_SIMD_indices,
                                              ARRAY<ARRAY<PAIR<T_INDEX,int> > >& duplicated_variables,const bool verbose=false);
    static void Identify_Duplicated_Variables_By_Level(const ARRAY<RANGE<T_INDEX> >& subdomains,VECTOR<ARRAY<ARRAY<PAIR<T_INDEX,int> > >,4>& duplicated_variables_by_level,
                                              ARRAY<int,T_INDEX> node_level, const bool verbose); 
    static void Reorder_Range(ARRAY<int,T_INDEX>& rank,ARRAY<int,T_INDEX>& level,const RANGE<T_INDEX>& coordinate_range,const VECTOR<int,2>& rank_range,
        const int current_level,const int min_level);
    template<class T,int SIMDw>
        static void SIMD_From_Unified(const ARRAY<RANGE<T_INDEX> >& subdomains,ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& SIMDx,const ARRAY<VECTOR<T,d>,T_INDEX>& x, const SIMD_CONV_MODE mode);
    template<class T,int SIMDw>
        static void SIMD_To_Unified(const ARRAY<RANGE<T_INDEX> >& subdomains,const ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& SIMDx,ARRAY<VECTOR<T,d>,T_INDEX>& x, const SIMD_CONV_MODE mode);
 template<class T,int SIMDw>
        static void SIMD_To_Unified2(const ARRAY<RANGE<T_INDEX> >& subdomains,
                                     const ARRAY<T_INDEX>& unified_index_list,
                                     const ARRAY<int,T_INDEX>& unified_index_id,
                                     const ARRAY<T_INDEX>& subdomain_index_list,
                                     const ARRAY<int,T_INDEX>& subdomain_index_id,
                                     const ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& SIMDx,ARRAY<VECTOR<T,d>,T_INDEX>& x, const SIMD_CONV_MODE mode);
    template<class T,int SIMDw>
        static void SIMD_OP_SIMD(const ARRAY<RANGE<T_INDEX> >& subdomains, ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& A, const ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& B, const SIMD_OP_MODE mode); 
    template<class T,int SIMDw>
        static void SIMD_UNARY_OP(const ARRAY<RANGE<T_INDEX> >& subdomains, ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& A, const SIMD_UNARY_OP_MODE mode); 
    template<class T,int SIMDw>
        static void SIMD_UNARY_OP(const ARRAY<RANGE<T_INDEX> >& subdomains, ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& A, const ARRAY<ARRAY<PAIR<T_INDEX,int> > >& duplicated_variables, const SIMD_UNARY_OP_MODE mode); 
    template<class T,int SIMDw>
        static void SIMD_UNARY_OP(const ARRAY<RANGE<T_INDEX> >& subdomains, ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& A, const int level,
               const VECTOR<ARRAY<ARRAY<PAIR<T_INDEX,int> > >,4>& duplicated_variables_by_level, const SIMD_UNARY_OP_MODE mode);
    template<class T,int SIMDw>
        static void SIMD_Prune_Duplicates(const ARRAY<RANGE<T_INDEX> >& subdomains,ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& SIMDx);
    template<class T,int SIMDw>
        static void SIMD_Dirty_Regions(const ARRAY<RANGE<T_INDEX> >& subdomains, const ARRAY< RANGE<T_INDEX>, VECTOR<int,1> >& regions, const ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& SIMDx);
//#####################################################################
};
#endif
