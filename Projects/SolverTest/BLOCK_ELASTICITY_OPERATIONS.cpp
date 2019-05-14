//#####################################################################
// Copyright 2015, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include "BLOCK_ELASTICITY_OPERATIONS.h"
#include "ARBITRARY_RANGE_ITERATOR.h"
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>

//#####################################################################
// Function Build_Association_Matrix
//#####################################################################
template<int d> void BLOCK_ELASTICITY_OPERATIONS<d>::
Build_Association_Matrix(ARRAY<ARRAY<T_INDEX>,T_INDEX>& association_matrix,const RANGE<T_INDEX>& node_range)
{
    LOG::SCOPE scope("BLOCK_ELASTICITY_OPERATIONS::Build_Association_Matrix");

    association_matrix.Resize(node_range);
    
    for(RANGE_ITERATOR<d> node_iterator(node_range);node_iterator.Valid();node_iterator.Next()){
        const T_INDEX& node_index=node_iterator.Index();
        association_matrix(node_index).Remove_All();
        for(RANGE_ITERATOR<d> other_node_iterator(RANGE<T_INDEX>::Intersect(RANGE<T_INDEX>(node_index-1,node_index+1),node_range));
            other_node_iterator.Valid();other_node_iterator.Next())
            association_matrix(node_index).Append(other_node_iterator.Index());}
}
//#####################################################################
// Function Initialize_Subdomains
//#####################################################################
namespace{
template<int d>
RANGE<VECTOR<int,d> > Duplicate_Range(const RANGE<VECTOR<int,d> >& range,const int axis)
{
    RANGE<VECTOR<int,d> > result(range);
    result.max_corner(axis) = range.max_corner(axis)*2;
    return result;
}
template<int d>
RANGE<VECTOR<int,d> > Reflect_Range(const RANGE<VECTOR<int,d> >& range,const RANGE<VECTOR<int,d> > aggregate_range,const int axis)
{
    RANGE<VECTOR<int,d> > result(range);
    result.min_corner(axis) = aggregate_range.max_corner(axis)*2-range.min_corner(axis);
    result.max_corner(axis) = aggregate_range.max_corner(axis)*2-range.max_corner(axis);
    return result;
}
}

template<int d> void BLOCK_ELASTICITY_OPERATIONS<d>::
Initialize_Subdomains(ARRAY<RANGE<T_INDEX> >& subdomains,const bool verbose)
{
    static_assert(d==3,"Initialize_Subdomains only supports the 3D case");

    LOG::SCOPE scope("BLOCK_ELASTICITY_OPERATIONS::Initialize_Subdomains");

    RANGE<T_INDEX> base_subdomain(T_INDEX(),T_INDEX::All_Ones_Vector()*4);
        
    ARRAY<int> reflection_sequence;
    // 8 cells per macroblock
    // reflection_sequence.Append_Elements(VECTOR<int,3>(3,2,1));
    // 16 cells per macroblock
    reflection_sequence.Append_Elements(VECTOR<int,4>(3,2,1,1));
        
    subdomains.Append(base_subdomain);
        
    RANGE<T_INDEX> aggregate_domain(base_subdomain);
    for(int reflection=1;reflection<=reflection_sequence.m;reflection++){
        int axis=reflection_sequence(reflection);
        ARRAY<RANGE<T_INDEX> > reflected_subdomains;

        for(int s=1;s<=subdomains.m;s++){
            reflected_subdomains.Append(subdomains(s));
            reflected_subdomains.Append(Reflect_Range(subdomains(s),aggregate_domain,axis));}

        subdomains.Exchange(reflected_subdomains);
        aggregate_domain=Duplicate_Range(aggregate_domain,axis);}

    if(verbose)
        LOG::cout<<"subdomains = "<<subdomains<<std::endl;
}
//#####################################################################
// Function Identify_Duplicated_Variables
//#####################################################################
template<int d> void BLOCK_ELASTICITY_OPERATIONS<d>::
Identify_Duplicated_Variables(
    const ARRAY<RANGE<T_INDEX> >& subdomains,
    HASHTABLE<T_INDEX, ARRAY<PAIR<T_INDEX,int> > >& real_index_to_canonical_s_SIMD_indices,
    ARRAY<ARRAY<PAIR<T_INDEX,int> > >& duplicated_variables,
    const bool verbose) 
{
    // entries in real_index_to_canonical_s_SIMD_indices map the real indices to an array of all 'multiplexed' SIMD variables that correspond to the same variable as pairs of canonical indices and subdomain id
    // each entry in duplicated_variables is a list of 'multiplexed' SIMD variables that correspond to the same unified variable (identified as a pair of canonical_index and subdomain_id)
    ARRAY<ARBITRARY_RANGE_ITERATOR<d>*> subdomain_iterators;
    for(int s=1;s<=subdomains.m;s++) subdomain_iterators.Append(new ARBITRARY_RANGE_ITERATOR<d>(subdomains(s)));
    while(subdomain_iterators(1)->Valid()){
        for(int s=1;s<=subdomains.m;s++) {
            const T_INDEX& real_index=subdomain_iterators(s)->Index();
            const T_INDEX& canonical_index=subdomain_iterators(1)->Index();
            real_index_to_canonical_s_SIMD_indices.Get_Or_Insert(real_index).Append(PAIR<T_INDEX,int>(canonical_index,s));}
        for(int s=1;s<=subdomains.m;s++) subdomain_iterators(s)->Next();
    }
    // Generate list of duplicated variables
    for (HASHTABLE_ITERATOR<T_INDEX, ARRAY<PAIR<T_INDEX,int> > > iterator(real_index_to_canonical_s_SIMD_indices);iterator.Valid();iterator.Next()){
        ARRAY<PAIR<T_INDEX,int> >& data=iterator.Data();
        if (data.m>1) duplicated_variables.Append(data);
    }    
}

//#####################################################################
// Function Identify_Duplicated_Variables_By_Level
//#####################################################################
template<int d> void BLOCK_ELASTICITY_OPERATIONS<d>::
Identify_Duplicated_Variables_By_Level(
    const ARRAY<RANGE<T_INDEX> >& subdomains,
    VECTOR<ARRAY<ARRAY<PAIR<T_INDEX,int> > >,4>& duplicated_variables_by_level,
    ARRAY<int,T_INDEX> node_level,
    const bool verbose) 
{
    

    HASHTABLE<T_INDEX, ARRAY<PAIR<T_INDEX,int> > > real_index_to_canonical_s_SIMD_indices;
    // entries in real_index_to_canonical_s_SIMD_indices map the real indices to an array of all 'multiplexed' SIMD variables that correspond to the same variable as pairs of canonical indices and subdomain id
    // each entry in duplicated_variables is a list of 'multiplexed' SIMD variables that correspond to the same unified variable (identified as a pair of canonical_index and subdomain_id)
    ARRAY<ARBITRARY_RANGE_ITERATOR<d>*> subdomain_iterators;
    
    for(int s=1;s<=subdomains.m;s++) subdomain_iterators.Append(new ARBITRARY_RANGE_ITERATOR<d>(subdomains(s)));
    while(subdomain_iterators(1)->Valid()){
        for(int s=1;s<=subdomains.m;s++) {
            const T_INDEX& real_index=subdomain_iterators(s)->Index();
            const T_INDEX& canonical_index=subdomain_iterators(1)->Index();
            real_index_to_canonical_s_SIMD_indices.Get_Or_Insert(real_index).Append(PAIR<T_INDEX,int>(canonical_index,s));}
        for(int s=1;s<=subdomains.m;s++) subdomain_iterators(s)->Next();
    }
    
    // Generate list of duplicated variables divided by level
    for (HASHTABLE_ITERATOR<T_INDEX, ARRAY<PAIR<T_INDEX,int> > > iterator(real_index_to_canonical_s_SIMD_indices);iterator.Valid();iterator.Next()){

        const ARRAY<PAIR<T_INDEX,int> >& data=iterator.Data();
        const T_INDEX& real_index=iterator.Key();
        const int level=node_level(real_index);
            if (data.m>1&&level>=1) duplicated_variables_by_level(level).Append(data);
    }    
    
}
//#####################################################################
// Function Verify_Subdomains
//#####################################################################
template<int d> void BLOCK_ELASTICITY_OPERATIONS<d>::
Verify_Subdomains(const ARRAY<RANGE<T_INDEX> >& subdomains,
    ARRAY<ARRAY<T_INDEX>,T_INDEX> global_association_matrix,
    ARRAY<ARRAY<T_INDEX>,T_INDEX> subdomain_association_matrix)
{
    LOG::SCOPE scope("BLOCK_ELASTICITY_OPERATIONS::Verify_Subdomains");

    ARRAY<ARBITRARY_RANGE_ITERATOR<d>*> range_iterators;
    for(int s=1;s<=subdomains.m;s++)
        range_iterators.Append(new ARBITRARY_RANGE_ITERATOR<d>(subdomains(s)));
        
    for(;range_iterators(1)->Valid();range_iterators(1)->Next())
        for(int s=2;s<=subdomains.m;s++){
            ARBITRARY_RANGE_ITERATOR<d> other_range_iterator_primary(subdomains(1));
            ARBITRARY_RANGE_ITERATOR<d> other_range_iterator_secondary(subdomains(s));
            for(;other_range_iterator_primary.Valid();other_range_iterator_primary.Next(),other_range_iterator_secondary.Next())
                if(subdomain_association_matrix(range_iterators(1)->Index()).Contains(other_range_iterator_primary.Index()))
                    PHYSBAM_ASSERT(global_association_matrix(range_iterators(s)->Index()).Contains(other_range_iterator_secondary.Index()));
            range_iterators(s)->Next();}

    for(int s=2;s<=subdomains.m;s++)
        PHYSBAM_ASSERT(!range_iterators(s)->Valid());

    range_iterators.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Reorder_Range
//#####################################################################
template<int d> void BLOCK_ELASTICITY_OPERATIONS<d>::
Reorder_Range(ARRAY<int,T_INDEX>& rank,ARRAY<int,T_INDEX>& level,const RANGE<VECTOR<int,d> >& coordinate_range,const VECTOR<int,2>& rank_range,
    const int current_level,const int min_level)
{
    // Split along largest axis (or first one if several equal)
    int axis=(abs(coordinate_range.max_corner-coordinate_range.min_corner)+1).Arg_Abs_Max();
    // If down to a single index, mark it and exit
    if(coordinate_range.Edge_Lengths()==T_INDEX()){
        PHYSBAM_ASSERT(rank_range.x==rank_range.y);
        PHYSBAM_ASSERT(current_level==min_level);
        rank(coordinate_range.min_corner)=rank_range(1);
        level(coordinate_range.min_corner)=current_level;
        return;}

    RANGE<T_INDEX> left_coordinate_range(coordinate_range),right_coordinate_range(coordinate_range),interface_coordinate_range(coordinate_range);

    int bias=(coordinate_range.min_corner(axis)<coordinate_range.max_corner(axis))?1:-1;

    left_coordinate_range.max_corner(axis)=coordinate_range.min_corner(axis);
    left_coordinate_range.min_corner(axis)=(coordinate_range.min_corner(axis)+coordinate_range.max_corner(axis))/2-bias;
    right_coordinate_range.max_corner(axis)=coordinate_range.max_corner(axis);
    right_coordinate_range.min_corner(axis)=(coordinate_range.min_corner(axis)+coordinate_range.max_corner(axis))/2+bias;
    interface_coordinate_range.min_corner(axis)=interface_coordinate_range.max_corner(axis)=(coordinate_range.min_corner(axis)+coordinate_range.max_corner(axis))/2;

    int left_size = (abs( left_coordinate_range.min_corner - left_coordinate_range.max_corner ) +1 ).Product();
    int right_size = (abs( right_coordinate_range.min_corner - right_coordinate_range.max_corner ) +1 ).Product();
    int interface_size = (abs( interface_coordinate_range.min_corner - interface_coordinate_range.max_corner ) +1 ).Product();

    VECTOR<int, 2> left_rank_range( rank_range(1), rank_range(1)+left_size-1 );
    VECTOR<int, 2> right_rank_range( rank_range(1)+left_size, rank_range(1)+left_size+right_size-1 );
    VECTOR<int, 2> interface_rank_range( rank_range(1)+left_size+right_size, rank_range(2) );

    Reorder_Range(rank,level,left_coordinate_range,left_rank_range,std::max(current_level-1,min_level),min_level);
    Reorder_Range(rank,level,right_coordinate_range,right_rank_range,std::max(current_level-1,min_level),min_level);
    Reorder_Range(rank,level,interface_coordinate_range,interface_rank_range,current_level,current_level);

#if 0

    //LOG::cout << "Reordering range (L) : " << left_coordinate_range << ", ( " << left_rank_range << " )" << std::endl;
    Reorder_Range( reordering, left_coordinate_range, left_rank_range, depth+1 );    

    if( left_size == reordering.subdomain_unit*reordering.subdomain_unit*reordering.subdomain_unit ) {
        reordering.subdomains.Append( HASHTABLE<int,VECTOR<int,d> >() );
        //LOG::cout << "Subdomain ("<< reordering.subdomains.m << ") : " << std::endl;
        //LOG::cout << "Range: " << left_coordinate_range << std::endl;
        for( ARBITRARY_RANGE_ITERATOR<d> iter(left_coordinate_range); iter.Valid(); iter.Next() ){
            //LOG::cout << iter.Index() << ":  " << rank(iter.Index()) << std::endl;
            reordering.subdomains.Last().Insert( rank(iter.Index()), iter.Index() );
        }
    }
    //LOG::cout << "Reordering range (R): " << right_coordinate_range << ", ( " << right_rank_range << " )" << std::endl;
    Reorder_Range( reordering, right_coordinate_range, right_rank_range, depth+1 );

    if( right_size == reordering.subdomain_unit*reordering.subdomain_unit*reordering.subdomain_unit ) {
        reordering.subdomains.Append( HASHTABLE<int,VECTOR<int,d> >() );
        //LOG::cout << "Subdomain ("<< reordering.subdomains.m << ") : " << std::endl;
        //LOG::cout << "Range: " << right_coordinate_range << std::endl;
        for( ARBITRARY_RANGE_ITERATOR<d> iter(right_coordinate_range); iter.Valid(); iter.Next() ){
            //LOG::cout << iter.Index() << ":  " << rank(iter.Index()) << std::endl;
            reordering.subdomains.Last().Insert( rank(iter.Index()), iter.Index() );
        }
    }

    //LOG::cout << "Reordering range (I): " << interface_coordinate_range << ", ( " << interface_rank_range << " )" << std::endl;
    Reorder_Range( reordering, interface_coordinate_range, interface_rank_range, depth+1 );

    if( left_size >= reordering.subdomain_unit*reordering.subdomain_unit*reordering.subdomain_unit && 
        interface_size >= reordering.subdomain_unit*reordering.subdomain_unit ) {       
        reordering.interfaces(depth).Append( HASHTABLE<int,VECTOR<int,d> >() );
        //LOG::cout << "Interface ("<< reordering.interfaces(depth).m << ") Level " << depth << " : " << std::endl;
        //LOG::cout << "Range: " << interface_coordinate_range << std::endl;
        for( ARBITRARY_RANGE_ITERATOR<d> iter(interface_coordinate_range); iter.Valid(); iter.Next() ){
            //LOG::cout << iter.Index() << ":  " << rank(iter.Index()) << std::endl;
            reordering.interfaces(depth).Last().Insert( rank(iter.Index()), iter.Index() );
        }
    }
    
#endif


}
//#####################################################################
// Function SIMD_From_Unified
//#####################################################################
template<int d> template<class T,int SIMDw>
void BLOCK_ELASTICITY_OPERATIONS<d>::
SIMD_From_Unified(const ARRAY<RANGE<T_INDEX> >& subdomains,ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& SIMDx,const ARRAY<VECTOR<T,d>,T_INDEX>& x,const SIMD_CONV_MODE mode)
{
    HASHTABLE<T_INDEX> value_consumed;
    ARRAY<ARBITRARY_RANGE_ITERATOR<d>*> subdomain_iterators;
    for(int s=1;s<=subdomains.m;s++) subdomain_iterators.Append(new ARBITRARY_RANGE_ITERATOR<d>(subdomains(s)));
    while(subdomain_iterators(1)->Valid()){

        // Check for consistency
        for(int s=2;s<=subdomains.m;s++) PHYSBAM_ASSERT(subdomain_iterators(s)->Valid());

        switch( mode ){
        case COPY_UNIQUE:
            {
                for(int s=1;s<=subdomains.m;s++){
                    if( !value_consumed.Contains(subdomain_iterators(s)->Index())){
                        for(int v=1;v<=d;v++)
                            SIMDx(subdomain_iterators(1)->Index())(v)(s)=x(subdomain_iterators(s)->Index())(v);
                        value_consumed.Set(subdomain_iterators(s)->Index());
                    }
                }
            }
            break;
        case COPY_UNIQUE_CLEAR:
            {
                for(int s=1;s<=subdomains.m;s++){
                    if( !value_consumed.Contains(subdomain_iterators(s)->Index())){
                        for(int v=1;v<=d;v++)
                            SIMDx(subdomain_iterators(1)->Index())(v)(s)=x(subdomain_iterators(s)->Index())(v);
                        value_consumed.Set(subdomain_iterators(s)->Index());
                    }
                    else{
                        for(int v=1;v<=d;v++)
                            SIMDx(subdomain_iterators(1)->Index())(v)(s)=T();
                    }
                }
            }
            break;
        case COPY_DUPLICATE:
            {
                for(int s=1;s<=subdomains.m;s++){
                    for(int v=1;v<=d;v++)
                        SIMDx(subdomain_iterators(1)->Index())(v)(s)=x(subdomain_iterators(s)->Index())(v);
                }
            }
            break;
        }      
        
        // Increment everything
        for(int s=1;s<=subdomains.m;s++) subdomain_iterators(s)->Next();

    }

    // Check for consistency
    for(int s=2;s<=subdomains.m;s++) PHYSBAM_ASSERT(!subdomain_iterators(s)->Valid());
}
//#####################################################################
// Function SIMD_To_Unified
//#####################################################################
template<int d> template<class T,int SIMDw>
void BLOCK_ELASTICITY_OPERATIONS<d>::
SIMD_To_Unified(const ARRAY<RANGE<T_INDEX> >& subdomains,const ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& SIMDx,ARRAY<VECTOR<T,d>,T_INDEX>& x,const SIMD_CONV_MODE mode)
{
    ARRAY<ARBITRARY_RANGE_ITERATOR<d>*> subdomain_iterators;
    for(int s=1;s<=subdomains.m;s++) subdomain_iterators.Append(new ARBITRARY_RANGE_ITERATOR<d>(subdomains(s)));
    while(subdomain_iterators(1)->Valid()){

        // Check for consistency
        for(int s=2;s<=subdomains.m;s++) PHYSBAM_ASSERT(subdomain_iterators(s)->Valid());

        // Copy to unified array
        for(int s=1;s<=subdomains.m;s++)
            for(int v=1;v<=d;v++){
                if(mode==ACCUM_SUB)
                    x(subdomain_iterators(s)->Index())(v) -= SIMDx(subdomain_iterators(1)->Index())(v)(s);
                if(mode==ACCUM_ADD)
                    x(subdomain_iterators(s)->Index())(v) += SIMDx(subdomain_iterators(1)->Index())(v)(s);
            }

        // Increment everything
        for(int s=1;s<=subdomains.m;s++) subdomain_iterators(s)->Next();

    }

    // Check for consistency
    for(int s=2;s<=subdomains.m;s++) PHYSBAM_ASSERT(!subdomain_iterators(s)->Valid());
}
//#####################################################################
// Function SIMD_To_Unified2
//#####################################################################
template<int d> template<class T,int SIMDw>
void BLOCK_ELASTICITY_OPERATIONS<d>::
SIMD_To_Unified2(const ARRAY<RANGE<T_INDEX> >& subdomains,
                 const ARRAY<T_INDEX>& unified_index_list,
                 const ARRAY<int,T_INDEX>& unified_index_id,
                 const ARRAY<T_INDEX>& subdomain_index_list,
                 const ARRAY<int,T_INDEX>& subdomain_index_id,
                 

                 const ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& SIMDx,
                 ARRAY<VECTOR<T,d>,T_INDEX>& x,
                 const SIMD_CONV_MODE mode)
{
    int counter = 0;
    ARRAY<ARBITRARY_RANGE_ITERATOR<d>*> subdomain_iterators;
    for(int s=1;s<=subdomains.m;s++) subdomain_iterators.Append(new ARBITRARY_RANGE_ITERATOR<d>(subdomains(s)));
    while(subdomain_iterators(1)->Valid()){

        // Check for consistency
        for(int s=2;s<=subdomains.m;s++) PHYSBAM_ASSERT(subdomain_iterators(s)->Valid());

        // Copy to unified array
        for(int s=1;s<=subdomains.m;s++)
            if( unified_index_id(subdomain_iterators(s)->Index()) != 0 ){
                LOG::cout << s <<"  S["<< std::setfill('0') << std::setw(2) <<unified_index_id(subdomain_iterators(s)->Index()) << "] += ";
                LOG::cout << "A["<<subdomain_index_id(subdomain_iterators(1)->Index())<<"]["<<s<<"];  "<</*subdomain_iterators(1)->Index() <<*/ std::endl;
                    
            }
        //for(int v=1;v<=d;v++){                
        //        x(subdomain_iterators(s)->Index())(v) -= SIMDx(subdomain_iterators(1)->Index())(v)(s);
        //    }

        // Increment everything
        for(int s=1;s<=subdomains.m;s++) subdomain_iterators(s)->Next();

    }

    // Check for consistency
    for(int s=2;s<=subdomains.m;s++) PHYSBAM_ASSERT(!subdomain_iterators(s)->Valid());
}
//#####################################################################
//  Function SIMD_OP_SIMD
//#####################################################################

template<int d> template<class T,int SIMDw>
void BLOCK_ELASTICITY_OPERATIONS<d>::
SIMD_OP_SIMD(const ARRAY<RANGE<T_INDEX> >& subdomains, ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& A, const ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& B, const SIMD_OP_MODE mode){

    typedef VECTOR<T,d> TV;
    for( ARBITRARY_RANGE_ITERATOR<d> subdomain_iterator( subdomains(1)); subdomain_iterator.Valid(); subdomain_iterator.Next() ){
        for(int v=1;v<=d;v++){
            switch( mode ){
            case PLUS_EQUAL:
                {
                    A(subdomain_iterator.Index())(v) += B(subdomain_iterator.Index())(v);
                }
                break;    
            case SUB_EQUAL:
                {
                    A(subdomain_iterator.Index())(v) -= B(subdomain_iterator.Index())(v);
                }
                break;  
            //case NEGATE:
            //    {
            //        A(subdomain_iterator.Index())(v) = -B(subdomain_iterator.Index())(v);
            //    }
            //    break; 
            }
        }
    }
}
//#####################################################################
//  Function SIMD_UNARY_OP with duplicated_variables_by_level
//#####################################################################
template<int d> template<class T,int SIMDw>
void BLOCK_ELASTICITY_OPERATIONS<d>::
SIMD_UNARY_OP(const ARRAY<RANGE<T_INDEX> >& subdomains, ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& A, const int level, const VECTOR<ARRAY<ARRAY<PAIR<T_INDEX,int> > >,4>& duplicated_variables_by_level, const SIMD_UNARY_OP_MODE mode)
{
    typedef VECTOR<T,d> TV;
    switch( mode ){
    case ACCUMULATE_AND_DISTRIBUTE:
        {   
            const ARRAY<ARRAY<PAIR<T_INDEX,int> > >& duplicated_variables=duplicated_variables_by_level(level);
            for (int i=1;i<=duplicated_variables.m;i++) {
                const ARRAY<PAIR<T_INDEX,int> >& array=duplicated_variables(i);
                // accumulate
                TV sum;
                for (int j=1;j<=array.m;j++) {
                    const T_INDEX& index=array(j).x;
                    const int s=array(j).y;
                    for (int v=1;v<=3;v++) {
                        sum(v)+=A(index)(v)(s);}}
                // distrubute
                for (int j=1;j<=array.m;j++) {
                    const T_INDEX& index=array(j).x;
                    const int s=array(j).y;
                    for (int v=1;v<=3;v++) {
                        A(index)(v)(s)=sum(v);}}
            }
        }
        break;
    case COLLAPSE:
        {   
            const ARRAY<ARRAY<PAIR<T_INDEX,int> > >& duplicated_variables=duplicated_variables_by_level(level);
            for (int i=1;i<=duplicated_variables.m;i++) {
                const ARRAY<PAIR<T_INDEX,int> >& array=duplicated_variables(i);
                // accumulate
                TV sum;
                for (int j=1;j<=array.m;j++) {
                    const T_INDEX& index=array(j).x;
                    const int s=array(j).y;
                    for (int v=1;v<=3;v++) {
                        sum(v)+=A(index)(v)(s);
                        A(index)(v)(s)=(T)0.;}}
                // assign uniquely
                    const T_INDEX& index=array(1).x;
                    const int s=array(1).y;
                    for (int v=1;v<=3;v++) {
                        A(index)(v)(s)=sum(v);}}
            }
        break;
    case NEGATE:
        {
            for( ARBITRARY_RANGE_ITERATOR<d> subdomain_iterator( subdomains(1)); subdomain_iterator.Valid(); subdomain_iterator.Next() ){
                for(int v=1;v<=d;v++){
                    A(subdomain_iterator.Index())(v)*=(T)-1.;
                }
            }
        }
        break;
    }
}




//#####################################################################
//  Function SIMD_UNARY_OP with duplicated_variables
//#####################################################################
template<int d> template<class T,int SIMDw>
void BLOCK_ELASTICITY_OPERATIONS<d>::
SIMD_UNARY_OP(const ARRAY<RANGE<T_INDEX> >& subdomains, ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& A, const ARRAY<ARRAY<PAIR<T_INDEX,int> > >& duplicated_variables, const SIMD_UNARY_OP_MODE mode){

    typedef VECTOR<T,d> TV;
    switch( mode ){
    case ACCUMULATE_AND_DISTRIBUTE:
        {
            for (int i=1;i<=duplicated_variables.m;i++) {
                const ARRAY<PAIR<T_INDEX,int> >& array=duplicated_variables(i);
                // accumulate
                TV sum;
                for (int j=1;j<=array.m;j++) {
                    const T_INDEX& index=array(j).x;
                    const int s=array(j).y;
                    for (int v=1;v<=3;v++) {
                        sum(v)+=A(index)(v)(s);}}
                // distrubute
                for (int j=1;j<=array.m;j++) {
                    const T_INDEX& index=array(j).x;
                    const int s=array(j).y;
                    for (int v=1;v<=3;v++) {
                        A(index)(v)(s)=sum(v);}}
            }
        }
        break;
    case NEGATE:
        {
            for( ARBITRARY_RANGE_ITERATOR<d> subdomain_iterator( subdomains(1)); subdomain_iterator.Valid(); subdomain_iterator.Next() ){
                for(int v=1;v<=d;v++){
                    A(subdomain_iterator.Index())(v)*=(T)-1.;
                }
            }
        }
        break;
        }
}
    //#####################################################################
    //  Function SIMD_UNARY_OP
    //#####################################################################
template<int d> template<class T,int SIMDw>
void BLOCK_ELASTICITY_OPERATIONS<d>::
SIMD_UNARY_OP(const ARRAY<RANGE<T_INDEX> >& subdomains, ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& A, const SIMD_UNARY_OP_MODE mode){

    typedef VECTOR<T,d> TV;
    switch( mode ){
        // this ACCUMULATE_AND_DISTRIBUTE is very slow.  It is only here for debugging purposes.
        case ACCUMULATE_AND_DISTRIBUTE:
            {
                HASHTABLE<T_INDEX, TV> real_index_to_sum;
                ARRAY<ARBITRARY_RANGE_ITERATOR<d>*> subdomain_iterators;
                for(int s=1;s<=subdomains.m;s++) subdomain_iterators.Append(new ARBITRARY_RANGE_ITERATOR<d>(subdomains(s)));
                // Accumulate
                while(subdomain_iterators(1)->Valid()){
                    for(int s=1;s<=subdomains.m;s++) {
                        TV value;
                        const T_INDEX& real_index=subdomain_iterators(s)->Index();
                        const T_INDEX& canonical_index=subdomain_iterators(1)->Index();
                        for(int v=1;v<=d;v++){
                            value(v) = A(canonical_index)(v)(s);
                        }
                        TV temp=real_index_to_sum.Get_Or_Insert(real_index, TV());
                        real_index_to_sum.Set(real_index, temp+value);
                    }
                    // Increment all iterators
                    for(int s=1;s<=subdomains.m;s++) subdomain_iterators(s)->Next();
                }
                for(int s=1;s<=subdomains.m;s++) subdomain_iterators(s)->Reset();
                // Distribute
                while(subdomain_iterators(1)->Valid()){
                    for(int s=1;s<=subdomains.m;s++) {
                        const T_INDEX& real_index=subdomain_iterators(s)->Index();
                        const T_INDEX& canonical_index=subdomain_iterators(1)->Index();
                        const TV& value=real_index_to_sum.Get(real_index);
                        for(int v=1;v<=d;v++){
                            A(canonical_index)(v)(s)=value(v);
                        }
                    }
                    // Increment all iterators
                    for(int s=1;s<=subdomains.m;s++) subdomain_iterators(s)->Next();
                }
            }
            break;
        case NEGATE:
        {
            for( ARBITRARY_RANGE_ITERATOR<d> subdomain_iterator( subdomains(1)); subdomain_iterator.Valid(); subdomain_iterator.Next() ){
                for(int v=1;v<=d;v++){
                    A(subdomain_iterator.Index())(v)*=(T)-1.;
                }
            }
        }
        break;
    }
}

//#####################################################################
// Function SIMD_Prune_Duplicates
//#####################################################################
template<int d> template<class T,int SIMDw>
void BLOCK_ELASTICITY_OPERATIONS<d>::
SIMD_Prune_Duplicates(const ARRAY<RANGE<T_INDEX> >& subdomains,ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& SIMDx)
{
    HASHTABLE<T_INDEX> value_consumed;
    ARRAY<ARBITRARY_RANGE_ITERATOR<d>*> subdomain_iterators;
    for(int s=1;s<=subdomains.m;s++) subdomain_iterators.Append(new ARBITRARY_RANGE_ITERATOR<d>(subdomains(s)));
    while(subdomain_iterators(1)->Valid()){

        // Check for consistency
        for(int s=2;s<=subdomains.m;s++) PHYSBAM_ASSERT(subdomain_iterators(s)->Valid());
        
        for(int s=1;s<=subdomains.m;s++){
            if( !value_consumed.Contains(subdomain_iterators(s)->Index()))
                value_consumed.Set(subdomain_iterators(s)->Index());
            else
                for(int v=1;v<=d;v++)
                    SIMDx(subdomain_iterators(1)->Index())(v)(s)=T();            
        }
        
        // Increment everything
        for(int s=1;s<=subdomains.m;s++) subdomain_iterators(s)->Next();

    }

    // Check for consistency
    for(int s=2;s<=subdomains.m;s++) PHYSBAM_ASSERT(!subdomain_iterators(s)->Valid());
}
//#####################################################################
// Function SIMD_Dirty_Regions
//#####################################################################
template<int d> template<class T,int SIMDw>
void BLOCK_ELASTICITY_OPERATIONS<d>::
SIMD_Dirty_Regions(const ARRAY<RANGE<T_INDEX> >& subdomains, const ARRAY< RANGE<T_INDEX>, VECTOR<int,1> >& regions, const ARRAY<VECTOR<VECTOR<T,SIMDw>,d>,T_INDEX>& SIMDx){

    typedef VECTOR<T,d> TV;

    LOG::cout << std::left << std::setw(12) << std::setfill(' ') << "";
    LOG::cout << std::left << std::setw(12) << std::setfill(' ') << "STATUS";
    LOG::cout << std::left << std::setw(12) << std::setfill(' ') << "UNIQUE";
    LOG::cout << std::left << std::setw(12) << std::setfill(' ') << "DUPLICATED";
    LOG::cout << std::left << std::setw(12) << std::setfill(' ') << "DIRTY" << std::endl;

    for( int region=0; region<5; region++){
        HASHTABLE<T_INDEX, TV> index_to_value;
        bool region_empty = true;
        bool region_duplicated = true;
        bool region_unique = true;
        
        ARBITRARY_RANGE_ITERATOR<d> region_iterator( regions(region) );
        while(region_iterator.Valid()){
            const T_INDEX& canonical_index = region_iterator.Index();
            for(int s=1;s<=subdomains.m;s++){

                // Find real index;
                T_INDEX real_index;
                ARRAY<ARBITRARY_RANGE_ITERATOR<d>*> subdomain_iterators;
                subdomain_iterators.Append(new ARBITRARY_RANGE_ITERATOR<d>(subdomains(1)));
                subdomain_iterators.Append(new ARBITRARY_RANGE_ITERATOR<d>(subdomains(s)));
                while( subdomain_iterators(1)->Valid() ){
                    if( subdomain_iterators(1)->Index() == canonical_index ){
                        real_index = subdomain_iterators(2)->Index();
                        break;
                    }
                    subdomain_iterators(1)->Next();
                    subdomain_iterators(2)->Next();
                }                

                TV value;
                for(int v=1;v<=d;v++){
                    value(v) = SIMDx(canonical_index)(v)(s);
                    if( SIMDx(canonical_index)(v)(s) != T() )
                        region_empty = false;
                }
                //if( region == 0 )
                //    LOG::cout << real_index << std::endl;

                if( value != TV() ){
                    if( index_to_value.Contains( real_index ) )
                        region_unique = false;
                    TV test_value = index_to_value.Get_Or_Insert( real_index, value );
                    if( test_value != value )
                        region_duplicated = false;                    
                }

            }  
            region_iterator.Next();
        }
        if(region_unique && region_duplicated)
            region_duplicated = false;
        if( region_empty ){
            region_unique = false;
            region_duplicated = false;
        }

        LOG::cout << "Region " << region << ":   ";
        LOG::cout << std::boolalpha << std::left << std::setw(12) << std::setfill(' ') << region_empty;
        LOG::cout << std::boolalpha << std::left << std::setw(12) << std::setfill(' ') << region_unique;
        LOG::cout << std::boolalpha << std::left << std::setw(12) << std::setfill(' ') << region_duplicated;
        LOG::cout << std::boolalpha << std::left << std::setw(12) << std::setfill(' ') << (!region_empty && !(region_unique || region_duplicated )) << std::endl;
        
    }
}

//#####################################################################

template class BLOCK_ELASTICITY_OPERATIONS<3>;
template void BLOCK_ELASTICITY_OPERATIONS<3>::SIMD_From_Unified<float,16>(const ARRAY<RANGE<VECTOR<int,3> > >&,
                                                                          ARRAY<VECTOR<VECTOR<float,16>,3>,VECTOR<int,3> >&,const ARRAY<VECTOR<float,3>,VECTOR<int,3> >&,const SIMD_CONV_MODE);
template void BLOCK_ELASTICITY_OPERATIONS<3>::SIMD_To_Unified<float,16>(const ARRAY<RANGE<VECTOR<int,3> > >&,
                                                                        const ARRAY<VECTOR<VECTOR<float,16>,3>,VECTOR<int,3> >&,ARRAY<VECTOR<float,3>,VECTOR<int,3> >&, const SIMD_CONV_MODE );
template void BLOCK_ELASTICITY_OPERATIONS<3>::SIMD_To_Unified2<float,16>(const ARRAY<RANGE<VECTOR<int,3> > >&,
                                                                         const ARRAY<VECTOR<int,3> >&,
                                                                         const ARRAY<int,VECTOR<int,3> >&,
                                                                         const ARRAY<VECTOR<int,3> >&,
                                                                         const ARRAY<int,VECTOR<int,3> >&,
                                                                         const ARRAY<VECTOR<VECTOR<float,16>,3>,VECTOR<int,3> >&,ARRAY<VECTOR<float,3>,VECTOR<int,3> >&, const SIMD_CONV_MODE );
template void BLOCK_ELASTICITY_OPERATIONS<3>::SIMD_OP_SIMD<float,16>(const ARRAY<RANGE<VECTOR<int,3> > >&,
                                                                        ARRAY<VECTOR<VECTOR<float,16>,3>,VECTOR<int,3> >&, const ARRAY<VECTOR<VECTOR<float,16>,3>,VECTOR<int,3> >&, const SIMD_OP_MODE );
template void BLOCK_ELASTICITY_OPERATIONS<3>::SIMD_UNARY_OP<float,16>(const ARRAY<RANGE<VECTOR<int,3> > >&,                                                         ARRAY<VECTOR<VECTOR<float,16>,3>,VECTOR<int,3> >&, const SIMD_UNARY_OP_MODE );
template void BLOCK_ELASTICITY_OPERATIONS<3>::SIMD_UNARY_OP<float,16>(const ARRAY<RANGE<VECTOR<int,3> > >&,ARRAY<VECTOR<VECTOR<float,16>,3>,VECTOR<int,3> >&, const ARRAY<ARRAY<PAIR<T_INDEX,int> > >&, const SIMD_UNARY_OP_MODE );
template void BLOCK_ELASTICITY_OPERATIONS<3>::SIMD_UNARY_OP<float,16>(const ARRAY<RANGE<VECTOR<int,3> > >&,ARRAY<VECTOR<VECTOR<float,16>,3>,VECTOR<int,3> >&,int, 
                                                                      const VECTOR<ARRAY<ARRAY<PAIR<T_INDEX,int> > >,4>&,const SIMD_UNARY_OP_MODE );
template void BLOCK_ELASTICITY_OPERATIONS<3>::SIMD_Prune_Duplicates<float,16>(const ARRAY<RANGE<VECTOR<int,3> > >&, ARRAY<VECTOR<VECTOR<float,16>,3>,VECTOR<int,3> >&);

template void  BLOCK_ELASTICITY_OPERATIONS<3>::SIMD_Dirty_Regions<float,16>(const ARRAY<RANGE<VECTOR<int,3> > >& subdomains, const ARRAY< RANGE<VECTOR<int,3> >, VECTOR<int,1> >& regions, const ARRAY<VECTOR<VECTOR<float,16>,3>,VECTOR<int,3> >& SIMDx);
