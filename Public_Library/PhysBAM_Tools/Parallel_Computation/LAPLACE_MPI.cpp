//#####################################################################
// Copyright 2005-2010, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE_PDE_Linear/LAPLACE_RLE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_RLE_GRID.h>
#endif
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/LAPLACE_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/LAPLACE_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#endif
#include <PhysBAM_Tools/Parallel_Computation/FLOOD_FILL_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_MPI_THREADED.h>
#endif
#ifdef USE_PTHREADS
#include <PhysBAM_Tools/Parallel_Computation/THREAD_PACKAGE.h>
#endif
using namespace PhysBAM;

#ifdef USE_MPI

//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> LAPLACE_MPI<T_GRID>::
LAPLACE_MPI(T_LAPLACE& laplace)
    :mpi_grid(laplace.mpi_grid),local_grid(laplace.grid),local_pcg(laplace.pcg),number_of_regions(laplace.number_of_regions),number_of_global_regions(-1),filled_region_colors(laplace.filled_region_colors),
    filled_region_touches_dirichlet(laplace.filled_region_touches_dirichlet),solve_neumann_regions(laplace.solve_neumann_regions),psi_N(laplace.psi_N)
{
    groups=new ARRAY<MPI::Group>;
    communicators=new ARRAY<MPI::Intracomm>;
#ifdef USE_PTHREADS
    buffers=new ARRAY<ARRAY<THREAD_PACKAGE>*>;
    locks=new ARRAY<pthread_mutex_t*>;
    barriers=new ARRAY<pthread_barrier_t*>;
#endif
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> LAPLACE_MPI<T_GRID>::
~LAPLACE_MPI()
{
    MPI_UTILITIES::Free_Elements_And_Clean_Memory(*groups);delete groups;
    MPI_UTILITIES::Free_Elements_And_Clean_Memory(*communicators);delete communicators;
#ifdef USE_PTHREADS
    delete buffers;
    delete locks;
    delete barriers;
#endif
}
//#####################################################################
// Function Synchronize_Solution_Regions
//#####################################################################
template<class T_GRID> void LAPLACE_MPI<T_GRID>::
Synchronize_Solution_Regions()
{
    FLOOD_FILL_MPI<T_GRID> flood_fill(*mpi_grid,local_grid,psi_N,number_of_regions,filled_region_colors,filled_region_ranks,&filled_region_touches_dirichlet);
    number_of_global_regions=flood_fill.Synchronize_Colors();
    global_colors=flood_fill.global_colors;
    max_global_color=flood_fill.max_global_color;
    
    if(!(mpi_grid->mpi_threaded_grid || mpi_grid->threaded_grid)){
        // allocate communicators for each color
        ARRAY<MPI::Group> new_groups(filled_region_ranks.m);
        ARRAY<MPI::Intracomm> new_communicators(filled_region_ranks.m);
        for(int color=1;color<=filled_region_ranks.m;color++){
            new_groups(color)=mpi_grid->group->Incl(filled_region_ranks(color).m,&filled_region_ranks(color)(1));
            int i;
            for(i=1;i<=groups->m;i++)if((*groups)(i)!=MPI::GROUP_NULL && MPI::Group::Compare((*groups)(i),new_groups(color))==MPI::IDENT){
                new_communicators(color)=(*communicators)(i);(*communicators)(i)=MPI::COMM_NULL;(*groups)(i).Free();break;}
            if(i>groups->m) new_communicators(color)=mpi_grid->comm->Create(new_groups(color));} // NOTE: returns null communicator for processes not in the group!
        MPI_UTILITIES::Free_Elements_And_Clean_Memory(*groups);MPI_UTILITIES::Free_Elements_And_Clean_Memory(*communicators);
        groups->Exchange(new_groups);communicators->Exchange(new_communicators);
        // allocate partitions and compute neighbor ranks
        partitions.Resize(filled_region_ranks.m);
        for(int color=1;color<=filled_region_ranks.m;color++){
            partitions(color).Set_Number_Of_Sides(T_PARALLEL_GRID::number_of_faces_per_cell);
            for(int s=1;s<=T_PARALLEL_GRID::number_of_faces_per_cell;s++){
                int global_rank=mpi_grid->side_neighbor_ranks(s);
                if(global_rank==MPI::PROC_NULL) partitions(color).neighbor_ranks(s)=MPI::PROC_NULL;
                else MPI::Group::Translate_ranks(*mpi_grid->group,1,&global_rank,(*groups)(color),&partitions(color).neighbor_ranks(s));}}}
#ifdef USE_PTHREADS
    else if(!mpi_grid->mpi_threaded_grid && mpi_grid->threaded_grid){
        ARRAY<ARRAY<THREAD_PACKAGE>*> new_buffers(filled_region_ranks.m);
        ARRAY<pthread_mutex_t*> new_locks(filled_region_ranks.m);
        ARRAY<pthread_barrier_t*> new_barriers(filled_region_ranks.m);
        int my_rank=mpi_grid->threaded_grid->rank;
        if(my_rank==0){
            for(int color=1;color<=filled_region_ranks.m;color++){
                new_buffers(color)=new ARRAY<THREAD_PACKAGE>;
                new_locks(color)=new pthread_mutex_t;
                pthread_mutex_init(new_locks(color),NULL);
                new_barriers(color)=new pthread_barrier_t;
                pthread_barrier_init(new_barriers(color),NULL,filled_region_ranks(color).m);}
            THREAD_PACKAGE pack(sizeof(pthread_mutex_t*)*filled_region_ranks.m+sizeof(pthread_barrier_t*)*filled_region_ranks.m+sizeof(ARRAY<THREAD_PACKAGE>*)*filled_region_ranks.m);
            int position=0;
            for(int color=1;color<=filled_region_ranks.m;color++){
                *(ARRAY<THREAD_PACKAGE>**)(&pack.buffer(position+1))=new_buffers(color);position+=sizeof(ARRAY<THREAD_PACKAGE>*);
                *(pthread_mutex_t**)(&pack.buffer(position+1))=new_locks(color);position+=sizeof(pthread_mutex_t*);
                *(pthread_barrier_t**)(&pack.buffer(position+1))=new_barriers(color);position+=sizeof(pthread_barrier_t*);}
            pthread_mutex_lock(mpi_grid->threaded_grid->lock);
            mpi_grid->threaded_grid->buffers.Append(pack);
            pthread_mutex_unlock(mpi_grid->threaded_grid->lock);
            pthread_barrier_wait(mpi_grid->threaded_grid->barr);
            pthread_barrier_wait(mpi_grid->threaded_grid->barr);
            mpi_grid->threaded_grid->buffers.m=0;
            pthread_barrier_wait(mpi_grid->threaded_grid->barr);}
        else{
            pthread_barrier_wait(mpi_grid->threaded_grid->barr);
            THREAD_PACKAGE& pack=mpi_grid->threaded_grid->buffers(1);int position=0;
            for(int color=1;color<=filled_region_ranks.m;color++){
                new_buffers(color)=*(ARRAY<THREAD_PACKAGE>**)(&pack.buffer(position+1));position+=sizeof(ARRAY<THREAD_PACKAGE>*);
                new_locks(color)=*(pthread_mutex_t**)(&pack.buffer(position+1));position+=sizeof(pthread_mutex_t*);
                new_barriers(color)=*(pthread_barrier_t**)(&pack.buffer(position+1));position+=sizeof(pthread_barrier_t*);}
            pthread_barrier_wait(mpi_grid->threaded_grid->barr);
            pthread_barrier_wait(mpi_grid->threaded_grid->barr);}
        buffers->Exchange(new_buffers);locks->Exchange(new_locks);barriers->Exchange(new_barriers);
        partitions.Resize(filled_region_ranks.m);
        for(int color=1;color<=filled_region_ranks.m;color++){
            bool region_contain_my_rank=filled_region_ranks(color).Contains(my_rank);
            partitions(color).Set_Number_Of_Sides(T_PARALLEL_GRID::number_of_faces_per_cell);
            for(int s=1;s<=T_PARALLEL_GRID::number_of_faces_per_cell;s++){
                int neighbor_rank=mpi_grid->threaded_grid->side_neighbor_ranks(s);
                if(neighbor_rank==-1 || !region_contain_my_rank || !filled_region_ranks(color).Contains(neighbor_rank)) partitions(color).neighbor_ranks(s)=-1;
                else partitions(color).neighbor_ranks(s)=neighbor_rank;}}}
    else{
        int my_rank=mpi_grid->mpi_threaded_grid->rank;
        partitions.Resize(filled_region_ranks.m);
        for(int color=1;color<=filled_region_ranks.m;color++){
            bool region_contain_my_rank=filled_region_ranks(color).Contains(my_rank);
            partitions(color).Set_Number_Of_Sides(T_PARALLEL_GRID::number_of_faces_per_cell);
            for(int s=1;s<=T_PARALLEL_GRID::number_of_faces_per_cell;s++){
                int neighbor_rank=mpi_grid->mpi_threaded_grid->side_neighbor_ranks(s);
                if(neighbor_rank==-1 || !region_contain_my_rank || !filled_region_ranks(color).Contains(neighbor_rank)) partitions(color).neighbor_ranks(s)=-1;
                else partitions(color).neighbor_ranks(s)=neighbor_rank;}}}
#endif
}
//#####################################################################
// Function Update_Solution_Regions_For_Solid_Fluid_Coupling
//#####################################################################
template<class T_GRID> void LAPLACE_MPI<T_GRID>::
Update_Solution_Regions_For_Solid_Fluid_Coupling(const T_MPI_GRID& mpi_grid)
{
    MPI_UTILITIES::Free_Elements_And_Clean_Memory(*groups);MPI_UTILITIES::Free_Elements_And_Clean_Memory(*communicators);
    groups->Resize(1);communicators->Resize(1);
    (*communicators)(1)=mpi_grid.comm->Dup();(*groups)(1)=*new MPI::Group((*communicators)(1).Get_group());
}
//#####################################################################
// Function Update_Solution_Regions_For_Solid_Fluid_Coupling
//#####################################################################
template<class T_GRID> int LAPLACE_MPI<T_GRID>::
Get_Total_Number_Of_Threads(const int input,const int color)
{
    int output;
    MPI_UTILITIES::Reduce(input,output,MPI::SUM,(*communicators)(color));
    return output;
}
//#####################################################################
// Function Solve
//#####################################################################
template<class T_GRID> void LAPLACE_MPI<T_GRID>::
Solve_Threaded(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,ARRAY<INTERVAL<int> >& interior_indices,ARRAY<ARRAY<INTERVAL<int> > >& ghost_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,VECTOR_ND<T>& q,VECTOR_ND<T>& s,VECTOR_ND<T>& r,VECTOR_ND<T>& k,VECTOR_ND<T>& z,const T tolerance,const int color,const int multi_proc_mode)
{
    if(color>filled_region_ranks.m){
        if(multi_proc_mode){local_pcg_threaded->Solve(domain,domain_index,interior_indices,ghost_indices,A,x,b,tolerance);return;}
        else{local_pcg.Solve(A,x,b,q,s,r,k,z,tolerance);return;}}
    else{
        PCG_SPARSE_MPI_THREADED<TV> pcg_mpi(*local_pcg_threaded,(*communicators)(color),partitions(color));
        assert(Use_Parallel_Solve());
        pcg_mpi.Solve(domain,domain_index,interior_indices,ghost_indices,A,x,b,tolerance);}
}
//#####################################################################
// Function Solve
//#####################################################################
template<class T_GRID> void LAPLACE_MPI<T_GRID>::
Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,VECTOR_ND<T>& q,VECTOR_ND<T>& s,VECTOR_ND<T>& r,VECTOR_ND<T>& k,VECTOR_ND<T>& z,const T tolerance,const int color)
{
    if(color>filled_region_ranks.m){local_pcg.Solve(A,x,b,q,s,r,k,z,tolerance);return;}
    else{
        MPI::Intracomm* comm;
        if(mpi_grid->mpi_threaded_grid) comm=mpi_grid->mpi_threaded_grid->mpi_grid->comm;
        else if(mpi_grid->threaded_grid) comm=0;
        else comm=&(*communicators)(color);
#ifdef USE_PTHREADS
        if(mpi_grid->mpi_threaded_grid) pthread_mutex_lock(mpi_grid->mpi_threaded_grid->lock);
        MPI::Intracomm threaded_comm=MPI::COMM_NULL;
        ARRAY<THREAD_PACKAGE>* threaded_buffer=mpi_grid->threaded_grid?(*buffers)(color):0;
        pthread_mutex_t* threaded_lock=mpi_grid->threaded_grid?(*locks)(color):0;
        pthread_barrier_t* threaded_barrier=mpi_grid->threaded_grid?(*barriers)(color):0;
        PCG_SPARSE_MPI<T_GRID> pcg_mpi(local_pcg,mpi_grid->threaded_grid?threaded_comm:*comm,threaded_buffer,threaded_lock,threaded_barrier,partitions(color));
        if(mpi_grid->mpi_threaded_grid) pthread_mutex_unlock(mpi_grid->mpi_threaded_grid->lock);
#else
        PCG_SPARSE_MPI<T_GRID> pcg_mpi(local_pcg,*comm,partitions(color));
#endif
        pcg_mpi.thread_grid=mpi_grid->threaded_grid;
        pcg_mpi.mpi_threaded_grid=mpi_grid->mpi_threaded_grid;
        if(Use_Parallel_Solve()) pcg_mpi.Parallel_Solve(A,x,b,tolerance);
        else pcg_mpi.Serial_Solve(A,x,b,q,s,r,k,z,1234,tolerance);}
}
//#####################################################################
// Function Solve
//#####################################################################
template<class T_GRID> void LAPLACE_MPI<T_GRID>::
Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance,const int color,const ARRAY<VECTOR<int,2> >& global_column_index_boundaries)
{
    SPARSE_MATRIX_PARTITION temp_partition;
    PCG_SPARSE_MPI<T_GRID> pcg_mpi(local_pcg,(*communicators)(color),temp_partition);
    pcg_mpi.thread_grid=mpi_grid->threaded_grid;
    pcg_mpi.Parallel_Solve(A,x,b,global_column_index_boundaries,tolerance,true);
}
//#####################################################################

#else

//#####################################################################
template<class T_GRID> LAPLACE_MPI<T_GRID>::LAPLACE_MPI(T_LAPLACE& laplace):mpi_grid(laplace.mpi_grid),local_grid(laplace.grid),local_pcg(laplace.pcg),
    number_of_regions(laplace.number_of_regions),filled_region_colors(laplace.filled_region_colors),filled_region_touches_dirichlet(laplace.filled_region_touches_dirichlet),
    solve_neumann_regions(laplace.solve_neumann_regions),psi_N(laplace.psi_N),groups(0),communicators(0){}
template<class T_GRID> LAPLACE_MPI<T_GRID>::~LAPLACE_MPI(){}
template<class T_GRID> void LAPLACE_MPI<T_GRID>::Synchronize_Solution_Regions(){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T_GRID> int LAPLACE_MPI<T_GRID>::Get_Total_Number_Of_Threads(const int input,const int color){PHYSBAM_FUNCTION_IS_NOT_DEFINED();return 1;}
template<class T_GRID> void LAPLACE_MPI<T_GRID>::Update_Solution_Regions_For_Solid_Fluid_Coupling(const T_MPI_GRID& mpi_grid){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T_GRID> void LAPLACE_MPI<T_GRID>::Solve(SPARSE_MATRIX_FLAT_NXN<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,
    const T,const int){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T_GRID> void LAPLACE_MPI<T_GRID>::Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance,const int color,
    const ARRAY<VECTOR<int,2> >& global_column_index_boundaries){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T_GRID> void LAPLACE_MPI<T_GRID>::Solve_Threaded(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,ARRAY<INTERVAL<int> >& interior_indices,
    ARRAY<ARRAY<INTERVAL<int> > >& ghost_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,VECTOR_ND<T>& q,VECTOR_ND<T>& s,VECTOR_ND<T>& r,VECTOR_ND<T>& k,
    VECTOR_ND<T>& z,const T tolerance,const int color,const int multi_proc_mode)
{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################

#endif

template<class T_GRID> bool& LAPLACE_MPI<T_GRID>::
Use_Parallel_Solve()
{
    static bool use_parallel_solve=true;
    return use_parallel_solve;
}

//#####################################################################
template class LAPLACE_MPI<GRID<VECTOR<float,1> > >;
template class LAPLACE_MPI<GRID<VECTOR<float,2> > >;
template class LAPLACE_MPI<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LAPLACE_MPI<GRID<VECTOR<double,1> > >;
template class LAPLACE_MPI<GRID<VECTOR<double,2> > >;
template class LAPLACE_MPI<GRID<VECTOR<double,3> > >;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template class LAPLACE_MPI<RLE_GRID_2D<float> >;
template class LAPLACE_MPI<RLE_GRID_3D<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LAPLACE_MPI<RLE_GRID_2D<double> >;
template class LAPLACE_MPI<RLE_GRID_3D<double> >;
#endif
#endif
