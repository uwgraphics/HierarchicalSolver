//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PCG_SPARSE_MPI
//#####################################################################
#ifndef __PCG_SPARSE_MPI__
#define __PCG_SPARSE_MPI__

#ifdef USE_MPI

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
namespace PhysBAM{

class SPARSE_MATRIX_PARTITION;

template<class T_GRID>
class PCG_SPARSE_MPI:public NONCOPYABLE
{
public:
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    PCG_SPARSE<T>& pcg;
    MPI::Intracomm& comm;
#ifdef USE_PTHREADS
    ARRAY<THREAD_PACKAGE>* buffers;
    pthread_mutex_t* lock;
    pthread_barrier_t* barr;
#endif
    THREADED_UNIFORM_GRID<T_GRID>* thread_grid;
    MPI_THREADED_UNIFORM_GRID<T_GRID>* mpi_threaded_grid;
    SPARSE_MATRIX_PARTITION& partition;
    ARRAY<MPI::Datatype> boundary_datatypes,ghost_datatypes;
    ARRAY<ARRAY<int> > columns_to_send;
    ARRAY<ARRAY<int> > columns_to_receive;

    PCG_SPARSE_MPI(PCG_SPARSE<T>& pcg_input,MPI::Intracomm& comm_input,SPARSE_MATRIX_PARTITION& partition_input)
        :pcg(pcg_input),comm(comm_input),thread_grid(0),mpi_threaded_grid(0),partition(partition_input)
    {}

#ifdef USE_PTHREADS
    PCG_SPARSE_MPI(PCG_SPARSE<T>& pcg_input,MPI::Intracomm& comm_input,ARRAY<THREAD_PACKAGE>* buffers_input,pthread_mutex_t* lock_input,pthread_barrier_t* barr_input,SPARSE_MATRIX_PARTITION& partition_input)
        :pcg(pcg_input),comm(comm_input),buffers(buffers_input),lock(lock_input),barr(barr_input),thread_grid(0),mpi_threaded_grid(0),partition(partition_input)
    {}
#endif

    virtual ~PCG_SPARSE_MPI()
    {MPI_UTILITIES::Free_Elements_And_Clean_Memory(boundary_datatypes);MPI_UTILITIES::Free_Elements_And_Clean_Memory(ghost_datatypes);}

    template<class TYPE> TYPE Global_Sum(const TYPE& input)
    {TYPE output;MPI_UTILITIES::Reduce(input,output,MPI::SUM,comm);return output;}

#ifdef USE_PTHREADS
    template<class TYPE> TYPE Global_Sum_Threaded(const TYPE& input)
    {TYPE output=TYPE(0);
    THREAD_PACKAGE pack(sizeof(TYPE));*(TYPE*)(&pack.buffer(1))=input;
    pthread_mutex_lock(lock);
    buffers->Append(pack);
    pthread_mutex_unlock(lock);
    pthread_barrier_wait(barr);
    for(int i=1;i<=buffers->m;i++) output+=*(TYPE*)(&((*buffers)(i).buffer(1)));
    pthread_barrier_wait(barr);
    pthread_mutex_lock(lock);
    buffers->m=0;
    pthread_mutex_unlock(lock);
    pthread_barrier_wait(barr);
    return output;}

    template<class TYPE> TYPE Global_Sum_MPI_Threaded(const TYPE& input)
    {TYPE output;TYPE input_sum_all_threads=TYPE(0);
    THREAD_PACKAGE pack(sizeof(TYPE));*(TYPE*)(&pack.buffer(1))=input;
    pthread_mutex_lock(mpi_threaded_grid->lock);
    mpi_threaded_grid->buffers->Append(pack);
    pthread_mutex_unlock(mpi_threaded_grid->lock);
    pthread_barrier_wait(mpi_threaded_grid->barr);
    if(mpi_threaded_grid->tid==1){ 
        for(int i=1;i<=mpi_threaded_grid->buffers->m;i++) input_sum_all_threads+=*(TYPE*)(&(*mpi_threaded_grid->buffers)(i).buffer(1));
        MPI_UTILITIES::Reduce(input_sum_all_threads,output,MPI::SUM,comm);
        THREAD_PACKAGE output_pack(sizeof(TYPE));*(TYPE*)(&output_pack.buffer(1))=output;
        mpi_threaded_grid->buffers->Append(output_pack);}
    pthread_barrier_wait(mpi_threaded_grid->barr);
    output=*(TYPE*)(&(*mpi_threaded_grid->buffers)(mpi_threaded_grid->buffers->m).buffer(1));
    pthread_barrier_wait(mpi_threaded_grid->barr);
    if(mpi_threaded_grid->tid==1) mpi_threaded_grid->buffers->m=0;
    pthread_barrier_wait(mpi_threaded_grid->barr);
    return output;}
#else
    template<class TYPE> TYPE Global_Sum_Threaded(const TYPE& input)
    {PHYSBAM_NOT_IMPLEMENTED();}

    template<class TYPE> TYPE Global_Sum_MPI_Threaded(const TYPE& input)
    {PHYSBAM_NOT_IMPLEMENTED();}
#endif

    template<class TYPE> TYPE Global_Max(const TYPE& input)
    {TYPE output;MPI_UTILITIES::Reduce(input,output,MPI::MAX,comm);return output;}

#ifdef USE_PTHREADS
    template<class TYPE> TYPE Global_Max_Threaded(const TYPE& input)
    {TYPE output=input;
    THREAD_PACKAGE pack(sizeof(TYPE));*(TYPE*)(&pack.buffer(1))=input;
    pthread_mutex_lock(lock);
    buffers->Append(pack);
    pthread_mutex_unlock(lock);
    pthread_barrier_wait(barr);
    for(int i=1;i<=buffers->m;i++) output=max(*(TYPE*)(&((*buffers)(i).buffer(1))),output);
    pthread_barrier_wait(barr);
    pthread_mutex_lock(lock);
    buffers->m=0;
    pthread_mutex_unlock(lock);
    pthread_barrier_wait(barr);
    return output;}

    template<class TYPE> TYPE Global_Max_MPI_Threaded(const TYPE& input)
    {TYPE output;TYPE input_max_all_threads=input;
    THREAD_PACKAGE pack(sizeof(TYPE));*(TYPE*)(&pack.buffer(1))=input;
    pthread_mutex_lock(mpi_threaded_grid->lock);
    mpi_threaded_grid->buffers->Append(pack);
    pthread_mutex_unlock(mpi_threaded_grid->lock);
    pthread_barrier_wait(mpi_threaded_grid->barr);
    if(mpi_threaded_grid->tid==1){ 
        for(int i=1;i<=mpi_threaded_grid->buffers->m;i++) if(input_max_all_threads<*(TYPE*)(&(*mpi_threaded_grid->buffers)(i).buffer(1)))
            input_max_all_threads=*(TYPE*)(&(*mpi_threaded_grid->buffers)(i).buffer(1));
        MPI_UTILITIES::Reduce(input_max_all_threads,output,MPI::MAX,comm);
        THREAD_PACKAGE output_pack(sizeof(TYPE));*(TYPE*)(&output_pack.buffer(1))=output;
        mpi_threaded_grid->buffers->Append(output_pack);}
    pthread_barrier_wait(mpi_threaded_grid->barr);
    output=*(TYPE*)(&(*mpi_threaded_grid->buffers)(mpi_threaded_grid->buffers->m).buffer(1));
    pthread_barrier_wait(mpi_threaded_grid->barr);
    if(mpi_threaded_grid->tid==1) mpi_threaded_grid->buffers->m=0;
    pthread_barrier_wait(mpi_threaded_grid->barr);
    return output;}
#else
    template<class TYPE> TYPE Global_Max_Threaded(const TYPE& input)
    {PHYSBAM_NOT_IMPLEMENTED();}

    template<class TYPE> TYPE Global_Max_MPI_Threaded(const TYPE& input)
    {PHYSBAM_NOT_IMPLEMENTED();}
#endif

    virtual void Fill_Ghost_Cells(VECTOR_ND<T>& v)
    {ARRAY<MPI::Request> requests;requests.Preallocate(2*partition.number_of_sides);
    for(int s=1;s<=partition.number_of_sides;s++)if(boundary_datatypes(s)!=MPI::DATATYPE_NULL) requests.Append(comm.Isend(v.x-1,1,boundary_datatypes(s),partition.neighbor_ranks(s),s));
    for(int s=1;s<=partition.number_of_sides;s++)if(ghost_datatypes(s)!=MPI::DATATYPE_NULL) requests.Append(comm.Irecv(v.x-1,1,ghost_datatypes(s),partition.neighbor_ranks(s),((s-1)^1)+1));
    MPI_UTILITIES::Wait_All(requests);}
    
//#####################################################################
    void Serial_Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,VECTOR_ND<T>& q,VECTOR_ND<T>& s,VECTOR_ND<T>& r,VECTOR_ND<T>& k,VECTOR_ND<T>& z,const int tag,const T tolerance=1e-7);
    void Parallel_Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance=1e-7,const bool recompute_preconditioner=true);
    void Parallel_Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x_local,VECTOR_ND<T>& b_local,const ARRAY<VECTOR<int,2> >& proc_column_index_boundaries,
        const T tolerance=1e-7,const bool recompute_preconditioner=true);
    void Parallel_Solve_Threaded(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance=1e-7,const bool recompute_preconditioner=true);
    void Parallel_Solve_MPI_Threaded(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance=1e-7,const bool recompute_preconditioner=true);
    void Find_Ghost_Regions(SPARSE_MATRIX_FLAT_NXN<T>& A,const ARRAY<VECTOR<int,2> >& proc_column_index_boundaries);
    void Find_Ghost_Regions_Threaded(SPARSE_MATRIX_FLAT_NXN<T>& A,const ARRAY<VECTOR<int,2> >& proc_column_index_boundaries);
    void Fill_Ghost_Cells_Far(VECTOR_ND<T>& x);
    void Fill_Ghost_Cells_Threaded(VECTOR_ND<T>& x); 
    void Fill_Ghost_Cells_MPI_Threaded(VECTOR_ND<T>& v);
    int Get_Send_Tag(int s);
    int Get_Recv_Tag(int s);
    void Fill_Ghost_Cells_Between_Threads(VECTOR_ND<T>& v);
    virtual void Initialize_Datatypes();
//#####################################################################
};
}
#endif
#endif
