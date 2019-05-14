//#####################################################################
// Copyright 2002-2009, Eftychios Sifakis, Michael Doescher, Qisi Wang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REORDERING
//#####################################################################
#ifndef __REORDERING__
#define __REORDERING__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <vector>
#include "ARRAY_CSC_MATRIX.h"

namespace PhysBAM{

    template<int d>
        class REORDERING {
        using T_INDEX = VECTOR<int,d>;
        using T_INDEX_ARRAYS = ARRAY<int,T_INDEX>;  // ARRAYS_ND
        using T_COORD_ARRAYS = ARRAY<T_INDEX>;
        using T_RANGE = RANGE<T_INDEX>;

    public:
        T_INDEX_ARRAYS linear_indices;
        T_COORD_ARRAYS coordinates;

        T_INDEX_ARRAYS level_map;
        T_INDEX_ARRAYS rank_map;
        T_INDEX_ARRAYS id_map;

        ARRAY<int> level_sizes;
        ARRAY<int> aggregate_sizes;

        std::vector<T_INDEX_ARRAYS> local_id_map;
        std::vector<T_INDEX> block_dims; // 0 is subdomain

        T_INDEX size;
        T_INDEX subdomain_size;
        int depth;

        std::vector<int> cut_order; // {0,..,depth-1} in range of [1,d], sbd = depth
        // temps for development only
        bool use_queue;
        ARRAY<int> cut_axes;

        REORDERING(const T_INDEX& size_input,const T_INDEX& subdomain_size_input);

        void Create_Linear_Indices();
        void Set_Cut_Array(std::string cut_string);

        int Get_First_Index(const int level, const int rank) const;
        std::pair<int, int> Get_Index_Range(const int level, const int id) const;
        static void print(const ARRAY<int, VECTOR<int,d> >& array, int offset = 4);

        int Get_Linear_Index_Stack(const int level, const int rank, const int id) const;
        int Get_Linear_Index_Stack(const int level, const int rank, const int id, const int local_level) const;
        int Get_Linear_Index_Queue(const int level, const int rank, const int id) const;

        void print() const {print(linear_indices);}
        void print_level() const {print(level_map);}
        void print_id() const {print(id_map);}
        void print_rank() const {print(rank_map);}



        template <class T>
            void To_CSC(ARRAY_CSC_MATRIX<T>& CSC) const;

    private:
        void Create_Linear_Indices_Stack(const RANGE<VECTOR<int,d> >& coordinate_range,const int starting_index, const int level, const int rank);
        void Create_Linear_Indices_Queue();
        int fill_maps(const T_RANGE& range,const int level,const int rank);
        void fill_local_id(const T_RANGE& range,const int level);
        void Number_Entry(const RANGE<VECTOR<int,d> >& coordinate_range,const int starting_index, const int starting_id, int current_cut);
        void Split_Range(const int split_axis, RANGE<VECTOR<int,d> >(&ranges)[3]) const;
        void Get_Start_Indices(RANGE<VECTOR<int,d> >(&ranges)[3], int (&start_indices)[3]) const;

        template<int strategy>
            int Get_Split_Axis(const VECTOR<int,d>& range_size);
        template<int strategy>
            int Get_Split_Axis(int &current_cut);

        template<class T>
            class Circular_Queue;
    };
}

#endif
