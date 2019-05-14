//#####################################################################
// Copyright 2002-2009, Eftychios Sifakis, Michael Doescher, Qisi Wang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REORDERING
//#####################################################################

#include "./REORDERING.h"

#include <utility>
#include <cstddef>
#include <iomanip>

#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include "ARBITRARY_RANGE_ITERATOR.h"

namespace PhysBAM{

    template<int d>
    REORDERING<d>::REORDERING(const T_INDEX& size_input,const T_INDEX& subdomain_size_input):
        size(size_input),
            subdomain_size(subdomain_size_input),
            linear_indices(T_RANGE(T_INDEX::All_Ones_Vector(),size_input)),
            coordinates(size_input.Product()),
            level_map(T_RANGE(T_INDEX::All_Ones_Vector(),size_input)),
            rank_map(T_RANGE(T_INDEX::All_Ones_Vector(),size_input)),
            id_map(T_RANGE(T_INDEX::All_Ones_Vector(),size_input)),
            depth(0),
            use_queue(false)
        {
            for (int v=1; v<=d; v++)
            {
                int ds=size(v),sds=subdomain_size(v);
                while (sds%2==1)
                {
                    ds=ds>>1;
                    sds=sds>>1;
                }
                while(ds%2==1)
                {
                    depth++;
                    ds=ds>>1;
                }
            }

            level_sizes.Resize(depth+1);
            aggregate_sizes.Resize(depth+1);
            cut_order.resize(depth);

            local_id_map.resize(depth+1);
            block_dims.resize(depth + 1);
        }

    template<int d>
    inline void REORDERING<d>::Split_Range(const int split_axis, RANGE<VECTOR<int,d> >(&ranges)[3]) const
        {
            PHYSBAM_ASSERT(ranges[0] == ranges[1]);
            PHYSBAM_ASSERT(ranges[1] == ranges[2]);

            const RANGE<VECTOR<int,d> > coordinate_range = ranges[0];

            int bias =  (coordinate_range.min_corner(split_axis)>coordinate_range.max_corner(split_axis))?-1:1;
            int midpoint;
            if ((coordinate_range.max_corner(split_axis)-coordinate_range.min_corner(split_axis))%2)
                midpoint=(coordinate_range.min_corner(split_axis)+coordinate_range.max_corner(split_axis)-bias)/2;
            else
                midpoint=(coordinate_range.min_corner(split_axis)+coordinate_range.max_corner(split_axis))/2;

            ranges[0].max_corner(split_axis)=midpoint-bias;
            //ranges[0].max_corner(split_axis)=coordinate_range.min_corner(split_axis);
            ranges[1].max_corner(split_axis)=midpoint+bias;
            ranges[1].min_corner(split_axis)=coordinate_range.max_corner(split_axis);
            ranges[2].min_corner(split_axis)=ranges[2].max_corner(split_axis)=midpoint;
        }

    template<int d>
    inline void REORDERING<d>::Get_Start_Indices(RANGE<VECTOR<int,d> >(&ranges)[3], int (&start_indices)[3]) const
        {

            PHYSBAM_ASSERT(start_indices[0] == start_indices[1]);
            PHYSBAM_ASSERT(start_indices[1] == start_indices[2]);
            const int starting_index = start_indices[0];
            start_indices[0]=starting_index;
            start_indices[1]=start_indices[0]+(abs(ranges[0].Edge_Lengths())+1).Product();
            start_indices[2]=start_indices[1]+(abs(ranges[1].Edge_Lengths())+1).Product();
        }

    template<int d>
    template<int strategy>
    inline int REORDERING<d>::Get_Split_Axis(const VECTOR<int,d>& range_size)
        {
            int split_axis = d;
            switch (strategy)
            {
            case 1:
                // TODO : Try different split strategies (maybe attempting to eliminate overlap between last 2 split interfaces)
                // Strategy #1 corrected : Split along the *first* of the largest axes.
                // Really need to compare all to get the correct largest dimension
                for (int v=1; v<d; v++)
                    if(range_size(v)>range_size(split_axis)) split_axis=v;
                break;
            case 2:
                // Need to find axis to split
                // Strategy #1 : Split along the *first* of the largest axes.
                // This strategy not working with interface on second dimension
                for(;split_axis>1;split_axis--)
                    if(range_size(split_axis-1)<range_size(split_axis)) break;
                break;
                // TODO : Work with making the coordinate ranges "directional", i.e. incorporating mirroring when necessary (be careful, since this breaks conventions)
            default:
                // Undefined strategy
                PHYSBAM_ASSERT(false);
            }
            return split_axis;
        }


    template<int d>
    template<int strategy>
        inline int REORDERING<d>::Get_Split_Axis(int &current_cut)
        {

            int split_axis = d;
            const int number_of_cuts = cut_axes.Size();
            switch (strategy)
            {
            case 3:
                split_axis=cut_axes(current_cut++);
                if (current_cut==number_of_cuts+1) current_cut = 1;
                break;
            default:
                // Undefined strategy
                PHYSBAM_ASSERT(false);
            }
            return split_axis;
        }

    template<int d>
    inline void REORDERING<d>::Number_Entry(const RANGE<VECTOR<int,d> >& coordinate_range,const int starting_index, const int starting_id, int current_cut)
        {

            const bool use_cut = (current_cut!=0);
            const T_INDEX range_size=abs(coordinate_range.Edge_Lengths())+1;

            int split_axis;

            if (use_cut)
                split_axis = Get_Split_Axis<3>(current_cut);
            else
                split_axis = Get_Split_Axis<1>(range_size);
            //LOG::cout<<"Number_Entry: "<<split_axis<<std::endl;
            if (range_size(split_axis) <= 2)
            {
                int index = starting_index;
                int id = starting_id;
                for (ARBITRARY_RANGE_ITERATOR<d> iterator(coordinate_range); iterator.Valid(); iterator.Next())
                    {
                        linear_indices(iterator.Index()) = index++;
                        id_map(iterator.Index()) = id++;
                    }
                return;
            }
#if 0
            if(coordinate_range.Edge_Lengths()==T_INDEX()) {
                // LOG::cout<<"Setting linear_indices("<<coordinate_range.min_corner<<") <- "<<starting_index<<std::endl;
                linear_indices(coordinate_range.min_corner)=starting_index;
                return;
            }
#endif
            T_RANGE ranges[3] = {coordinate_range,coordinate_range,coordinate_range};
            int start_indices[3] = {starting_index,starting_index,starting_index};
            int start_ids[3] = {starting_id,starting_id,starting_id};

            Split_Range(split_axis, ranges);
            Get_Start_Indices(ranges, start_indices);
            Get_Start_Indices(ranges, start_ids);

            Number_Entry(ranges[0],start_indices[0], start_ids[0], current_cut);
            Number_Entry(ranges[1],start_indices[1], start_ids[1], current_cut);
            bool sbd = true;
            for (int v=1; v<=d; v++)
              if (coordinate_range.min_corner(v)==coordinate_range.max_corner(v)) {
                sbd=false;
                break;
              }
            if (sbd)
              Number_Entry(ranges[2],start_indices[2], start_ids[2], current_cut);
            else
              Number_Entry(ranges[2],start_indices[2], start_ids[2], 0);
        }

    template<int d>
    void REORDERING<d>::Create_Linear_Indices()
        {
            if (use_queue)
                Create_Linear_Indices_Queue();
            else
                Create_Linear_Indices_Stack(T_RANGE(T_INDEX(1),size), 1, 0, 0);

            // reverse mapping
            for (RANGE_ITERATOR<d> iterator(T_INDEX(1), size); iterator.Valid(); iterator.Next())
            {
                const T_INDEX& coord = iterator.Index();
                const size_t index = linear_indices(coord);
                coordinates(index) = coord;
            }

            // normalize sizes
            for (int l=0; l<=depth; l++)
            {
                level_sizes(l+1) /= (1<<l);
            }

            // compute aggregated_sizes
            aggregate_sizes(depth+1) = level_sizes(depth+1);
            for (int l = depth; l>=1; l--)
                aggregate_sizes(l) = 2*aggregate_sizes(l+1)+level_sizes(l);

            // compute block dimensions
            block_dims[0] = subdomain_size;
            for (int l=1; l<=depth; l++)
            {
                block_dims[l] = block_dims[l-1];
                const int axis = cut_order[depth-l];
                // LOG::cout<<"level = "<<l<<std::endl;
                // LOG::cout<<axis<<std::endl;
                block_dims[l](axis) *= 2;
                block_dims[l](axis) += 1;
                //LOG::cout<<block_dims[l]<<std::endl;
            }
            //LOG::cout<<level_sizes<<std::endl;
            //LOG::cout<<aggregate_sizes<<std::endl;
        }

    template<int d>
    inline void REORDERING<d>::Create_Linear_Indices_Stack(const RANGE<VECTOR<int,d> >& coordinate_range,const int starting_index, const int level, const int rank)
        {
            const bool use_cut = (cut_axes.Size() != 0);
            const T_INDEX range_size=abs(coordinate_range.Edge_Lengths())+1;

            int split_axis;
            int current_cut = 0;

            if (use_cut)
            {
                current_cut = level%cut_axes.Size()+1;
                split_axis = Get_Split_Axis<3>(current_cut);
            }
            else
                split_axis = Get_Split_Axis<1>(range_size);

            if(level==depth) {
                // number the subdomain
                PHYSBAM_ASSERT(range_size==subdomain_size);
                //LOG::cout<<range_size<<subdomain_size<<std::endl;
                fill_maps(coordinate_range, level, rank);

                Number_Entry(coordinate_range,starting_index,0,current_cut);
                if (rank == 0)
                    fill_local_id(coordinate_range,level);
            }
            else
            {
                T_RANGE ranges[3] = {coordinate_range,coordinate_range,coordinate_range};
                int start_indices[3] = {starting_index,starting_index,starting_index};



                if (rank==0){
                    cut_order[level] = split_axis;
                }

                Split_Range(split_axis, ranges);
                Get_Start_Indices(ranges, start_indices);

                Create_Linear_Indices_Stack(ranges[0],start_indices[0], level+1, rank<<1);
                Create_Linear_Indices_Stack(ranges[1],start_indices[1], level+1, (rank<<1)+1);
                fill_maps(ranges[2],level,rank);
                Number_Entry(ranges[2],start_indices[2], 0, 0);
                if (rank==0)
                    fill_local_id(ranges[2],level);
            }
            // if (level == 0 && level != depth){
            //     for (int i=0; i<depth; i++)
            //         LOG::cout<<(char)(cut_order[i]+'w');
            //     LOG::cout<<std::endl;}
        }



    template<int d>
    inline void REORDERING<d>::Create_Linear_Indices_Queue()
        {
            int index = size.Product();
            int rank = 1;

            const bool use_cut = (cut_axes.Size() != 0);

            using T_QENTRY = std::pair<T_RANGE, int>;
            Circular_Queue<T_QENTRY> queue(1<<depth);
            queue.enqueue(T_QENTRY(T_RANGE(T_INDEX(1), size), 0));

            while (index)
            {
                const auto current = queue.dequeue();

                const T_RANGE coordinate_range = current.first;
                const int level = current.second;
                const T_INDEX range_size=abs(coordinate_range.Edge_Lengths())+1;

                if (rank == 0)
                {
                    rank = 1<<level;
                }

                int current_cut = 0;
                if (use_cut)
                    current_cut = level%cut_axes.Size()+1;

                if(range_size==subdomain_size) {
                    // number the subdomain
                    PHYSBAM_ASSERT(level==depth);
                    index -= fill_maps(coordinate_range, level, --rank);
                    Number_Entry(coordinate_range,index+1,0,current_cut);
                    if (rank==0)
                        fill_local_id(coordinate_range,level);
                }
                else
                {
                    int split_axis;
                    if (use_cut)
                        split_axis=Get_Split_Axis<3>(current_cut);
                    else
                        split_axis = Get_Split_Axis<1>(range_size);

                    if (rank-1==0){
                        cut_order[level] = split_axis;
                    }

                    //   LOG::cout<<"Linear_Indices_Queue: "<<split_axis<<std::endl;
                    T_RANGE ranges[3] = {coordinate_range,coordinate_range,coordinate_range};
                    Split_Range(split_axis, ranges);

                    // number the interface
                    index -= fill_maps(ranges[2], level, --rank);
                    Number_Entry(ranges[2],index+1,0, 0);
                    if (rank==0)
                        fill_local_id(ranges[2],level);

                    // enqueue second half
                    queue.enqueue(T_QENTRY(ranges[1], level+1));
                    // enqueue first half
                    queue.enqueue( T_QENTRY(ranges[0], level+1));
                }
            }


        }
    template<int d>
    inline int REORDERING<d>::fill_maps(const T_RANGE& range,const int level,const int rank)
        {

            const T_INDEX range_size=abs(range.Edge_Lengths())+1;
            //LOG::cout<<range_size<<std::endl;
            const int N_entry = range_size.Product();
            level_sizes(level+1)+=N_entry;
            //int id = 0;
            for (ARBITRARY_RANGE_ITERATOR<d> iterator(range); iterator.Valid(); iterator.Next())
            {
                level_map(iterator.Index()) = level;
                rank_map(iterator.Index()) = rank;
                //        id_map(iterator.Index()) = id++;
            }


            return N_entry;
        }

    template<int d>
    inline void REORDERING<d>::fill_local_id(const T_RANGE& range,const int level)
    {
        PHYSBAM_ASSERT(range.min_corner.All_Less_Equal(range.max_corner));
        local_id_map[level].Resize(range);
        for (RANGE_ITERATOR<d> iterator(range); iterator.Valid(); iterator.Next())
        {
            local_id_map[level](iterator.Index()) = id_map(iterator.Index());
        }

    }

    template<int d>
    int REORDERING<d>::Get_First_Index(const int level, const int rank) const
        {
            PHYSBAM_ASSERT(level <= depth);
            PHYSBAM_ASSERT(rank < (1<<level));
            return use_queue?
                Get_Linear_Index_Queue(level,rank,0):
                Get_Linear_Index_Stack(level,rank,0);
        }

    template<int d>
    std::pair<int, int> REORDERING<d>::Get_Index_Range(const int level, const int id) const
        {
            PHYSBAM_ASSERT(id<(1<<level));
            //LOG::cout<<level<<" "<<id<<std::endl;
	    const int start_index = Get_First_Index(level,id);
	    const int end_index = start_index+level_sizes(level+1)-1;
            return std::pair<int, int>(start_index,end_index);
        }

    template<int d>
    void REORDERING<d>::print(const ARRAY<int, VECTOR<int,d> >& array, int offset)
        {
            for (RANGE_ITERATOR<d> iterator(array.domain); iterator.Valid(); iterator.Next())
            {
                const T_INDEX& index = iterator.Index();
                LOG::cout<<std::setw(max(offset,4))<<array(index);
                if (index(3)==array.domain.max_corner(3))
                {
                    LOG::cout<<std::endl;
                    for (int i=0; i<(offset-1)*((index(2)-array.domain.min_corner(2)+1)%((array.domain.Edge_Lengths())(2)+1)); i++)
                        LOG::cout<<" ";
                }
            }
        }



    template<int d>
    void REORDERING<d>::Set_Cut_Array(std::string cut_string)
    {
        const int number_of_cuts(cut_string.length());
        int total_cut[d];
        for (int v=0; v<d; v++)
            total_cut[v] = 0;
        cut_axes.Resize(number_of_cuts);

        for (int i=1; i<=number_of_cuts; i++) {
            const char cut_c = cut_string[i-1];
            const int cut = (int)(cut_c-'w');
            cut_axes(i)= cut;
            PHYSBAM_ASSERT(cut>0 && cut<=d);
            total_cut[cut-1]++;
        }

        for (int v=1; v<=d; v++)
        {
            int ds=size(v)>>1;
            while(ds!=0){
                total_cut[v-1]--;
                ds = ds>>1;
            }
            PHYSBAM_ASSERT(total_cut[v-1]==0);
        }
    }

    template<int d>
    int REORDERING<d>::Get_Linear_Index_Stack(const int level, const int rank, const int id) const
    {
        const int depth = level_sizes.Size() - 1;
        PHYSBAM_ASSERT(level <= depth);
        PHYSBAM_ASSERT(rank < (1<<level));
        PHYSBAM_ASSERT(id < level_sizes(level+1));
        int linear_index = id + 1;

        for (int l = 0; l<level; l++)
        {
            linear_index += (rank>>l) * level_sizes(level-l+1);
        }
        for (int l = level+1; l<=depth; l++)
        {
            linear_index += ((rank+1)<<(l-level)) * level_sizes(l+1);
        }

        return linear_index;
    }

    template<int d>
    int REORDERING<d>::Get_Linear_Index_Stack(const int level, const int rank, const int id, const int local_level) const
    {
        return Get_Linear_Index_Stack(level, rank&((1<<(level-local_level))-1),id);
    }

    template <int d>
    int REORDERING<d>::Get_Linear_Index_Queue(const int level, const int rank, const int id) const
    {
        const int depth = level_sizes.Size() - 1;
        PHYSBAM_ASSERT(level <= depth);
        PHYSBAM_ASSERT(rank < (1<<level));
        PHYSBAM_ASSERT(id < level_sizes(level+1));

        int linear_index = rank*level_sizes(level+1) + id + 1;

        for (int l = level+1; l<=depth; l++)
        {
            linear_index += level_sizes(l+1) * (1<<l);
        }
        return linear_index;
    }

    // template<int d>
    // RANGE<VECTOR<int,d> > REORDERING<d>::Level_Range(int level, int rank)
    // {
    //     const auto start = starting_coord(level+1)(rank+1);
    //     return T_RANGE(start,start+level_dims(level+1)-1);
    // }

    // template<int d>
    // RANGE<VECTOR<int,d> > Aggregate_Range(level);
     template<int d>
     template< class T>
    class REORDERING<d>::Circular_Queue
        {
            ARRAY<T> queue;
            int head;
            int tail;
        public:
        Circular_Queue(int q_length):
            queue(q_length), head(1), tail(1)
                {}

            inline const T& enqueue(const T& entry)
            {
                const int q_length =queue.Size();
                queue(tail++) = entry;
                if (tail == q_length+1) {tail = 1;}
                return entry;
            }

            inline const T& dequeue()
            {
                const int q_length =queue.Size();
                const T& current = queue(head++);
                if (head == q_length+1) {head = 1;}
                return current;
            }

            inline bool isEmpty() const
            {
                return head == tail;
            }
        };
    template class REORDERING<3>;
}
