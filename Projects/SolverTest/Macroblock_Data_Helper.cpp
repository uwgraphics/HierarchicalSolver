#include "./Macroblock_Data_Helper.h"


namespace Macroblock_Utilities {

    void Clear_DOF(MACROBLOCK::DOF_INTERIOR_DATA* x) {
        constexpr int mb_level = 4;
// exchange x and z axis for macroblocks
        // clear out x
        for (int l=0; l<=mb_level; l++) {
            float* ptr = x->Level(l);
            std::fill(ptr, ptr+(x->LevelSize(l)/sizeof(float)), 0);
        }
        {
            float* ptr = &(x->L5[0][0][0]);
            std::fill(ptr, ptr+45*3*16, 0);
        }
    }

std::vector<int> Adjacent_Subdomain_Ranks(const int level, const int rank, const std::vector<int>& axes) {
        constexpr int d = 3;
        // axes - agg order of subdomains
        // construct masks and count bit
        const int axis = axes[level-1];
        int bit_count = 0;
        int axis_bit_count = 0;
        int add_carries = 0;
        int last_level = 0;
        int last_non_axis_level = 0;

        for (int l=1, s=1; l<level; l++, s<<=1) {
            const int dir = axes[l-1];
            if (dir != axis) {
                if (!bit_count)
                    add_carries|=s;
                bit_count++;
                last_non_axis_level = l;
            } else {
                last_level = l;
                if (bit_count)
                    add_carries|=s;
                axis_bit_count++;
            }
        }
        for (int l=level-1; l>last_non_axis_level; l--)
            add_carries&=(1<<l-1)-1;
        std::vector<int> adj_ranks(1<<(bit_count+1));
        //std::cout<<"add_carries = "<<add_carries<<std::endl;
#if 1
        for (int b=0,index=0; b<2; b++) {
            int adj_rank = rank<<level;
            adj_rank|=b<<(level-1);
            // be careful with level 1
            if (axis_bit_count)
                adj_rank|=1<<(last_level-1);
            for (int c=0; c<(1<<bit_count); c++) {
                adj_ranks[index++] = adj_rank;
                adj_rank+=add_carries;
            }
        }
#endif
        return adj_ranks;
    }

    int Primary_Cube(const std::array<int,3>& canonical_coord, const int cube_n) {
        /*
          Given the canonical coord of a node in cube_n, return the cube number of the 5x5x5 block that owns the data of the node
        */
        constexpr int d = 3;
        constexpr int mb_level = 4;
        constexpr std::array<int,mb_level> mb_axes{2,1,0,0}; // the agg order (in canonical coord) within the macroblock
        constexpr std::array<int,d> cube_size{5,5,5};

        int primary_n = cube_n;
        for (int v=0; v<d; v++)
            if (canonical_coord[v]==cube_size[v]-1) {
                for (unsigned int l=0, s=1<<(mb_level-1); l<mb_level; l++,s>>=1)
                    if (mb_axes[l]==v) {
                        // anil bit at s
                        primary_n &= ((1<<mb_axes.size())-1)^s;
                        break;
                    }

            } else if (canonical_coord[v] == 0) {
                bool found = false;
                for (unsigned int l=0, s=1<<(mb_level-1); l<mb_level; l++, s>>=1) {
                    if (!found){
                        if (mb_axes[l]==v && (primary_n&s))
                            found = true;
                    } else {
                        if (mb_axes[l] == v) {
                            // anil bit at l
                            primary_n &= ((1<<mb_level)-1)^s;
                            break;
                        }
                    }
                }

            }
        return primary_n;
    }

        void Cube_Coord_Mapping(std::array<int,3>& step, std::array<int,3>& min_corner, const int rank, const int cube_number, const std::vector<int> axes) {
        /* Given a 5x5x5 cube (identified by its cube# and the rank of the containing macroblock), return its global coord mapping
           axes - agg order of macroblocks should be yxzyxzyx...
        */
        // std::cout<<"cube_number = "<<cube_number<<std::endl;
        // std::cout<<"rank = "<<rank<<std::endl;
        constexpr int d=3;
        constexpr int mb_level = 4;
        constexpr std::array<int,d> sbd_size{3,3,3};
        constexpr std::array<int,mb_level> mb_axes{0,1,2,2}; // the agg order (in global coord) within the macroblock
        step={1,1,1};
        min_corner={0,0,0};
        std::array<int,d> box_size={sbd_size[0]+1,sbd_size[1]+1, sbd_size[2]+1};
        for (int s=1<<(mb_level-1), b=0; s; s>>=1, b++) {
            const int axis = mb_axes[b];
            box_size[axis]<<=1;
            if (s&cube_number) {
                min_corner[axis] = box_size[axis]-min_corner[axis];
                step[axis] =-step[axis];
            }
        }
        //   std::cout<<" box_size : ["<<box_size[0]<<","<<box_size[1]<<","<<box_size[2]<<"]"<<std::endl;
//        std::cout<<axes.size()<<std::endl;
        for (int r=rank, b=0; r; r>>=1, b++) {
            const int axis = axes[b];
            //std::cout<<"axis = "<<axis<<std::endl;
            box_size[axis]<<=1;
            if (r&1) {
                min_corner[axis] = box_size[axis]-min_corner[axis];
                step[axis] =-step[axis];
            }
        }
    }
    void Macroblock_Rank_Coord_Mapping(std::array<int,3>& step, std::array<int,3>& min_corner, std::array<int,3>& box_size, const int cube_n)
    {
        constexpr int d=3;
        constexpr int mb_level = 4;
        constexpr std::array<int,d> sbd_size{3,3,3};
        constexpr std::array<int,mb_level> mb_axes{0,1,2,2}; // the agg order (in global coord) within the macroblock
        step={1,1,1};
        min_corner={0,0,0};
        box_size={sbd_size[0]+1,sbd_size[1]+1, sbd_size[2]+1};
        for (int s=1<<(mb_level-1), b=0; s; s>>=1, b++) {
            const int axis = mb_axes[b];
            box_size[axis]<<=1;
            if (s&cube_n) {
                min_corner[axis] = box_size[axis]-min_corner[axis];
                step[axis] =-step[axis];
            }
        }
    }
    int Convert_Hierarchical_Sbd_and_Macroblock_Sbd(int sbd_number) {
        // the inverse of this function is itself
        int rank = 0;
        for (int s=4;s>0;s--,rank|=sbd_number&1,rank<<=1,sbd_number>>=1);
        rank>>=1;
        return rank;
    }

    void Inplace_Conjugate_L(Matrix4<float,float>& diag) {
        // in-place multiplication for inv(D) = inv(L)'*inv(L)
        for (int j=0; j<4; j++) {
            for (int i=0; i<j; i++)
                diag.data[j*4+i] = diag.data[i*4+j];
            for (int i=j; i<4; i++) {
                float tmp = 0;
                for (int v=i; v<4; v++)
                    tmp+=diag.data[i*4+v]*diag.data[j*4+v];
                diag.data[j*4+i] = tmp;
            }
        }
    }
}
