#include "./Diag_Solve.h"
// used by macroblock solve
uint64_t master_ticks = 0;
uint64_t master_iterations = 0;

uint64_t level_4_ticks = 0;
uint64_t level_3_ticks = 0;
uint64_t level_2_ticks = 0;
uint64_t level_1_ticks = 0;
uint64_t level_0_ticks = 0;

uint64_t level_4_interface_ticks = 0;

uint64_t level_4_iterations = 0;
uint64_t level_3_iterations = 0;
uint64_t level_2_iterations = 0;
uint64_t level_1_iterations = 0;
uint64_t level_0_iterations = 0;

uint64_t xfer_ticks = 0;
uint64_t interface_ticks = 0;

namespace Hierarchical_Solve {
    void Print_Out_Result(MACROBLOCK::DOF_INTERIOR_DATA* x) {
        constexpr int d=3;
        // print out result
#if 1
        for (int rank=0; rank<16; rank++) {
            std::cout<<"rank "<<rank<<", L0 : "<<std::endl;
            for (int i=0; i<27; i++) {
                for (int v=0; v<d; v++)
                    std::cout<<x->L0[i][v][rank]<<" ";
                std::cout<<std::endl;
            }
            std::cout<<"rank "<<rank<<", L1 : "<<std::endl;
            for (int i=0; i<9; i++) {
                for (int v=0; v<d; v++)
                    std::cout<<x->L1[i][v][rank]<<" ";
                std::cout<<std::endl;
            }
            std::cout<<"rank "<<rank<<", L2 : "<<std::endl;
            for (int i=0; i<12; i++) {
                for (int v=0; v<d; v++)
                    std::cout<<x->L2[i][v][rank]<<" ";
                std::cout<<std::endl;
            }
            std::cout<<"rank "<<rank<<", L3 : "<<std::endl;
            for (int i=0; i<16; i++) {
                for (int v=0; v<d; v++)
                    std::cout<<x->L3[i][v][rank]<<" ";
                std::cout<<std::endl;
            }
            std::cout<<"rank "<<rank<<", L4 : "<<std::endl;
            for (int i=0; i<16; i++) {
                for (int v=0; v<d; v++)
                    std::cout<<x->L4[i][v][rank]<<" ";
                std::cout<<std::endl;
            }
            std::cout<<"rank "<<rank<<", L5 : "<<std::endl;
            for (int i=0; i<45; i++) {
                for (int v=0; v<d; v++)
                    std::cout<<x->L5[i][v][rank]<<" ";
                std::cout<<std::endl;
            }
        }
#endif
    }

}
