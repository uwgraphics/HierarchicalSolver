#include "macro_block_solver.h"
#include <array>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <utility>

#include "Macroblock_Utilities/Common/KernelCommon.h"
#include "Macroblock_Utilities/Kernels/Weighted_Gradient/Weighted_Gradient.h"
#include "Macroblock_Utilities/Kernels/Weighted_Accumulation/Weighted_Accumulation.h"
#include "Macroblock_Utilities/Kernels/Singular_Value_Decomposition/Singular_Value_Decomposition.h"
#include "Macroblock_Utilities/Kernels/Add_Force_Single_QPoint/Add_Force_Single_QPoint.h"
#include "Macroblock_Utilities/Kernels/Update_Position_Based_State/Update_Position_Based_State.h"
#include "Macroblock_Utilities/Kernels/Matrix_Times_Transpose/Matrix_Times_Transpose.h"
#include "Macroblock_Utilities/Kernels/Compute_Cell_Matrix/Compute_Cell_Matrix.h"
#include "Macroblock_Utilities/Sparse_Blocked_CSR_Multiply_Metaprogram.h"

#include "Macroblock_Utilities/Matrix3.h"
#include "Macroblock_Utilities/SymmetricMatrix3.h"

#include "signal.h"

//#include <Kernels/Weighted_Acculmuation.h>
// Some other stuff.

void Sparse_Blocked_CSR_Masked_Multiply_Flat(const BlockedCSRMatrix3<float[16]> &U,
                                             const typename std::remove_extent<float[16]>::type *b,
                                             typename std::remove_extent<float[16]>::type *x,
                                             const float (&mask)[16],
                                             bool transpose);

extern uint64_t master_ticks;
extern uint64_t master_iterations;

extern uint64_t level_4_ticks;
extern uint64_t level_3_ticks;
extern uint64_t level_2_ticks;
extern uint64_t level_1_ticks;
extern uint64_t level_0_ticks;

extern uint64_t level_4_interface_ticks;

extern uint64_t level_4_iterations;
extern uint64_t level_3_iterations;
extern uint64_t level_2_iterations;
extern uint64_t level_1_iterations;
extern uint64_t level_0_iterations;

extern uint64_t xfer_ticks;
extern uint64_t interface_ticks;


struct COROTATED_TAG;

namespace std
{
    template<>
    struct hash<std::pair<int, int>>
    {
        size_t operator () (std::pair<int, int> const& p) const
        {
            return (std::hash<int>()(p.first)*2397917 + std::hash<int>()(p.second)*199777 );
        }
    };
}


namespace MACROBLOCK {

    void inline CollapseLevel(DOF_INTERIOR_DATA& data, int level);
    void inline PruneLevel(DOF_INTERIOR_DATA& data, int level);
    void inline DistributeLevel(DOF_INTERIOR_DATA& data, int level);
    void inline NegateLevel(DOF_INTERIOR_DATA& data, int level, bool copydown);
    void inline PlusEqualsLevel(DOF_INTERIOR_DATA& data, int level, bool copydown);

    template<int level, bool in_place>
    void Solve_Level(const MATRIX_STRUCTURE& matrix_structure, const MATRIX_DATA& matrix_data,  DOF_INTERIOR_DATA& y, DOF_INTERIOR_DATA& x, DOF_SCRATCH_DATA& scratch)
    {

        uint64_t start_ticks = __rdtsc();

        // Setup

        const BlockedCSRMatrix3< T_DATA > TransferMatrix = {matrix_structure.TransferRows(level),
                                                            matrix_structure.TransferColumns(level),
                                                            (void*)(matrix_data.TransferPointer(level)),
                                                            matrix_structure.TransferOffsets(level),
                                                            matrix_structure.TransferColumn(level)};

        // Algorithm
        // Notes: Interface variables are all variables at the current level
        //        Non-Interface variables are all variables at less than the current level
        //        We can ignore all variables greater than the current level
        //        Original refers to input y variables
        //        RHS refers to Interface+NonInterface variables


        // M5
        //
        // In: Unmodified Original Non-Interface variables
        //
        // Out: Modified Non-Interface Variables,
        //      Original Interface Variables
        //

        Solve_Level<level-1,true>(matrix_structure, matrix_data, y, scratch.Level(level-1), scratch);
        for(int l = level-1; l>=0; l--)
            DistributeLevel( scratch.Level(level-1), l );

        // M4
        //
        // In: Modified Non-Interface variables from M5, Original RHS
        //
        // Out: Original NonInterface variables, Modified Interface + Original Interface
        //
        //
        uint64_t xfer_start_ticks = __rdtsc();

        for( int i=0; i<VECTOR_WIDTH; i++){
            scratch.mask[i]=1.0f;
            if(level==4 && ( (i / 2) % 2 == 0 ) )
                scratch.mask[i] = 0.0f;
        }

        for( int i = 0; i<VECTOR_WIDTH; i+=INSTRUCTION_WIDTH){
                Sparse_Blocked_CSR_Masked_Multiply<INSTRUCTION_TYPE,T_DATA,I_DATA>( TransferMatrix.Offset(i),
                                                                                    scratch.Level(level-1).Level(0)+i,
                                                                                    scratch.Level(level).Level(level)+i,
                                                                                    reinterpret_cast<C_T_DATA_SCALAR>(scratch.mask[i]),
                                                                                    false);
        }

        uint64_t xfer_stop_ticks = __rdtsc();
        xfer_ticks += (xfer_stop_ticks - xfer_start_ticks);

        // Negate the results of the Transfer matrix (instead of having a negative matrix )
        NegateLevel(scratch.Level(level), level, false);
        // Add to results the saved copy of the orginal interface
        PlusEqualsLevel(scratch.Level(level), y, level, false);
        // Collapse the now modified and combined interface values
        CollapseLevel( scratch.Level(level), level );


        // M3
        //
        // In: Modified Interface from M4
        //
        // Out: Passed through Non-Interface from M4, Modified Interface
        //

        uint64_t interface_start_ticks = __rdtsc();

        for( int matrix_num=0;matrix_num<matrix_structure.SigmaCount(level);matrix_num++ ){
            const BlockedMatrixNXN<T> SigmaMatrix = {matrix_structure.SigmaN(level),
                                                     (void*)matrix_data.SigmaDiagonal(level,matrix_num),
                                                     (void*)matrix_data.SigmaOffDiagonal(level,matrix_num) };
            CopyInSigmaData<level>( scratch, scratch.Level(level), matrix_num );
            Block_Forward_Substitution<INSTRUCTION_TYPE,float>(SigmaMatrix, scratch.SigmaScratch, true);
            Block_Backward_Substitution<INSTRUCTION_TYPE,float>(SigmaMatrix, scratch.SigmaScratch, false);
            CopyOutSigmaData<level>( scratch, scratch.Level(level), matrix_num );
        }
        uint64_t interface_stop_ticks = __rdtsc();
        if(level == 4) level_4_interface_ticks += (interface_stop_ticks - interface_start_ticks);
        interface_ticks += (interface_stop_ticks - interface_start_ticks);

        // Collapse the now modified and combined interface values
        DistributeLevel( scratch.Level(level), level );



        // M2
        //
        // In: Interface from M3, Passed through Non-Interface from M4
        //
        // Out: Interface from M3, Modified Non-Interface + Non-Interface from M4
        //
        //
        xfer_start_ticks = __rdtsc();

        for( int i = 0; i<VECTOR_WIDTH; i+=INSTRUCTION_WIDTH){
            Sparse_Blocked_CSR_Multiply<INSTRUCTION_TYPE,T_DATA,I_DATA>( TransferMatrix.Offset(i),
                                                                         scratch.Level(level).Level(level)+i,
                                                                         scratch.Level(level).Level(0)+i,
                                                                         true);
        }

        xfer_stop_ticks = __rdtsc();
        xfer_ticks += (xfer_stop_ticks - xfer_start_ticks);

        // Negate the results of the Transfer matrix (instead of having a negative matrix )
        NegateLevel(scratch.Level(level),level-1,true);
        // Add to results the saved copy of the M4 non-interface
        PlusEqualsLevel(scratch.Level(level),y,level-1,true);
        // Collapse
        for(int l = level-1; l>=0; l--)
            CollapseLevel( scratch.Level(level), l );
        PruneLevel( scratch.Level(level), level );

        // M1
        //
        // In: Modified Non-Interface from M2, Passed through Interface from M3
        //
        // Out: Solved LHS
        //
        Solve_Level<level-1,false>(matrix_structure, matrix_data, scratch.Level(level), scratch.Level(level-1), scratch);

        // Finally, copy LHS to x
        CopyDown( scratch.Level(level-1), x, level-1 );
        CopyAt( scratch.Level(level), x, level );

        uint64_t end_ticks = __rdtsc();
        switch( level ){
        case 1:
            level_1_ticks += (end_ticks - start_ticks);
            level_1_iterations++;
            break;
        case 2:
            level_2_ticks += (end_ticks - start_ticks);
            level_2_iterations++;
            break;
        case 3:
            level_3_ticks += (end_ticks - start_ticks);
            level_3_iterations++;
            break;
        case 4:
            level_4_ticks += (end_ticks - start_ticks);
            level_4_iterations++;
            break;
        }

    }

    template void Solve_Level<1,false>(const MATRIX_STRUCTURE& matrix_structure, const MATRIX_DATA& matrix_data,  DOF_INTERIOR_DATA& y, DOF_INTERIOR_DATA& x, DOF_SCRATCH_DATA& scratch);
    template void Solve_Level<1,true>(const MATRIX_STRUCTURE& matrix_structure, const MATRIX_DATA& matrix_data,  DOF_INTERIOR_DATA& y, DOF_INTERIOR_DATA& x, DOF_SCRATCH_DATA& scratch);

    template void Solve_Level<2,false>(const MATRIX_STRUCTURE& matrix_structure, const MATRIX_DATA& matrix_data,  DOF_INTERIOR_DATA& y, DOF_INTERIOR_DATA& x, DOF_SCRATCH_DATA& scratch);
    template void Solve_Level<2,true>(const MATRIX_STRUCTURE& matrix_structure, const MATRIX_DATA& matrix_data,  DOF_INTERIOR_DATA& y, DOF_INTERIOR_DATA& x, DOF_SCRATCH_DATA& scratch);

    template void Solve_Level<3,false>(const MATRIX_STRUCTURE& matrix_structure, const MATRIX_DATA& matrix_data,  DOF_INTERIOR_DATA& y, DOF_INTERIOR_DATA& x, DOF_SCRATCH_DATA& scratch);
    template void Solve_Level<3,true>(const MATRIX_STRUCTURE& matrix_structure, const MATRIX_DATA& matrix_data,  DOF_INTERIOR_DATA& y, DOF_INTERIOR_DATA& x, DOF_SCRATCH_DATA& scratch);

    template void Solve_Level<4,true>(const MATRIX_STRUCTURE& matrix_structure, const MATRIX_DATA& matrix_data,  DOF_INTERIOR_DATA& y, DOF_INTERIOR_DATA& x, DOF_SCRATCH_DATA& scratch);

    template<>
    void Solve_Level<0, true>(const MATRIX_STRUCTURE& matrix_structure, const MATRIX_DATA& matrix_data,  DOF_INTERIOR_DATA& y, DOF_INTERIOR_DATA& x, DOF_SCRATCH_DATA& scratch)
    {

        //*************************************
        //      Solve Subdomains exactly
        //
        //*************************************
        // Solve_Subdomains only changes values corresponding to subdomain variables
        // and retains values for interface variables (representing multiplication
        // by identity for interface values.

        uint64_t start_ticks = __rdtsc();

        CopyAt( y, x, 0 );
        BlockedCSCSymmetricMatrix3< T_DATA > Subdomain_Factorization = { MATRIX_STRUCTURE::L0_n,
                                                                         (void*)(matrix_data.L0_diagonal_entries),
                                                                         (void*)(matrix_data.L0_offdiagonal_entries),
                                                                         (int*)(matrix_structure.L0_offsets),
                                                                         (int*)(matrix_structure.L0_rows)};





        for( int i = 0; i<VECTOR_WIDTH; i+=INSTRUCTION_WIDTH){
            Sparse_Blocked_Forwards_Substitution<INSTRUCTION_TYPE,T_DATA,I_DATA>( Subdomain_Factorization.Offset(i),
                                                                                  x.Level(0)+i,
                                                                                  true);

            Sparse_Blocked_Backwards_Substitution<INSTRUCTION_TYPE,T_DATA,I_DATA>( Subdomain_Factorization.Offset(i),
                                                                                   x.Level(0)+i,
                                                                                   false);
        }

        uint64_t end_ticks = __rdtsc();
        level_0_ticks += (end_ticks - start_ticks);
        level_0_iterations++;

    }

    template<>
    void Solve_Level<0, false>(const MATRIX_STRUCTURE& matrix_structure, const MATRIX_DATA& matrix_data,  DOF_INTERIOR_DATA& y, DOF_INTERIOR_DATA& x, DOF_SCRATCH_DATA& scratch)
    {
        Solve_Level<0, true>( matrix_structure, matrix_data, y, x, scratch);
    }

/****************************************************************************


                       Initialize_Matrix_Structure


*///*************************************************************************

    void Initialize_Matrix_Structure( MATRIX_STRUCTURE& matrix_structure ){

        std::array< std::array<int,3>, 8 > flatnode_to_index = {{ { 0,0,0 },
                                                                  { 0,0,1 },
                                                                  { 0,1,0 },
                                                                  { 0,1,1 },
                                                                  { 1,0,0 },
                                                                  { 1,0,1 },
                                                                  { 1,1,0 },
                                                                  { 1,1,1 } }};

        std::array< std::array< std::array<int,3>, 2>, 7> canonical_levels = {{ {{ {1,1,1}, {3,3,3} }},     // L0 3x3x3
                                                                                {{ {1,1,4}, {3,3,4} }},     // L1 3x3x1
                                                                                {{ {1,4,1}, {3,4,4} }},     // L2 3x1x4
                                                                                {{ {4,1,1}, {4,4,4} }},     // L3 1x4x4
                                                                                {{ {0,1,1}, {0,4,4} }},     // L4 1x4x4
                                                                                {{ {0,0,1}, {4,0,4} }},     // L5 5x1x4
                                                                                {{ {0,0,0}, {4,4,0} }} }};  // L6 5x5x1

        std::array< std::array< std::array<int,3>, 2>, 16> subdomains = {{ {{ {0, 0, 0}, {3, 3, 3} }},
                                                                           {{ {15, 0, 0}, {12, 3, 3} }},
                                                                           {{ {7, 0, 0}, {4, 3, 3} }},
                                                                           {{ {8, 0, 0}, {11, 3, 3} }},
                                                                           {{ {0, 7, 0}, {3, 4, 3} }},
                                                                           {{ {15, 7, 0}, {12, 4, 3} }},
                                                                           {{ {7, 7, 0}, {4, 4, 3} }},
                                                                           {{ {8, 7, 0}, {11, 4, 3} }},
                                                                           {{ {0, 0, 7}, {3, 3, 4} }},
                                                                           {{ {15, 0, 7}, {12, 3, 4} }},
                                                                           {{ {7, 0, 7}, {4, 3, 4} }},
                                                                           {{ {8, 0, 7}, {11, 3, 4} }},
                                                                           {{ {0, 7, 7}, {3, 4, 4} }},
                                                                           {{ {15, 7, 7}, {12, 4, 4} }},
                                                                           {{ {7, 7, 7}, {4, 4, 4} }},
                                                                           {{ {8, 7, 7}, {11, 4, 4} }} }};


        for( int subdomain = 0; subdomain < 16; subdomain++ ){
            int flat_index = 0;
            std::array<int,3> dir = { subdomains[subdomain][0][0]<subdomains[subdomain][1][0]?1:-1,
                                      subdomains[subdomain][0][1]<subdomains[subdomain][1][1]?1:-1,
                                      subdomains[subdomain][0][2]<subdomains[subdomain][1][2]?1:-1 };
            for( int i=subdomains[subdomain][0][0];i!=subdomains[subdomain][1][0]+dir[0]; i+=dir[0])
                for( int j=subdomains[subdomain][0][1];j!=subdomains[subdomain][1][1]+dir[1]; j+=dir[1])
                    for( int k=subdomains[subdomain][0][2];k!=subdomains[subdomain][1][2]+dir[2]; k+=dir[2]){
                        matrix_structure.cell_index[i][j][k] = flat_index;
                        matrix_structure.cell_subdomain[i][j][k] = subdomain;
                        flat_index++;
                    }
        }


        int dof_to_level [5][5][5];
        std::array< int, 3> mirrored_reordering_helper = { 1, 3, 2 };
        std::vector< std::array<int,3> > index_lists;
        for( int level = 0; level <= 6; level++ )
            for( int i=canonical_levels[level][0][0];i<=canonical_levels[level][1][0];i++)
                for( int j=canonical_levels[level][0][1];j<=canonical_levels[level][1][1];j++)
                    for( int k=canonical_levels[level][0][2];k<=canonical_levels[level][1][2];k++){
                        std::array<int, 3> index;
                        index[0] = i; index[1] = j; index[2] = k;
                        if( level == 0 ){
                            index[0] = mirrored_reordering_helper[i-1];
                            index[1] = mirrored_reordering_helper[j-1];
                            index[2] = mirrored_reordering_helper[k-1];
                        }
                        index_lists.push_back( index );
                        dof_to_level[i][j][k] = level;
                        if( level == 6 )
                            dof_to_level[i][j][k] = 5;
                    }
        for( std::vector< std::array<int,3> >::const_iterator iter=index_lists.begin();
             iter != index_lists.end();
             iter++ ){
            const std::array<int, 3>& index=*iter;
            matrix_structure.dof_index[index[0]][index[1]][index[2]] = (iter-index_lists.begin());
        }

        int pairing_sentinal = -1;
        int pairing_count = 0;
        int pairing_index[125][125];
        int done_pairing_index[125][125];

        for( int i=0;i<125;i++)
            for( int j=0;j<125;j++)
                pairing_index[i][j] = pairing_sentinal;

        // Block construction

        // Need this later...
        std::unordered_set< std::pair<int,int> > L0_sparcity;

        for( int block_x =0; block_x < 6; block_x++){
            int block_y = block_x;

            for( int node_v1 = 0; node_v1 < 125; node_v1++){
                for( int node_v2 = 0; node_v2 < node_v1 ; node_v2++){
                    int node_v1_level = dof_to_level[index_lists[node_v1][0]]
                                                    [index_lists[node_v1][1]]
                                                    [index_lists[node_v1][2]];
                    int node_v2_level = dof_to_level[index_lists[node_v2][0]]
                                                    [index_lists[node_v2][1]]
                                                    [index_lists[node_v2][2]];
                    if( node_v1_level != block_x || node_v2_level != block_y ) continue;

                    std::array<int,3> node_v1_index = index_lists[node_v1];
                    std::array<int,3> node_v2_index = index_lists[node_v2];
                    std::array<int,3> node_index_difference;
                    bool connected = true;
                    for( int w = 0; w < 3; w ++){
                        node_index_difference[w] = abs( node_v1_index[w] - node_v2_index[w] );
                        if( node_index_difference[w] > 1 )
                            connected = false;
                    }
                    if(connected){
                        pairing_index[node_v1][node_v2] = pairing_count++;
                        if( block_x == 0 )
                            L0_sparcity.insert( std::pair<int,int>( node_v1, node_v2 ) );
                    }
                }

            }
        }

        for( int block_x=1; block_x < 6; block_x++){
            int* Transfer_Offsets = matrix_structure.TransferOffsets(block_x);
            int* Transfer_Column = matrix_structure.TransferColumn(block_x);

            int row = 0;
            int offset = 0;
            for( int node_v1 = 0; node_v1 < 125; node_v1++){
                Transfer_Offsets[row] = offset;
                int node_v1_level = dof_to_level[index_lists[node_v1][0]]
                    [index_lists[node_v1][1]]
                    [index_lists[node_v1][2]];
                int column = 0;
                for( int block_y = 0; block_y < block_x; block_y++){
                    for( int node_v2 = 0; node_v2 < node_v1 ; node_v2++){
                        int node_v2_level = dof_to_level[index_lists[node_v2][0]]
                            [index_lists[node_v2][1]]
                            [index_lists[node_v2][2]];
                        if( node_v1_level != block_x || node_v2_level != block_y ) continue;

                        std::array<int,3> node_v1_index = index_lists[node_v1];
                        std::array<int,3> node_v2_index = index_lists[node_v2];
                        std::array<int,3> node_index_difference;
                        bool connected = true;
                        for( int w = 0; w < 3; w ++){
                            node_index_difference[w] = abs( node_v1_index[w] - node_v2_index[w] );
                            if( node_index_difference[w] > 1 )
                                connected = false;
                        }
                        if(connected){
                            pairing_index[node_v1][node_v2] = pairing_count++;
                            Transfer_Column[offset] = column;
                            offset++;
                        }

                        if(node_v2_level == block_y) column++;
                    }
                }
                if(node_v1_level == block_x) row++;
            }
            Transfer_Offsets[row] = offset;
        }

        // element_coeff Construction

        for( int ci = 0; ci < 4; ci ++)
            for( int cj = 0; cj < 4; cj ++)
                for( int ck = 0; ck < 4; ck ++){
                    for( int v1 = 0; v1<8; v1++ )
                        for( int v2 = 0; v2 < 8; v2++){
                            std::array<int,3> index_v1 = {{ci,cj,ck}};
                            std::array<int,3> index_v2 = {{ci,cj,ck}};
                            for( int w=0; w<3; w++ ){
                                index_v1[w] += flatnode_to_index[v1][w];
                                index_v2[w] += flatnode_to_index[v2][w];
                            }
                            int dof_index_v1 = matrix_structure.dof_index[index_v1[0]][index_v1[1]][index_v1[2]];
                            int dof_index_v2 = matrix_structure.dof_index[index_v2[0]][index_v2[1]][index_v2[2]];
                            matrix_structure.element_coeff[ci][cj][ck][v1][v2] = pairing_index[dof_index_v1][dof_index_v2];
                        }
                }


        // Subdomain factorization structure
        {
            // Extend Sparcity
            for( int j = 0; j < 27; j++ ) // For each column
                for( int l = j+1; l < 27; l++ ){ // Examine each potential row
                    if( L0_sparcity.find( std::pair<int,int>(l,j) )  == L0_sparcity.end() )
                        continue; // Skip over rows without entries
                    for( int r = l+1, k = l+1; k < 27; k++ ){ // Now we look at pairs in the columns of j and l
                        if( L0_sparcity.find( std::pair<int,int>(k,j) ) == L0_sparcity.end() )
                            continue; // Skip over rows without entries
                        for(;r!=k;r++) // Real work here.
                            assert(  r < 27 ); // We can't go over the size of the matrix...
                        if( L0_sparcity.find( std::pair<int,int>(r,l) ) == L0_sparcity.end() )
                            L0_sparcity.insert( std::pair<int,int>(r,l) );
                    }
                }
            // Generate Offsets and Rows...
            int* L0_offsets = matrix_structure.L0_offsets;
            int* L0_rows = matrix_structure.L0_rows;
            int count = 0;
            for( int j = 0; j < 27; j++ ){ // For each column
                L0_offsets[j] = count;
                for( int l = j+1; l < 27; l++ ){ // Examine each potential row
                    if( L0_sparcity.find( std::pair<int,int>(l,j) )  == L0_sparcity.end() )
                        continue; // Skip over rows without entries
                    L0_rows[count] = l;
                    count++;
                }
            }
            L0_offsets[27] = count;
        }
    }

    #if 0
/****************************************************************************


                        Update Position Based State


*///*************************************************************************

    void Update_Position_Based_State_And_Add_Force( const MATRIX_STRUCTURE& matrix_structure, const MATRIX_PARAMETERS& parameters, MATRIX_DATA& matrix_data, DOF_INTERIOR_DATA& X, DOF_INTERIOR_DATA& forces ){

        std::array< std::array<int,3>, 8 > flatnode_to_index = {{ { 0,0,0 },
                                                                  { 0,0,1 },
                                                                  { 0,1,0 },
                                                                  { 0,1,1 },
                                                                  { 1,0,0 },
                                                                  { 1,0,1 },
                                                                  { 1,1,0 },
                                                                  { 1,1,1 } }};
        const float one_over_root_three = 0.57735027;
        enum {d=3};
        typedef std::array< int, 3> T_INDEX;
        typedef float (&T_DATA_ref)[16];
        typedef Number<INSTRUCTION_TYPE> Tn;

        typedef float (&fArr_W)[VECTOR_WIDTH];
        typedef const float (&cfArr_W)[VECTOR_WIDTH];

        typedef int (&iArr_W)[VECTOR_WIDTH];
        typedef const int (&ciArr_W)[VECTOR_WIDTH];

        typedef float (&fArr_3_8_W)[3][8][VECTOR_WIDTH];
        typedef const float (&cfArr_3_8_W)[3][8][VECTOR_WIDTH];

        typedef float (&fArr_3_W)[3][VECTOR_WIDTH];
        typedef const float (&cfArr_3_W)[3][VECTOR_WIDTH];

        typedef float (&fArr_9_W)[9][VECTOR_WIDTH];
        typedef const float (&cfArr_9_W)[9][VECTOR_WIDTH];

        typedef float (&fArr_6_W)[6][VECTOR_WIDTH];
        typedef const float (&cfArr_6_W)[6][VECTOR_WIDTH];

        typedef float (&fArr_12_W)[12][VECTOR_WIDTH];
        typedef const float (&cfArr_12_W)[12][VECTOR_WIDTH];

#define KCAST(T_ARCH_TYPE, VAR) reinterpret_cast<T_ARCH_TYPE>(VAR)

        alignas(64) T_DATA Du [3][8];
        alignas(64) T_DATA Hu [3][8];
        alignas(64) T_DATA Weights [8][3];
        alignas(64) T_DATA P [9];
        alignas(64) T_DATA one_over_h[3];
        alignas(64) T_DATA scale;
        alignas(64) T_DATA one;
        alignas(64) T_DATA cell_volume;
        alignas(64) T_DATA cutoff;
        alignas(64) T_DATA U[9];
        alignas(64) T_DATA V[9];
        alignas(64) T_DATA Sigma[3];
        alignas(64) T_DATA dPdF[12];
        alignas(64) T_DATA mu;
        alignas(64) T_DATA lambda;

        // Initialize 'Global' Variables
        {int qp=0;
            for(int _i=0; _i<2; _i++)
                for(int _j=0; _j<2; _j++)
                    for(int _k=0; _k<2; _k++, qp++){
                        const T_INDEX index = {_i, _j, _k};
                        for(int i=0;i<d;i++)
                            for( int c=0;c<16;c++ )
                                if(index[i]==0)  Weights[qp][i][c]=(1-one_over_root_three)*.5;
                                else             Weights[qp][i][c]=(1+one_over_root_three)*.5;
                    }
        }
        for( int c=0;c<16;c++){
            for(int q=0;q<3;q++)
                one_over_h[q][c] = 1.0f/parameters.dx;
            if(c&0x1)
                one_over_h[0][c] *= -1.f;
            if(c&0x2)
                one_over_h[0][c] *= -1.f;
            if(c&0x4)
                one_over_h[1][c] *= -1.f;
            if(c&0x8)
                one_over_h[2][c] *= -1.f;
            scale[c] = -1.0f * parameters.dx*parameters.dx*parameters.dx / 8;
            cell_volume[c] = parameters.dx*parameters.dx*parameters.dx;
            cutoff[c] = 1e-4f;
            one[c] = 1.0f;
       }
        memset( Hu, 0, sizeof( Hu ) );
        memset( P, 0, sizeof( P ) );
        memset( U, 0, sizeof( U ) );
        memset( V, 0, sizeof( V ) );
        memset( Sigma, 0, sizeof( Sigma ) );
        memset( dPdF, 0, sizeof( dPdF ) );

        //raise (SIGUSR1);

        for( int cell_i=0;cell_i<4; cell_i++)
            for( int cell_j=0;cell_j<4; cell_j++)
                for( int cell_k=0;cell_k<4; cell_k++){

                    memcpy( mu, parameters.mu[cell_i*16+cell_j*4+cell_k], sizeof( T_DATA ) );
                    memcpy( lambda, parameters.lambda[cell_i*16+cell_j*4+cell_k], sizeof( T_DATA ) );

                    // Initialize Du
                    for( int v1 = 0; v1 < 8; v1++){
                        T_INDEX node_index = {cell_i,cell_j,cell_k};
                        for( int w= 0; w< 3; w++ )
                            node_index[w] += flatnode_to_index[v1][w];
                        int node_flatindex = matrix_structure.dof_index[node_index[0]][node_index[1]][node_index[2]];

                        for( int w= 0; w< 3; w++ ){
                            memcpy( &(Du[w][v1][0]), X.Level(0)+( (node_flatindex*16*3)+(w*16) ), sizeof(float)*VECTOR_WIDTH  );
                        }
                    }

                    // Generate forces for each QP
                    for( int qp=0; qp<8; qp++) {
                        memset( Hu, 0, sizeof( Hu ) ); // Super important!!! - can't accumulate onto old H variable from before...

                        // Build U, V, Sigma, dPdF
                        for( int off = 0; off<VECTOR_WIDTH; off+=INSTRUCTION_WIDTH){
                            ::Update_Position_Based_State<COROTATED_TAG,INSTRUCTION_TYPE,T_DATA,I_DATA>::Run( KCAST(cfArr_3_8_W,Du[0][0][off]),
                                                                                                              KCAST(cfArr_W,mu[off]),
                                                                                                              KCAST(cfArr_W,lambda[off]),
                                                                                                              KCAST(cfArr_3_W,Weights[qp][0][off]),
                                                                                                              KCAST(cfArr_W,cutoff[off]),
                                                                                                              KCAST(cfArr_3_W,one_over_h[0][off]),
                                                                                                              KCAST(cfArr_W,cell_volume[off]),
                                                                                                              KCAST(fArr_9_W,U[0][off]),
                                                                                                              KCAST(fArr_9_W,V[0][off]),
                                                                                                              KCAST(fArr_3_W,Sigma[0][off]),
                                                                                                              KCAST(fArr_12_W,dPdF[0][off]));

                        // Compute H
                            Transpose<INSTRUCTION_TYPE,T_DATA,I_DATA>(KCAST(fArr_9_W, P[0][off]),
                                                                      KCAST(cfArr_9_W, V[0][off]));
                            Weighted_Accumulation<INSTRUCTION_TYPE,T_DATA,I_DATA>( KCAST(fArr_3_8_W, Hu[0][0][off]),
                                                                                   KCAST(cfArr_9_W, P[0][off]),
                                                                                   KCAST(cfArr_3_W, Weights[qp][0][off]),
                                                                                   KCAST(cfArr_3_W, one_over_h[0][off]),
                                                                                   KCAST(cfArr_W, one[off]) );

                            for( int i = 0; i < 8; i++){
                                for( int k = 0; k < 8; k++){
                                    T_INDEX node_index = {cell_i,cell_j,cell_k};
                                    for( int w= 0; w< 3; w++ ){
                                        node_index[w] += flatnode_to_index[i][w];
                                    }
                                    int node_flatindex = matrix_structure.dof_index[node_index[0]][node_index[1]][node_index[2]];
                                    int location = matrix_structure.element_coeff[cell_i][cell_j][cell_k][i][k];
                                    if( location != -1 || i == k ){
                                        Build_M<INSTRUCTION_TYPE,T_DATA,I_DATA>(KCAST(fArr_9_W, P[0][off]),
                                                                                KCAST(cfArr_12_W,dPdF[0][off]),
                                                                                KCAST(cfArr_3_8_W,Hu[0][0][off]),
                                                                                KCAST(cfArr_W,scale[off]),
                                                                                i, k);

                                        Matrix_Times_Transpose<INSTRUCTION_TYPE,T_DATA,I_DATA>(KCAST(cfArr_9_W,U[0][off]),
                                                                                               KCAST(cfArr_9_W,P[0][off]),
                                                                                               KCAST(fArr_9_W,P[0][off]));
                                        Matrix_Times_Transpose<INSTRUCTION_TYPE,T_DATA,I_DATA>(KCAST(cfArr_9_W,U[0][off]),
                                                                                               KCAST(cfArr_9_W,P[0][off]),
                                                                                               KCAST(fArr_9_W,P[0][off]));

                                        if(location != -1){
                                            // Perform a plus equals here...
                                            Matrix3<Tn> old;
                                            Matrix3<Tn> update;
                                            old.Load_Aligned( KCAST(fArr_9_W, *(matrix_data.A_Offdiagonals() + location*(9*VECTOR_WIDTH) + off )) );
                                            update.Load_Aligned( KCAST(fArr_9_W, P[0][off] ) );
                                            old += update;
                                            old.Store( KCAST(fArr_9_W, *(matrix_data.A_Offdiagonals() + location*(9*VECTOR_WIDTH)  + off )) );
                                        }
                                        else{
                                            // Perform a plus equals here... But its a little more complicated because we have a M9 -> M6
                                            SymmetricMatrix3<Tn> old;
                                            SymmetricMatrix3<Tn> update;
                                            old.Load_Aligned( KCAST(fArr_6_W, *(matrix_data.A_Diagonals() + node_flatindex*(6*VECTOR_WIDTH)  + off )) );
                                            update.Load_Aligned( KCAST(fArr_9_W, P[0][off] ) );
                                            old += update;
                                            old.Store( KCAST(fArr_6_W, *(matrix_data.A_Diagonals() + node_flatindex*(6*VECTOR_WIDTH)  + off )) );
                                        }

                                    }

                                }
                            }


                        // Compute P
                            Add_Force_Single_QPoint<COROTATED_TAG,INSTRUCTION_TYPE,T_DATA,I_DATA>::Run( KCAST(cfArr_W, mu[off]),
                                                                                                        KCAST(cfArr_W, lambda[off]),
                                                                                                        KCAST(cfArr_9_W, U[0][off]),
                                                                                                        KCAST(cfArr_9_W, V[0][off]),
                                                                                                        KCAST(cfArr_3_W, Sigma[0][off]),
                                                                                                        KCAST(fArr_9_W, P[0][off]) );
                            memset( Hu, 0, sizeof( Hu ) ); // Super important!!! - can't accumulate onto old H variable from before...
                        // Compute Forces
                            Weighted_Accumulation<INSTRUCTION_TYPE,T_DATA,I_DATA>( KCAST(fArr_3_8_W, Hu[0][0][off]),
                                                                                   KCAST(cfArr_9_W, P[0][off]),
                                                                                   KCAST(cfArr_3_W, Weights[qp][0][off]),
                                                                                   KCAST(cfArr_3_W, one_over_h[0][off]),
                                                                                   KCAST(cfArr_W, scale[off]) );



                        // Accumulate
                            for( int v1 = 0; v1 < 8; v1++){
                                T_INDEX node_index = {cell_i,cell_j,cell_k};
                                for( int w= 0; w< 3; w++ )
                                    node_index[w] += flatnode_to_index[v1][w];
                                int node_flatindex = matrix_structure.dof_index[node_index[0]][node_index[1]][node_index[2]];

                                for( int w= 0; w< 3; w++ ){
                                    for( int c=off; c<off+INSTRUCTION_WIDTH; c++ )
                                        forces.Level(0)[( (node_flatindex*16*3)+(w*16)+c )] += Hu[w][v1][c];
                                }
                            }


                        }

                    }

                }
    }
#endif





/****************************************************************************


                                 Utilities


*///*************************************************************************



    template<>
    void CopyInSigmaData<1>( DOF_SCRATCH_DATA& output,
                             const DOF_INTERIOR_DATA& input,
                             const int matrix_num ){
        std::array<int,8> matrix_map = {2,6,0,4,3,7,1,5};
        for( int w=0; w<3; w++){
            output.SigmaScratch[0*3+w] = input.L1[ 8][w][matrix_map[matrix_num]];
            output.SigmaScratch[1*3+w] = input.L1[ 6][w][matrix_map[matrix_num]];
            output.SigmaScratch[2*3+w] = input.L1[ 7][w][matrix_map[matrix_num]];
            output.SigmaScratch[3*3+w] = input.L1[ 2][w][matrix_map[matrix_num]];
            output.SigmaScratch[4*3+w] = input.L1[ 0][w][matrix_map[matrix_num]];
            output.SigmaScratch[5*3+w] = input.L1[ 1][w][matrix_map[matrix_num]];
            output.SigmaScratch[6*3+w] = input.L1[ 5][w][matrix_map[matrix_num]];
            output.SigmaScratch[7*3+w] = input.L1[ 3][w][matrix_map[matrix_num]];
            output.SigmaScratch[8*3+w] = input.L1[ 4][w][matrix_map[matrix_num]];
        }
    }

    template<>
    void CopyInSigmaData<2>( DOF_SCRATCH_DATA& output,
                             const DOF_INTERIOR_DATA& input,
                             const int matrix_num ){
        std::array<int,4> matrix_mapA = {2,0,3,1};
        std::array<int,4> matrix_mapB = {10,8,11,9};
        for( int w=0; w<3; w++){
            output.SigmaScratch[ 0*3+w] = input.L2[10][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 1*3+w] = input.L2[ 8][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 2*3+w] = input.L2[ 9][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 3*3+w] = input.L2[ 2][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 4*3+w] = input.L2[ 0][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 5*3+w] = input.L2[ 1][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 6*3+w] = input.L2[ 6][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 7*3+w] = input.L2[ 4][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 8*3+w] = input.L2[ 5][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 9*3+w] = input.L2[10][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[10*3+w] = input.L2[ 8][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[11*3+w] = input.L2[ 9][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[12*3+w] = input.L2[ 2][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[13*3+w] = input.L2[ 0][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[14*3+w] = input.L2[ 1][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[15*3+w] = input.L2[ 6][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[16*3+w] = input.L2[ 4][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[17*3+w] = input.L2[ 5][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[18*3+w] = input.L2[11][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[19*3+w] = input.L2[ 3][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[20*3+w] = input.L2[ 7][w][matrix_mapA[matrix_num]];
        }
    }


    template<> void CopyInSigmaData<3>( DOF_SCRATCH_DATA& output,
                                        const DOF_INTERIOR_DATA& input,
                                        const int matrix_num ){
        std::array<int,2> matrix_mapA = {0,1};
        std::array<int,2> matrix_mapB = {4,5};
        std::array<int,2> matrix_mapC = {8,9};
        std::array<int,2> matrix_mapD = {12,13};
        for( int w=0; w<3; w++){
            output.SigmaScratch[ 0*3+w] = input.L3[10][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 1*3+w] = input.L3[ 8][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 2*3+w] = input.L3[ 9][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 3*3+w] = input.L3[ 2][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 4*3+w] = input.L3[ 0][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 5*3+w] = input.L3[ 1][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 6*3+w] = input.L3[ 6][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 7*3+w] = input.L3[ 4][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 8*3+w] = input.L3[ 5][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 9*3+w] = input.L3[10][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[10*3+w] = input.L3[ 8][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[11*3+w] = input.L3[ 9][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[12*3+w] = input.L3[ 2][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[13*3+w] = input.L3[ 0][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[14*3+w] = input.L3[ 1][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[15*3+w] = input.L3[ 6][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[16*3+w] = input.L3[ 4][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[17*3+w] = input.L3[ 5][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[18*3+w] = input.L3[11][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[19*3+w] = input.L3[ 3][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[20*3+w] = input.L3[ 7][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[21*3+w] = input.L3[10][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[22*3+w] = input.L3[ 8][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[23*3+w] = input.L3[ 9][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[24*3+w] = input.L3[ 2][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[25*3+w] = input.L3[ 0][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[26*3+w] = input.L3[ 1][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[27*3+w] = input.L3[ 6][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[28*3+w] = input.L3[ 4][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[29*3+w] = input.L3[ 5][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[30*3+w] = input.L3[10][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[31*3+w] = input.L3[ 8][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[32*3+w] = input.L3[ 9][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[33*3+w] = input.L3[ 2][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[34*3+w] = input.L3[ 0][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[35*3+w] = input.L3[ 1][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[36*3+w] = input.L3[ 6][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[37*3+w] = input.L3[ 4][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[38*3+w] = input.L3[ 5][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[39*3+w] = input.L3[11][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[40*3+w] = input.L3[ 3][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[41*3+w] = input.L3[ 7][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[42*3+w] = input.L3[14][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[43*3+w] = input.L3[12][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[44*3+w] = input.L3[13][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[45*3+w] = input.L3[14][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[46*3+w] = input.L3[12][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[47*3+w] = input.L3[13][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[48*3+w] = input.L3[15][w][matrix_mapA[matrix_num]] ;
        }
    }

    template<> void CopyInSigmaData<4>( DOF_SCRATCH_DATA& output,
                                        const DOF_INTERIOR_DATA& input,
                                        const int matrix_num ){
        std::array<int,1> matrix_mapA = {2};
        std::array<int,1> matrix_mapB = {6};
        std::array<int,1> matrix_mapC = {10};
        std::array<int,1> matrix_mapD = {14};
        for( int w=0; w<3; w++){
            output.SigmaScratch[ 0*3+w] = input.L4[10][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 1*3+w] = input.L4[ 8][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 2*3+w] = input.L4[ 9][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 3*3+w] = input.L4[ 2][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 4*3+w] = input.L4[ 0][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 5*3+w] = input.L4[ 1][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 6*3+w] = input.L4[ 6][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 7*3+w] = input.L4[ 4][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 8*3+w] = input.L4[ 5][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[ 9*3+w] = input.L4[10][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[10*3+w] = input.L4[ 8][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[11*3+w] = input.L4[ 9][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[12*3+w] = input.L4[ 2][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[13*3+w] = input.L4[ 0][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[14*3+w] = input.L4[ 1][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[15*3+w] = input.L4[ 6][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[16*3+w] = input.L4[ 4][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[17*3+w] = input.L4[ 5][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[18*3+w] = input.L4[11][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[19*3+w] = input.L4[ 3][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[20*3+w] = input.L4[ 7][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[21*3+w] = input.L4[10][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[22*3+w] = input.L4[ 8][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[23*3+w] = input.L4[ 9][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[24*3+w] = input.L4[ 2][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[25*3+w] = input.L4[ 0][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[26*3+w] = input.L4[ 1][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[27*3+w] = input.L4[ 6][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[28*3+w] = input.L4[ 4][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[29*3+w] = input.L4[ 5][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[30*3+w] = input.L4[10][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[31*3+w] = input.L4[ 8][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[32*3+w] = input.L4[ 9][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[33*3+w] = input.L4[ 2][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[34*3+w] = input.L4[ 0][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[35*3+w] = input.L4[ 1][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[36*3+w] = input.L4[ 6][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[37*3+w] = input.L4[ 4][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[38*3+w] = input.L4[ 5][w][matrix_mapD[matrix_num]];
            output.SigmaScratch[39*3+w] = input.L4[11][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[40*3+w] = input.L4[ 3][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[41*3+w] = input.L4[ 7][w][matrix_mapB[matrix_num]];
            output.SigmaScratch[42*3+w] = input.L4[14][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[43*3+w] = input.L4[12][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[44*3+w] = input.L4[13][w][matrix_mapA[matrix_num]];
            output.SigmaScratch[45*3+w] = input.L4[14][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[46*3+w] = input.L4[12][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[47*3+w] = input.L4[13][w][matrix_mapC[matrix_num]];
            output.SigmaScratch[48*3+w] = input.L4[15][w][matrix_mapA[matrix_num]] ;
        }
    }



    template<>
    void CopyOutSigmaData<1>(const DOF_SCRATCH_DATA& input,
                             DOF_INTERIOR_DATA& output,
                             const int matrix_num ){
        //std::array<int,8> matrix_map = {3,7,1,5,4,8,2,6};
        std::array<int,8> matrix_map = {2,6,0,4,3,7,1,5};
        for( int w=0; w<3; w++){
            output.L1[ 0][w][matrix_map[matrix_num]] = input.SigmaScratch[4*3+w];
            output.L1[ 1][w][matrix_map[matrix_num]] = input.SigmaScratch[5*3+w];
            output.L1[ 2][w][matrix_map[matrix_num]] = input.SigmaScratch[3*3+w];
            output.L1[ 3][w][matrix_map[matrix_num]] = input.SigmaScratch[7*3+w];
            output.L1[ 4][w][matrix_map[matrix_num]] = input.SigmaScratch[8*3+w];
            output.L1[ 5][w][matrix_map[matrix_num]] = input.SigmaScratch[6*3+w];
            output.L1[ 6][w][matrix_map[matrix_num]] = input.SigmaScratch[1*3+w];
            output.L1[ 7][w][matrix_map[matrix_num]] = input.SigmaScratch[2*3+w];
            output.L1[ 8][w][matrix_map[matrix_num]] = input.SigmaScratch[0*3+w];
        }
    }

    template<>
    void CopyOutSigmaData<2>(const DOF_SCRATCH_DATA& input,
                             DOF_INTERIOR_DATA& output,
                             const int matrix_num ){
        std::array<int,4> matrix_mapA = {2,0,3,1};
        std::array<int,4> matrix_mapB = {10,8,11,9};
        for( int w=0; w<3; w++){
            output.L2[ 0][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 4*3+w];
            output.L2[ 0][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[13*3+w];
            output.L2[ 1][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 5*3+w];
            output.L2[ 1][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[14*3+w];
            output.L2[ 2][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 3*3+w];
            output.L2[ 2][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[12*3+w];
            output.L2[ 3][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[19*3+w];
            output.L2[ 4][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 7*3+w];
            output.L2[ 4][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[16*3+w];
            output.L2[ 5][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 8*3+w];
            output.L2[ 5][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[17*3+w];
            output.L2[ 6][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 6*3+w];
            output.L2[ 6][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[15*3+w];
            output.L2[ 7][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[20*3+w];
            output.L2[ 8][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 1*3+w];
            output.L2[ 8][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[10*3+w];
            output.L2[ 9][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 2*3+w];
            output.L2[ 9][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[11*3+w];
            output.L2[10][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 0*3+w];
            output.L2[10][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[ 9*3+w];
            output.L2[11][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[18*3+w];
        }
    }

    template<>
    void CopyOutSigmaData<3>(const DOF_SCRATCH_DATA& input,
                             DOF_INTERIOR_DATA& output,
                             const int matrix_num ){
        std::array<int,2> matrix_mapA = {0,1};
        std::array<int,2> matrix_mapB = {4,5};
        std::array<int,2> matrix_mapC = {8,9};
        std::array<int,2> matrix_mapD = {12,13};
        for( int w=0; w<3; w++){
            output.L3[ 0][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 4*3+w];
            output.L3[ 0][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[25*3+w];
            output.L3[ 0][w][matrix_mapC[matrix_num]]  = input.SigmaScratch[13*3+w];
            output.L3[ 0][w][matrix_mapD[matrix_num]]  = input.SigmaScratch[34*3+w];
            output.L3[ 1][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 5*3+w];
            output.L3[ 1][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[26*3+w];
            output.L3[ 1][w][matrix_mapC[matrix_num]]  = input.SigmaScratch[14*3+w];
            output.L3[ 1][w][matrix_mapD[matrix_num]]  = input.SigmaScratch[35*3+w];
            output.L3[ 2][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 3*3+w];
            output.L3[ 2][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[24*3+w];
            output.L3[ 2][w][matrix_mapC[matrix_num]]  = input.SigmaScratch[12*3+w];
            output.L3[ 2][w][matrix_mapD[matrix_num]]  = input.SigmaScratch[33*3+w];
            output.L3[ 3][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[19*3+w];
            output.L3[ 3][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[40*3+w];
            output.L3[ 4][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 7*3+w];
            output.L3[ 4][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[28*3+w];
            output.L3[ 4][w][matrix_mapC[matrix_num]]  = input.SigmaScratch[16*3+w];
            output.L3[ 4][w][matrix_mapD[matrix_num]]  = input.SigmaScratch[37*3+w];
            output.L3[ 5][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 8*3+w];
            output.L3[ 5][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[29*3+w];
            output.L3[ 5][w][matrix_mapC[matrix_num]]  = input.SigmaScratch[17*3+w];
            output.L3[ 5][w][matrix_mapD[matrix_num]]  = input.SigmaScratch[38*3+w];
            output.L3[ 6][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 6*3+w];
            output.L3[ 6][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[27*3+w];
            output.L3[ 6][w][matrix_mapC[matrix_num]]  = input.SigmaScratch[15*3+w];
            output.L3[ 6][w][matrix_mapD[matrix_num]]  = input.SigmaScratch[36*3+w];
            output.L3[ 7][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[20*3+w];
            output.L3[ 7][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[41*3+w];
            output.L3[ 8][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 1*3+w];
            output.L3[ 8][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[22*3+w];
            output.L3[ 8][w][matrix_mapC[matrix_num]]  = input.SigmaScratch[10*3+w];
            output.L3[ 8][w][matrix_mapD[matrix_num]]  = input.SigmaScratch[31*3+w];
            output.L3[ 9][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 2*3+w];
            output.L3[ 9][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[23*3+w];
            output.L3[ 9][w][matrix_mapC[matrix_num]]  = input.SigmaScratch[11*3+w];
            output.L3[ 9][w][matrix_mapD[matrix_num]]  = input.SigmaScratch[32*3+w];
            output.L3[10][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[ 0*3+w];
            output.L3[10][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[21*3+w];
            output.L3[10][w][matrix_mapC[matrix_num]]  = input.SigmaScratch[ 9*3+w];
            output.L3[10][w][matrix_mapD[matrix_num]]  = input.SigmaScratch[30*3+w];
            output.L3[11][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[18*3+w];
            output.L3[11][w][matrix_mapB[matrix_num]]  = input.SigmaScratch[39*3+w];
            output.L3[12][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[43*3+w];
            output.L3[12][w][matrix_mapC[matrix_num]]  = input.SigmaScratch[46*3+w];
            output.L3[13][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[44*3+w];
            output.L3[13][w][matrix_mapC[matrix_num]]  = input.SigmaScratch[47*3+w];
            output.L3[14][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[42*3+w];
            output.L3[14][w][matrix_mapC[matrix_num]]  = input.SigmaScratch[45*3+w];
            output.L3[15][w][matrix_mapA[matrix_num]]  = input.SigmaScratch[48*3+w];
        }
    }

    template<>
    void CopyOutSigmaData<4>(const DOF_SCRATCH_DATA& input,
                             DOF_INTERIOR_DATA& output,
                             const int matrix_num ){
        std::array<int,1> matrix_mapA = {2};
        std::array<int,1> matrix_mapB = {6};
        std::array<int,1> matrix_mapC = {10};
        std::array<int,1> matrix_mapD = {14};
        for( int w=0; w<3; w++){
            output.L4[ 0][w][matrix_mapA[matrix_num]] = input.SigmaScratch[ 4*3+w];
            output.L4[ 0][w][matrix_mapB[matrix_num]] = input.SigmaScratch[25*3+w];
            output.L4[ 0][w][matrix_mapC[matrix_num]] = input.SigmaScratch[13*3+w];
            output.L4[ 0][w][matrix_mapD[matrix_num]] = input.SigmaScratch[34*3+w];
            output.L4[ 1][w][matrix_mapA[matrix_num]] = input.SigmaScratch[ 5*3+w];
            output.L4[ 1][w][matrix_mapB[matrix_num]] = input.SigmaScratch[26*3+w];
            output.L4[ 1][w][matrix_mapC[matrix_num]] = input.SigmaScratch[14*3+w];
            output.L4[ 1][w][matrix_mapD[matrix_num]] = input.SigmaScratch[35*3+w];
            output.L4[ 2][w][matrix_mapA[matrix_num]] = input.SigmaScratch[ 3*3+w];
            output.L4[ 2][w][matrix_mapB[matrix_num]] = input.SigmaScratch[24*3+w];
            output.L4[ 2][w][matrix_mapC[matrix_num]] = input.SigmaScratch[12*3+w];
            output.L4[ 2][w][matrix_mapD[matrix_num]] = input.SigmaScratch[33*3+w];
            output.L4[ 3][w][matrix_mapA[matrix_num]] = input.SigmaScratch[19*3+w];
            output.L4[ 3][w][matrix_mapB[matrix_num]] = input.SigmaScratch[40*3+w];
            output.L4[ 4][w][matrix_mapA[matrix_num]] = input.SigmaScratch[ 7*3+w];
            output.L4[ 4][w][matrix_mapB[matrix_num]] = input.SigmaScratch[28*3+w];
            output.L4[ 4][w][matrix_mapC[matrix_num]] = input.SigmaScratch[16*3+w];
            output.L4[ 4][w][matrix_mapD[matrix_num]] = input.SigmaScratch[37*3+w];
            output.L4[ 5][w][matrix_mapA[matrix_num]] = input.SigmaScratch[ 8*3+w];
            output.L4[ 5][w][matrix_mapB[matrix_num]] = input.SigmaScratch[29*3+w];
            output.L4[ 5][w][matrix_mapC[matrix_num]] = input.SigmaScratch[17*3+w];
            output.L4[ 5][w][matrix_mapD[matrix_num]] = input.SigmaScratch[38*3+w];
            output.L4[ 6][w][matrix_mapA[matrix_num]] = input.SigmaScratch[ 6*3+w];
            output.L4[ 6][w][matrix_mapB[matrix_num]] = input.SigmaScratch[27*3+w];
            output.L4[ 6][w][matrix_mapC[matrix_num]] = input.SigmaScratch[15*3+w];
            output.L4[ 6][w][matrix_mapD[matrix_num]] = input.SigmaScratch[36*3+w];
            output.L4[ 7][w][matrix_mapA[matrix_num]] = input.SigmaScratch[20*3+w];
            output.L4[ 7][w][matrix_mapB[matrix_num]] = input.SigmaScratch[41*3+w];
            output.L4[ 8][w][matrix_mapA[matrix_num]] = input.SigmaScratch[ 1*3+w];
            output.L4[ 8][w][matrix_mapB[matrix_num]] = input.SigmaScratch[22*3+w];
            output.L4[ 8][w][matrix_mapC[matrix_num]] = input.SigmaScratch[10*3+w];
            output.L4[ 8][w][matrix_mapD[matrix_num]] = input.SigmaScratch[31*3+w];
            output.L4[ 9][w][matrix_mapA[matrix_num]] = input.SigmaScratch[ 2*3+w];
            output.L4[ 9][w][matrix_mapB[matrix_num]] = input.SigmaScratch[23*3+w];
            output.L4[ 9][w][matrix_mapC[matrix_num]] = input.SigmaScratch[11*3+w];
            output.L4[ 9][w][matrix_mapD[matrix_num]] = input.SigmaScratch[32*3+w];
            output.L4[10][w][matrix_mapA[matrix_num]] = input.SigmaScratch[ 0*3+w];
            output.L4[10][w][matrix_mapB[matrix_num]] = input.SigmaScratch[21*3+w];
            output.L4[10][w][matrix_mapC[matrix_num]] = input.SigmaScratch[ 9*3+w];
            output.L4[10][w][matrix_mapD[matrix_num]] = input.SigmaScratch[30*3+w];
            output.L4[11][w][matrix_mapA[matrix_num]] = input.SigmaScratch[18*3+w];
            output.L4[11][w][matrix_mapB[matrix_num]] = input.SigmaScratch[39*3+w];
            output.L4[12][w][matrix_mapA[matrix_num]] = input.SigmaScratch[43*3+w];
            output.L4[12][w][matrix_mapC[matrix_num]] = input.SigmaScratch[46*3+w];
            output.L4[13][w][matrix_mapA[matrix_num]] = input.SigmaScratch[44*3+w];
            output.L4[13][w][matrix_mapC[matrix_num]] = input.SigmaScratch[47*3+w];
            output.L4[14][w][matrix_mapA[matrix_num]] = input.SigmaScratch[42*3+w];
            output.L4[14][w][matrix_mapC[matrix_num]] = input.SigmaScratch[45*3+w];
            output.L4[15][w][matrix_mapA[matrix_num]] = input.SigmaScratch[48*3+w];
        }
    }


    template<> void collapse<float,0x8>(float (&data)[16])
    {
        data[ 0]+=data[ 8];data[ 8]=0.f;
        data[ 1]+=data[ 9];data[ 9]=0.f;
        data[ 2]+=data[10];data[10]=0.f;
        data[ 3]+=data[11];data[11]=0.f;
        data[ 4]+=data[12];data[12]=0.f;
        data[ 5]+=data[13];data[13]=0.f;
        data[ 6]+=data[14];data[14]=0.f;
        data[ 7]+=data[15];data[15]=0.f;
    }

    template<> void collapse<float,0x4>(float (&data)[16])
    {
        data[ 0]+=data[ 4];data[ 4]=0.f;
        data[ 1]+=data[ 5];data[ 5]=0.f;
        data[ 2]+=data[ 6];data[ 6]=0.f;
        data[ 3]+=data[ 7];data[ 7]=0.f;
        data[ 8]+=data[12];data[12]=0.f;
        data[ 9]+=data[13];data[13]=0.f;
        data[10]+=data[14];data[14]=0.f;
        data[11]+=data[15];data[15]=0.f;
    }

    template<> void collapse<float,0x2>(float (&data)[16])
    {
        data[ 0]+=data[ 2];data[ 2]=0.f;
        data[ 1]+=data[ 3];data[ 3]=0.f;
        data[ 4]+=data[ 6];data[ 6]=0.f;
        data[ 5]+=data[ 7];data[ 7]=0.f;
        data[ 8]+=data[10];data[10]=0.f;
        data[ 9]+=data[11];data[11]=0.f;
        data[12]+=data[14];data[14]=0.f;
        data[13]+=data[15];data[15]=0.f;
    }

    template<> void collapse<float,0x1>(float (&data)[16])
    {
        data[ 0]+=data[ 1];data[ 1]=0.f;
        data[ 2]+=data[ 3];data[ 3]=0.f;
        data[ 4]+=data[ 5];data[ 5]=0.f;
        data[ 6]+=data[ 7];data[ 7]=0.f;
        data[ 8]+=data[ 9];data[ 9]=0.f;
        data[10]+=data[11];data[11]=0.f;
        data[12]+=data[13];data[13]=0.f;
        data[14]+=data[15];data[15]=0.f;
    }

    void inline CollapseLevel(DOF_INTERIOR_DATA& data, int level){
        switch( level ){
        case 4:
            {
                std::array<int,4> z_stick = {3,7,11,15};
                std::array<int,4> y_stick = {12,13,14,15};
                std::array<int,16> plane = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
                for(int i=0; i<4; i++){
                    collapse<T,0x8>( data.L4[z_stick[i]][0] );
                    collapse<T,0x8>( data.L4[z_stick[i]][1] );
                    collapse<T,0x8>( data.L4[z_stick[i]][2] );
                }
                for(int i=0; i<4; i++){
                    collapse<T,0x4>( data.L4[y_stick[i]][0] );
                    collapse<T,0x4>( data.L4[y_stick[i]][1] );
                    collapse<T,0x4>( data.L4[y_stick[i]][2] );
                }
                for(int i=0; i<16; i++){
                    collapse<T,0x1>( data.L4[plane[i]][0] );
                    collapse<T,0x1>( data.L4[plane[i]][1] );
                    collapse<T,0x1>( data.L4[plane[i]][2] );
                }
            }
            break;
        case 3:
            {
                std::array<int,4> z_stick = {3,7,11,15};
                std::array<int,4> y_stick = {12,13,14,15};
                std::array<int,16> plane = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
                for(int i=0; i<4; i++){
                    collapse<T,0x8>( data.L3[z_stick[i]][0] );
                    collapse<T,0x8>( data.L3[z_stick[i]][1] );
                    collapse<T,0x8>( data.L3[z_stick[i]][2] );
                }
                for(int i=0; i<4; i++){
                    collapse<T,0x4>( data.L3[y_stick[i]][0] );
                    collapse<T,0x4>( data.L3[y_stick[i]][1] );
                    collapse<T,0x4>( data.L3[y_stick[i]][2] );
                }
                for(int i=0; i<16; i++){
                    collapse<T,0x2>( data.L3[plane[i]][0] );
                    collapse<T,0x2>( data.L3[plane[i]][1] );
                    collapse<T,0x2>( data.L3[plane[i]][2] );
                }
            }
            break;
        case 2:
            {
                std::array<int,3> stick = {3,7,11};
                for(int i=0; i<3; i++){
                    collapse<T,0x8>( data.L2[stick[i]][0] );
                    collapse<T,0x8>( data.L2[stick[i]][1] );
                    collapse<T,0x8>( data.L2[stick[i]][2] );
                }
                std::array<int,12> plane = {0,1,2,3,4,5,6,7,8,9,10,11};
                for(int i=0; i<12; i++){
                    collapse<T,0x4>( data.L2[plane[i]][0] );
                    collapse<T,0x4>( data.L2[plane[i]][1] );
                    collapse<T,0x4>( data.L2[plane[i]][2] );
                }
            }
            break;
        case 1:
            {
                for(int i=0; i<9; i++){
                    collapse<T,0x8>( data.L1[i][0] );
                    collapse<T,0x8>( data.L1[i][1] );
                    collapse<T,0x8>( data.L1[i][2] );
                }
            }
            break;
        case 0:
            break;
        }
    }


    template<> void prune<float,0x8>(float (&data)[16])
    {
        data[ 8]=0.f;
        data[ 9]=0.f;
        data[10]=0.f;
        data[11]=0.f;
        data[12]=0.f;
        data[13]=0.f;
        data[14]=0.f;
        data[15]=0.f;
    }

    template<> void prune<float,0x4>(float (&data)[16])
    {
        data[ 4]=0.f;
        data[ 5]=0.f;
        data[ 6]=0.f;
        data[ 7]=0.f;
        data[12]=0.f;
        data[13]=0.f;
        data[14]=0.f;
        data[15]=0.f;
    }

    template<> void prune<float,0x2>(float (&data)[16])
    {
        data[ 2]=0.f;
        data[ 3]=0.f;
        data[ 6]=0.f;
        data[ 7]=0.f;
        data[10]=0.f;
        data[11]=0.f;
        data[14]=0.f;
        data[15]=0.f;
    }

    template<> void prune<float,0x1>(float (&data)[16])
    {
        data[ 1]=0.f;
        data[ 3]=0.f;
        data[ 5]=0.f;
        data[ 7]=0.f;
        data[ 9]=0.f;
        data[11]=0.f;
        data[13]=0.f;
        data[15]=0.f;
    }

    void inline PruneLevel(DOF_INTERIOR_DATA& data, int level){
        switch( level ){
        case 4:
            {
                std::array<int,4> z_stick = {3,7,11,15};
                for(int i=0; i<4; i++){
                    prune<T,0x8>( data.L4[z_stick[i]][0] );
                    prune<T,0x8>( data.L4[z_stick[i]][1] );
                    prune<T,0x8>( data.L4[z_stick[i]][2] );
                }
                std::array<int,4> y_stick = {12,13,14,15};
                for(int i=0; i<4; i++){
                    prune<T,0x4>( data.L4[y_stick[i]][0] );
                    prune<T,0x4>( data.L4[y_stick[i]][1] );
                    prune<T,0x4>( data.L4[y_stick[i]][2] );
                }
                std::array<int,16> plane = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
                for(int i=0; i<16; i++){
                    prune<T,0x1>( data.L4[plane[i]][0] );
                    prune<T,0x1>( data.L4[plane[i]][1] );
                    prune<T,0x1>( data.L4[plane[i]][2] );
                }
            }
            break;
        case 3:
            {
                std::array<int,4> z_stick = {3,7,11,15};
                for(int i=0; i<4; i++){
                    prune<T,0x8>( data.L3[z_stick[i]][0] );
                    prune<T,0x8>( data.L3[z_stick[i]][1] );
                    prune<T,0x8>( data.L3[z_stick[i]][2] );
                }
                std::array<int,4> y_stick = {12,13,14,15};
                for(int i=0; i<4; i++){
                    prune<T,0x4>( data.L3[y_stick[i]][0] );
                    prune<T,0x4>( data.L3[y_stick[i]][1] );
                    prune<T,0x4>( data.L3[y_stick[i]][2] );
                }
                std::array<int,16> plane = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
                for(int i=0; i<16; i++){
                    prune<T,0x2>( data.L3[plane[i]][0] );
                    prune<T,0x2>( data.L3[plane[i]][1] );
                    prune<T,0x2>( data.L3[plane[i]][2] );
                }
            }
            break;
        case 2:
            {
                std::array<int,3> stick = {3,7,11};
                for(int i=0; i<3; i++){
                    prune<T,0x8>( data.L2[stick[i]][0] );
                    prune<T,0x8>( data.L2[stick[i]][1] );
                    prune<T,0x8>( data.L2[stick[i]][2] );
                }
                std::array<int,12> plane = {0,1,2,3,4,5,6,7,8,9,10,11};
                for(int i=0; i<12; i++){
                    prune<T,0x4>( data.L2[plane[i]][0] );
                    prune<T,0x4>( data.L2[plane[i]][1] );
                    prune<T,0x4>( data.L2[plane[i]][2] );
                }
            }
            break;
        case 1:
            {
                for(int i=0; i<9; i++){
                    prune<T,0x8>( data.L1[i][0] );
                    prune<T,0x8>( data.L1[i][1] );
                    prune<T,0x8>( data.L1[i][2] );
                }
            }
            break;
        case 0:
            break;
        }
    }

    template<> void distribute<float,0x8>(float (&data)[16])
    {
        data[ 8]=data[ 0];
        data[ 9]=data[ 1];
        data[10]=data[ 2];
        data[11]=data[ 3];
        data[12]=data[ 4];
        data[13]=data[ 5];
        data[14]=data[ 6];
        data[15]=data[ 7];
    }

    template<> void distribute<float,0x4>(float (&data)[16])
    {
        data[ 4]=data[ 0];
        data[ 5]=data[ 1];
        data[ 6]=data[ 2];
        data[ 7]=data[ 3];
        data[12]=data[ 8];
        data[13]=data[ 9];
        data[14]=data[10];
        data[15]=data[11];
    }

    template<> void distribute<float,0x2>(float (&data)[16])
    {
        data[ 2]=data[ 0];
        data[ 3]=data[ 1];
        data[ 6]=data[ 4];
        data[ 7]=data[ 5];
        data[10]=data[ 8];
        data[11]=data[ 9];
        data[14]=data[12];
        data[15]=data[13];
    }

    template<> void distribute<float,0x1>(float (&data)[16])
    {
        data[ 1]=data[ 0];
        data[ 3]=data[ 2];
        data[ 5]=data[ 4];
        data[ 7]=data[ 6];
        data[ 9]=data[ 8];
        data[11]=data[10];
        data[13]=data[12];
        data[15]=data[14];
    }

    void inline DistributeLevel(DOF_INTERIOR_DATA& data, int level){
        switch( level ){
        case 4:
            {
                std::array<int,4> z_stick = {3,7,11,15};
                for(int i=0; i<4; i++){
                    distribute<T,0x8>( data.L4[z_stick[i]][0] );
                    distribute<T,0x8>( data.L4[z_stick[i]][1] );
                    distribute<T,0x8>( data.L4[z_stick[i]][2] );
                }
                std::array<int,4> y_stick = {12,13,14,15};
                for(int i=0; i<4; i++){
                    distribute<T,0x4>( data.L4[y_stick[i]][0] );
                    distribute<T,0x4>( data.L4[y_stick[i]][1] );
                    distribute<T,0x4>( data.L4[y_stick[i]][2] );
                }
                std::array<int,16> plane = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
                for(int i=0; i<16; i++){
                    distribute<T,0x1>( data.L4[plane[i]][0] );
                    distribute<T,0x1>( data.L4[plane[i]][1] );
                    distribute<T,0x1>( data.L4[plane[i]][2] );
                }
            }
            break;
        case 3:
            {
                std::array<int,4> z_stick = {3,7,11,15};
                for(int i=0; i<4; i++){
                    distribute<T,0x8>( data.L3[z_stick[i]][0] );
                    distribute<T,0x8>( data.L3[z_stick[i]][1] );
                    distribute<T,0x8>( data.L3[z_stick[i]][2] );
                }
                std::array<int,4> y_stick = {12,13,14,15};
                for(int i=0; i<4; i++){
                    distribute<T,0x4>( data.L3[y_stick[i]][0] );
                    distribute<T,0x4>( data.L3[y_stick[i]][1] );
                    distribute<T,0x4>( data.L3[y_stick[i]][2] );
                }
                std::array<int,16> plane = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
                for(int i=0; i<16; i++){
                    distribute<T,0x2>( data.L3[plane[i]][0] );
                    distribute<T,0x2>( data.L3[plane[i]][1] );
                    distribute<T,0x2>( data.L3[plane[i]][2] );
                }
            }
            break;
        case 2:
            {
                std::array<int,3> stick = {3,7,11};
                for(int i=0; i<3; i++){
                    distribute<T,0x8>( data.L2[stick[i]][0] );
                    distribute<T,0x8>( data.L2[stick[i]][1] );
                    distribute<T,0x8>( data.L2[stick[i]][2] );
                }
                std::array<int,12> plane = {0,1,2,3,4,5,6,7,8,9,10,11};
                for(int i=0; i<12; i++){
                    distribute<T,0x4>( data.L2[plane[i]][0] );
                    distribute<T,0x4>( data.L2[plane[i]][1] );
                    distribute<T,0x4>( data.L2[plane[i]][2] );
                }
            }
            break;
        case 1:
            {
                for(int i=0; i<9; i++){
                    distribute<T,0x8>( data.L1[i][0] );
                    distribute<T,0x8>( data.L1[i][1] );
                    distribute<T,0x8>( data.L1[i][2] );
                }
            }
            break;
        case 0:
            break;
        }
    }

    template<class T> void negate(T (&data)[16]){
        data[ 0] = -data[ 0];
        data[ 1] = -data[ 1];
        data[ 2] = -data[ 2];
        data[ 3] = -data[ 3];
        data[ 4] = -data[ 4];
        data[ 5] = -data[ 5];
        data[ 6] = -data[ 6];
        data[ 7] = -data[ 7];
        data[ 8] = -data[ 8];
        data[ 9] = -data[ 9];
        data[10] = -data[10];
        data[11] = -data[11];
        data[12] = -data[12];
        data[13] = -data[13];
        data[14] = -data[14];
        data[15] = -data[15];
    }

    template<class T> void plus_equals(T (&data)[16], const T (&input)[16]){
        data[ 0] += input[ 0];
        data[ 1] += input[ 1];
        data[ 2] += input[ 2];
        data[ 3] += input[ 3];
        data[ 4] += input[ 4];
        data[ 5] += input[ 5];
        data[ 6] += input[ 6];
        data[ 7] += input[ 7];
        data[ 8] += input[ 8];
        data[ 9] += input[ 9];
        data[10] += input[10];
        data[11] += input[11];
        data[12] += input[12];
        data[13] += input[13];
        data[14] += input[14];
        data[15] += input[15];
    }


    void inline CopyDown( const DOF_INTERIOR_DATA& in, DOF_INTERIOR_DATA& out, int level){
        switch( level ){
        case 4:
            memcpy( out.L4, in.L4, sizeof( in.L4 ) );
        case 3:
            memcpy( out.L3, in.L3, sizeof( in.L3 ) );
        case 2:
            memcpy( out.L2, in.L2, sizeof( in.L2 ) );
        case 1:
            memcpy( out.L1, in.L1, sizeof( in.L1 ) );
        case 0:
            memcpy( out.L0, in.L0, sizeof( in.L0 ) );
        }
    }

    void inline CopyUp( const DOF_INTERIOR_DATA& in, DOF_INTERIOR_DATA& out, int level){
        switch( level ){
        case 0:
            memcpy( out.L0, in.L0, sizeof( in.L0 ) );
        case 1:
            memcpy( out.L1, in.L1, sizeof( in.L1 ) );
        case 2:
            memcpy( out.L2, in.L2, sizeof( in.L2 ) );
        case 3:
            memcpy( out.L3, in.L3, sizeof( in.L3 ) );
        case 4:
            memcpy( out.L4, in.L4, sizeof( in.L4 ) );
        }
    }

    void inline CopyAt( const DOF_INTERIOR_DATA& in, DOF_INTERIOR_DATA& out, int level){
        switch( level ){
        case 0:
            memcpy( out.L0, in.L0, sizeof( in.L0 ) );
            break;
        case 1:
            memcpy( out.L1, in.L1, sizeof( in.L1 ) );
            break;
        case 2:
            memcpy( out.L2, in.L2, sizeof( in.L2 ) );
            break;
        case 3:
            memcpy( out.L3, in.L3, sizeof( in.L3 ) );
            break;
        case 4:
            memcpy( out.L4, in.L4, sizeof( in.L4 ) );
            break;
        }
    }

    void inline NegateLevel(DOF_INTERIOR_DATA& data, int level, bool copydown){
        int term = copydown ? DOF_INTERIOR_DATA::LevelSizeAndBelow(level)/4 : DOF_INTERIOR_DATA::LevelSize(level)/4;
        float *data_ptr = copydown ? data.Level(0) : data.Level(level);
        for(int i=0; i<term; i+=VECTOR_WIDTH){
            negate<T>( reinterpret_cast<T_DATA_SCALAR>(*(data_ptr+i)) );
        }
    }

    void inline PlusEqualsLevel(DOF_INTERIOR_DATA& data, const DOF_INTERIOR_DATA& other, int level, bool copydown){
        int term = copydown ? DOF_INTERIOR_DATA::LevelSizeAndBelow(level)/4 : DOF_INTERIOR_DATA::LevelSize(level)/4;
        float *data_ptr = copydown ? data.Level(0) : data.Level(level);
        float *other_ptr = copydown ? other.Level(0) : other.Level(level);
        for(int i=0; i<term; i+=VECTOR_WIDTH){
            plus_equals<T>( reinterpret_cast<T_DATA_SCALAR>(*(data_ptr+i)),
                            reinterpret_cast<T_DATA_SCALAR>(*(other_ptr+i)) );
        }
    }

}
