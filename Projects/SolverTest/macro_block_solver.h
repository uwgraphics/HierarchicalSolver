#ifndef __MACRO_BLOCK_SOLVER_H__
#define __MACRO_BLOCK_SOLVER_H__

#include "Macroblock_Utilities/BlockedCSCSymmetricMatrix3.h"
#include "Macroblock_Utilities/BlockedCSCSymmetricMatrix3_Reference.h"
#include "Macroblock_Utilities/BlockedCSCSymmetricMatrix3_Utilities.h"
#include "Macroblock_Utilities/Sparse_Cholesky.h"
#include "Macroblock_Utilities/Sparse_Blocked_Forwards_Substitution.h"
#include "Macroblock_Utilities/Sparse_Blocked_Backwards_Substitution.h"

#include "Macroblock_Utilities/BlockedCSRMatrix3.h"
#include "Macroblock_Utilities/BlockedCSRMatrix3_Reference.h"
#include "Macroblock_Utilities/BlockedCSRMatrix3_Utilities.h"
#include "Macroblock_Utilities/Sparse_Blocked_CSR_Multiply.h"

#include "Macroblock_Utilities/BlockedMatrixNXN.h"
#include "Macroblock_Utilities/BlockedMatrixNXN_Utilities.h"
#include "Macroblock_Utilities/BlockedMatrixNXN_Cholesky.h"

#include "immintrin.h"

#include <tuple>

namespace MACROBLOCK{

    enum{ VECTOR_WIDTH=16, INSTRUCTION_WIDTH=16 };
    typedef __m512 INSTRUCTION_TYPE;
    typedef float T;
    typedef T T_DATA[VECTOR_WIDTH];
    typedef int I_DATA[VECTOR_WIDTH];
    typedef const float (&C_T_DATA_SCALAR)[VECTOR_WIDTH];    typedef float (&T_DATA_SCALAR)[VECTOR_WIDTH];


    struct MATRIX_PARAMETERS {
        float lambda [64][16];
        float mu [64][16];

        float dx;
        int seed;
        int max_level;
    };


    struct MATRIX_DATA {
        // Subdomain Factorization Matrix data (BlockedCSCSymmetricMatrix3)
        // [27 dof in subdomain matrices][6 entries in symmetric 3x3 matrix][16 subdomains / vector width]
        // [exactly 188 entries in the SPARSE off diagonal part][][]
        float L0_diagonal_entries    [27][ 6][16];
        float L0_offdiagonal_entries [188][9][16];

        // Sigma Factorization Matrix data (BlockedMatrixNXN)
        // each matrix block is 4x4 - need the whole block 16 entries
        // Level 1
        float L1_diagonal_entries    [8][  8][16]; // 7 = 3x3x1 x3/4   // [8 blocks][8 matrices][vector width]
        float L1_offdiagonal_entries [8][ 28][16]; // 21 = 7x6/2
        // Level 2
        float L2_diagonal_entries    [4][ 16][16]; // 16 = 3x1x7 x3/4
        float L2_offdiagonal_entries [4][120][16]; // 120 = 16*15/2
        // Level 3
        float L3_diagonal_entries    [2][ 40][16]; // 37 = 1x7x7 x3/4
        float L3_offdiagonal_entries [2][780][16]; // 703 why? 780 = 40*39/2
        // Level 4
        float L4_diagonal_entries    [1][ 40][16]; // 37 = 1x7x7 x3/4
        float L4_offdiagonal_entries [1][780][16]; // 703

        /* float A0_diagonal_entries     [27][6][16]; */
        /* float A1_diagonal_entries     [ 9][6][16]; */
        /* float A2_diagonal_entries     [12][6][16]; */
        /* float A3_diagonal_entries     [16][6][16]; */
        /* float A4_diagonal_entries     [16][6][16]; */
        /* float A5_diagonal_entries     [45][6][16]; */

        /* float A0_offdiagonal_entries [158][9][16]; */
        /* float A1_offdiagonal_entries [ 20][9][16]; */
        /* float A2_offdiagonal_entries [ 29][9][16]; */
        /* float A3_offdiagonal_entries [ 42][9][16]; */
        /* float A4_offdiagonal_entries [ 42][9][16]; */
        /* float A5_offdiagonal_entries [153][9][16]; */

        // Transfer Matricies (BlockedCSRMatrix3)
        float L1_to_L0               [ 49][9][16]; //
        float L2_to_L1L0             [ 70][9][16];
        float L3_to_L2L1L0           [100][9][16];
        float L4_to_L3L2L1L0         [100][9][16];
        // this is not used
        //    float L5_to_L4L3L2L1L0       [273][9][16];  // Everything should add up to 1036


        // Accessor Methods
        constexpr float* SigmaDiagonal(int level, int matrix) const {
            return level == 4 ?
                (float*)L4_diagonal_entries[matrix]:
                level == 3 ?
                (float*)L3_diagonal_entries[matrix]:
                level == 2 ?
                (float*)L2_diagonal_entries[matrix]:
                level == 1 ?
                (float*)L1_diagonal_entries[matrix]:
                0;
        }

        constexpr float* SigmaOffDiagonal(int level, int matrix) const {
            return level == 4 ?
                (float*)L4_offdiagonal_entries[matrix]:
                level == 3 ?
                (float*)L3_offdiagonal_entries[matrix]:
                level == 2 ?
                (float*)L2_offdiagonal_entries[matrix]:
                level == 1 ?
                (float*)L1_offdiagonal_entries[matrix]:
                0;
        }

        static constexpr int SigmaDiagonalSize(int level){
            return level == 4 ?
                sizeof(L4_diagonal_entries[0]):
                level == 3 ?
                sizeof(L3_diagonal_entries[0]):
                level == 2 ?
                sizeof(L2_diagonal_entries[0]):
                level == 1 ?
                sizeof(L1_diagonal_entries[0]):
                0;
        }

        static constexpr int SigmaOffDiagonalSize(int level){
            return level == 4 ?
                sizeof(L4_offdiagonal_entries[0]):
                level == 3 ?
                sizeof(L3_offdiagonal_entries[0]):
                level == 2 ?
                sizeof(L2_offdiagonal_entries[0]):
                level == 1 ?
                sizeof(L1_offdiagonal_entries[0]):
                0;
        }

        static constexpr int TransferSize(int level){
            return level == 4 ?
                sizeof( L4_to_L3L2L1L0 ) :
                level == 3 ?
                sizeof( L3_to_L2L1L0 ) :
                level == 2 ?
                sizeof( L2_to_L1L0 ):
                level == 1 ?
                sizeof( L1_to_L0 ):
                0;
        }

        constexpr float* TransferPointer(int level) const {
            return level == 4 ?
                (float*)L4_to_L3L2L1L0 :
                level == 3 ?
                (float*)L3_to_L2L1L0 :
                level == 2 ?
                (float*)L2_to_L1L0:
                level == 1 ?
                (float*)L1_to_L0:
                0;
        }

        /* constexpr float* A_Offdiagonals() const { */
        /*     return (float*)( A0_offdiagonal_entries  ); */
        /* } */

        /* constexpr float* A_Diagonals() const { */
        /*     return (float*)( A0_diagonal_entries ); */
        /* } */

    };


    struct MATRIX_STRUCTURE {
        // Subdomain Factorization Matrix structure
        enum {L0_n=27};
      int L0_offsets               [ 28]; // 3x3x3 +1
      int L0_rows                  [188];

      // Transfer Matricies structure
      // Level 1 to Level 0
      int L1_to_L0_offsets           [10]; // 3x3x1 +1 = TransferRows(1)+1
      int L1_to_L0_column            [49];
      // Level 2 to Levels 1,0
      int L2_to_L1L0_offsets         [13]; // 3x1x4 +1
      int L2_to_L1L0_column          [70];
      // Level 3 to Levels 2,1,0
      int L3_to_L2L1L0_offsets       [17]; // 1x4x4 +1
      int L3_to_L2L1L0_column        [100];
      // Level 4 to Levels 3,2,1,0
      int L4_to_L3L2L1L0_offsets     [17]; // 1x4x4 +1
      int L4_to_L3L2L1L0_column      [100];
      // Level 5 to Levels 4,3,2,1,0
      int L5_to_L4L3L2L1L0_offsets   [46]; // 5x1x4 + 5x5x1 +1
      int L5_to_L4L3L2L1L0_column    [273];

        int cell_index            [16][8][8];
        int cell_subdomain        [16][8][8];
        int dof_index              [5][5][5];
        int element_coeff    [4][4][4][8][8];

        // Accessor Functions
        static constexpr int SigmaN(int level){
            return level == 4 ?
                40:
                level == 3 ?
                40:
                level == 2 ?
                16:
                level == 1 ?
                8:
                0;
        };

        static constexpr int SigmaCount(int level){
            return level == 4 ?
                1:
                level == 3 ?
                2:
                level == 2 ?
                4:
                level == 1 ?
                8:
                0;
        };

        static constexpr int TransferRows(int level){
            return level == 5 ?
                45:
                level == 4 ?
                16:
                level == 3 ?
                16:
                level == 2 ?
                12:
                level == 1 ?
              9:
                0;
        }

        static constexpr int TransferColumns(int level){
            return level == 5 ?
              80: //4x4x5
              level == 4 ?
              64: //4x4x4
              level == 3 ?
              48: //3x4x4
              level == 2 ?
              36: //3x3x4
              level == 1 ?
              27: // 3x3x3
              0;
        }

        constexpr int* TransferOffsets(int level) const {
            return level == 5 ?
                (int*)L5_to_L4L3L2L1L0_offsets :
                level == 4 ?
                (int*)L4_to_L3L2L1L0_offsets :
                level == 3 ?
                (int*)L3_to_L2L1L0_offsets :
                level == 2 ?
                (int*)L2_to_L1L0_offsets:
                level == 1 ?
                (int*)L1_to_L0_offsets:
                0;
        }

        constexpr int* TransferColumn(int level) const {
            return level == 5 ?
                (int*)L5_to_L4L3L2L1L0_column :
                level == 4 ?
                (int*)L4_to_L3L2L1L0_column :
                level == 3 ?
                (int*)L3_to_L2L1L0_column :
                level == 2 ?
                (int*)L2_to_L1L0_column:
                level == 1 ?
                (int*)L1_to_L0_column:
                0;
        }

        static constexpr int TransferOffsetsSize(int level){
            return level == 5 ?
                sizeof(L5_to_L4L3L2L1L0_offsets) :
                level == 4 ?
                sizeof(L4_to_L3L2L1L0_offsets):
                level == 3 ?
                sizeof(L3_to_L2L1L0_offsets):
                level == 2 ?
                sizeof(L2_to_L1L0_offsets):
                level == 1 ?
                sizeof(L1_to_L0_offsets):
                0;
        }

        static constexpr int TransferColumnSize(int level){
            return level == 5 ?
                sizeof(L5_to_L4L3L2L1L0_column):
                level == 4 ?
                sizeof(L4_to_L3L2L1L0_column):
                level == 3 ?
                sizeof(L3_to_L2L1L0_column):
                level == 2 ?
                sizeof(L2_to_L1L0_column):
                level == 1 ?
                sizeof(L1_to_L0_column):
                0;
        }

    };


    struct DOF_INTERIOR_DATA {
        // Interior + plus some extra
      float L0       [27][3][16]; // 3x3x3
      float L1       [ 9][3][16]; // 3x3x1
      float L2       [12][3][16]; // 3x1x4
      float L3       [16][3][16]; // 1x4x4
        float L4       [16][3][16]; // 1x4x4
        float L5       [45][3][16]; // 5x1x4 + 5x5x1

        // Accessors
        constexpr float* Level(int level) const {
            return level == 4 ?
                (float*)L4:
                level == 3 ?
                (float*)L3:
                level == 2 ?
                (float*)L2:
                level == 1 ?
                (float*)L1:
                level == 0 ?
                (float*)L0:
                0;
        }

        inline float* Level_(int level) const {
            return level == 4 ?
                (float*)L4:
                level == 3 ?
                (float*)L3:
                level == 2 ?
                (float*)L2:
                level == 1 ?
                (float*)L1:
                level == 0 ?
                (float*)L0:
                0;
        }

        static constexpr int LevelSize(int level) {
            return level == 4 ?
                sizeof(L4):
                level == 3 ?
                sizeof(L3):
                level == 2 ?
                sizeof(L2):
                level == 1 ?
                sizeof(L1):
                level == 0 ?
                sizeof(L0):
                0;
        }

        static constexpr int LevelSizeAndBelow(int level) {
            return level == 4 ?
                sizeof(L4)+sizeof(L3)+sizeof(L2)+sizeof(L1)+sizeof(L0):
                level == 3 ?
                sizeof(L3)+sizeof(L2)+sizeof(L1)+sizeof(L0):
                level == 2 ?
                sizeof(L2)+sizeof(L1)+sizeof(L0):
                level == 1 ?
                sizeof(L1)+sizeof(L0):
                level == 0 ?
                sizeof(L0):
                0;
        }

    };

    struct DOF_SCRATCH_DATA {
        DOF_INTERIOR_DATA L0;
        DOF_INTERIOR_DATA L1;
        DOF_INTERIOR_DATA L2;
        DOF_INTERIOR_DATA L3;
        DOF_INTERIOR_DATA L4;
        float mask[16];
        float SigmaScratch[160];

        DOF_INTERIOR_DATA& Level(int level) {
            return *( level == 4 ?
                      &L4:
                      level == 3 ?
                      &L3:
                      level == 2 ?
                      &L2:
                      level == 1 ?
                      &L1:
                      level == 0 ?
                      &L0:
                      0);
        }

    };

    /****************************************************************************


                               Utilities


    *///*************************************************************************

 template<int d> static void CopyInSigmaData( DOF_SCRATCH_DATA& output, const DOF_INTERIOR_DATA& input, const int matrix_num );
     template<int d> static void CopyOutSigmaData(const DOF_SCRATCH_DATA& input, DOF_INTERIOR_DATA& output, const int matrix_num );
    void CopyDown( const DOF_INTERIOR_DATA& in, DOF_INTERIOR_DATA& out, int level);
    void CopyUp( const DOF_INTERIOR_DATA& in, DOF_INTERIOR_DATA& out, int level);
    void CopyAt( const DOF_INTERIOR_DATA& in, DOF_INTERIOR_DATA& out, int level);
    template<class T,int mask> void collapse(T (&data)[16]);
    template<class T,int mask> void prune(T (&data)[16]);
    template<class T,int mask> void distribute(T (&data)[16]);
    template<class T> void negate(T (&data)[16]);
    template<class T> void plus_equals(T (&data)[16], const T (&input)[16]);

    /****************************************************************************


                         Solve Level Routines


    *///*************************************************************************

    template<int level, bool in_place>
        void Solve_Level(const MATRIX_STRUCTURE& matrix_structure, const MATRIX_DATA& matrix_data, DOF_INTERIOR_DATA& y, DOF_INTERIOR_DATA& x, DOF_SCRATCH_DATA& scratch);

    void Initialize_Matrix_Structure( MATRIX_STRUCTURE& matrix_structure );

    void Update_Position_Based_State_And_Add_Force( const MATRIX_STRUCTURE& matrix_structure, const MATRIX_PARAMETERS& parameters, MATRIX_DATA& matrix_data, DOF_INTERIOR_DATA& X, DOF_INTERIOR_DATA& forces );

}

#endif
