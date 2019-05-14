//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class Tw,class T_DATA=void, class I_DATA=void> 
    void Transpose(T_DATA (&At)[9], const T_DATA (&A)[9]);

template<class Tw,class T_DATA=void, class I_DATA=void>
    void Build_M(T_DATA (&M)[9], const T_DATA (&dPdF)[12], const T_DATA (&H)[3][8], const T_DATA (&scale), int i, int k );

template<class Tw,class T_DATA=void, class I_DATA=void> 
    void Compute_Per_Vertex_Cell_Matrix(const int v1, const int v2,
                                        const T_DATA (&one_over_h),
                                        const T_DATA (&cell_volume),
                                        const T_DATA (&Weights)[3],
                                        const T_DATA (&U)[9],
                                        const T_DATA (&V)[9],
                                        const T_DATA (&dPdF)[12],
                                        
                                        T_DATA (&ij_matrix)[9]);
