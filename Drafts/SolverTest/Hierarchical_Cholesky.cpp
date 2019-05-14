//#####################################################################
// Copyright (c) 2018, Eftychios Sifakis, Qisi Wang
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include "Hierarchical_Cholesky.h"
#include "./Macroblock_Utilities/Sparse_Cholesky.h"
#include "./DenseBlockOperations.h"
#include "./SymmetricMatrix.h"
#include <iostream>
#include <algorithm>
#include <mkl.h>
//#include <typeinfo>


#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
namespace Macroblock_Utilities
{
  //#####################################################################
  // Function Hierarchical_Cholesky::Initialize
  //#####################################################################
  template <class T>
  void Hierarchical_Cholesky<T>::Initialize()
  {
    sbd_topology.Initialize_Index_3x3x3();
    sbd_topology.Initialize_Coordinates();
    sbd_topology.Initialize_Matrix_Topology();

    depth=0;
    {
      const int number_of_sbd = (domain_size[0]+1)*(domain_size[1]+1)*(domain_size[2]+1)/(sbd_topology.subdomain_size[0]+1)/(sbd_topology.subdomain_size[1]+1)/(sbd_topology.subdomain_size[2]+1);
      for (int r=number_of_sbd-1; r; r>>=1, depth++);
    }

    agg_topologies.resize(depth);

    for (int l=1; l<=depth; l++) {
      if (l==1)
        agg_topologies[l-1]=Aggregate_Topology::Aggregate(sbd_topology,depth);
      else
        agg_topologies[l-1]=Aggregate_Topology::Aggregate(*(agg_topologies[l-2]), depth);
      agg_topologies[l-1]->Initialize_Coordinates();

    }
    sbd_csc_topology.n = sbd_topology.Nsubdomain;
    sbd_csc_topology.offsets.resize(sbd_csc_topology.n+1);
    sbd_csc_topology.offsets[0] = 0;
    int count = 0;
    for (int j=0; j<sbd_csc_topology.n; j++) {
      for (int ii=sbd_topology.L.offsets[j]; ii<sbd_topology.L.offsets[j+1]; ii++) {
        int i = sbd_topology.L.row[ii];
        if (i<sbd_csc_topology.n) {
          sbd_csc_topology.row.push_back(i);
          count++;
        }
        sbd_csc_topology.offsets[j+1] = count;
      }}
  }

  //#####################################################################
  // static Function Hierarchical_Cholesky::Extract_Schur
  //#####################################################################
  template<class T>
  void Hierarchical_Cholesky<T>::Extract_Schur(const int level,  const int rank) {
    //PhysBAM::LOG::SCOPE scope("Extract_Schur");
    if (level ==0) {
      BlockedCSCSymmetricMatrix3<T> schur{sbd_csc_topology.n,sbd_schur[rank].L_diagonal,sbd_schur[rank].L_lower,&sbd_csc_topology.offsets[0],&sbd_csc_topology.row[0]};
      BlockedCSCSymmetricMatrix3<T> workspace{sbd_topology.L.n,sbd_workspace.L_diagonal,sbd_workspace.L_lower,&sbd_topology.L.offsets[0],&sbd_topology.L.row[0]};

      for (int j=0; j<schur.n; j++) {
        {
          auto Mat_S = schur.D(j);
          auto Mat_W = workspace.D(j);
          for (int v=0; v<6; v++)
            Mat_S[v] = Mat_W[v];
        }
        for (int ii=schur.offsets[j], kk=workspace.offsets[j]; ii<schur.offsets[j+1]; ii++, kk++) {
          assert(schur.row[ii]==workspace.row[kk]);
          auto Mat_S = schur.L(ii);
          auto Mat_W = workspace.L(kk);
          for (int v=0; v<d*d; v++)
            Mat_S[v] = Mat_W[v];
        }}
    } else {
      SymmetricMatrix<T> schur{agg_topologies[level-1]->Ninterface*3, agg_schur[level-1][rank]};
      Aggregate_Data<T>& workspace=agg_workspace[level-1];
      for (int j=0; j<schur.n; j++)
        for (int i=j; i<schur.n; i++)
          schur(i,j) = workspace(i,j);
    }
  }

  //#####################################################################
  // static Function Hierarchical_Cholesky::Accumulate_to_Next_Level
  //#####################################################################
  template<class T>
  void Hierarchical_Cholesky<T>::Accumulate_to_Next_Level(const int level,  const int child_index) {
    //PhysBAM::LOG::SCOPE scope("Accumulate_to_Next_Level");
    Aggregate_Data<T>& parent=agg_workspace[level];
    const std::vector<std::array<int, 2>>& mapping = agg_topologies[level]->map;
    if (level == 0) {
      BlockedCSCSymmetricMatrix3<T> child{sbd_topology.L.n, sbd_workspace.L_diagonal, sbd_workspace.L_lower, &sbd_topology.L.offsets[0],&sbd_topology.L.row[0]};
      // assume data of parent already initialized to all 0's
      for (int j=0; j<child.n; j++) {
        const int j_p = mapping[j][child_index]*3;

        if (j_p>=0) {
          // copy the diagonal terms over
          auto Mat_C = child.D(j);
          parent(j_p,   j_p)   += Mat_C[0];
          parent(j_p+1, j_p)   += Mat_C[1];
          parent(j_p+2, j_p)   += Mat_C[2];
          parent(j_p+1, j_p+1) += Mat_C[3];
          parent(j_p+2, j_p+1) += Mat_C[4];
          parent(j_p+2, j_p+2) += Mat_C[5];
          for (int ii = child.offsets[j]; ii<child.offsets[j+1]; ii++) {
            const int i = child.row[ii];

            const int i_p = mapping[i][child_index]*3;
            if (i_p >= 0)
              {
                //assert(i_p != j_p && i_p>=0)
                auto Mat_C = child.L(ii);
                if (i_p>j_p) {
                  for (int v=0; v<9; v++)
                    parent(i_p+v%3,j_p+v/3) += Mat_C[v];
                } else {
                  // add the transpose
                  for (int v=0; v<9; v++)
                    parent(j_p+v/3,i_p+v%3) += Mat_C[v];
                }}}}}
    } else {
	  Aggregate_Data<T>& child=agg_workspace[level-1];
      for (int j=0; j<child.n/d; j++) {
        const int j_p = mapping[j][child_index]*d;
        if (j_p>=0) {
          for (int w=0; w<d; w++)
            for (int v=w; v<d; v++)
              parent(j_p+v, j_p+w) += child(j*d+v, j*d+w);
          for (int i = j+1; i<child.n/d; i++) {
            const int i_p = mapping[i][child_index]*d;
            if (i_p>=0) {
              if (i_p>j_p) {
                for (int w=0; w<d; w++)
                  for (int v=0; v<d; v++)
                    parent(i_p+v, j_p+w) += child(i*d+v, j*d+w);
              } else {
                // add the transpose
                for (int w=0; w<d; w++)
                  for (int v=0; v<d; v++)
                    parent(j_p+w, i_p+v) += child(i*d+v, j*d+w);
              }}}}}
    }
  }

  template<class T>
  void Hierarchical_Cholesky<T>::Accumulate_to_Next_Level_Naive(const int level,  const int child_index) {
    //PhysBAM::LOG::SCOPE scope("Accumulate_to_Next_Level_Naive");
    Aggregate_Data<T>& parent=agg_workspace[level];
    const std::vector<std::array<int, 2>>& mapping = agg_topologies[level]->map;
    if (level == 0) {
      BlockedCSCSymmetricMatrix3<T> child{sbd_topology.L.n, sbd_workspace.L_diagonal, sbd_workspace.L_lower, &sbd_topology.L.offsets[0],&sbd_topology.L.row[0]};
      // assume data of parent already initialized to all 0's
      for (int j=0; j<child.n; j++) {
        const int j_p = mapping[j][child_index]*3;

        if (j_p>=0) {
          // copy the diagonal terms over
          auto Mat_C = child.D(j);
          parent(j_p,   j_p)   += Mat_C[0];
          parent(j_p+1, j_p)   += Mat_C[1];
          parent(j_p+2, j_p)   += Mat_C[2];
          parent(j_p+1, j_p+1) += Mat_C[3];
          parent(j_p+2, j_p+1) += Mat_C[4];
          parent(j_p+2, j_p+2) += Mat_C[5];
          for (int ii = child.offsets[j]; ii<child.offsets[j+1]; ii++) {
            const int i = child.row[ii];

            const int i_p = mapping[i][child_index]*3;
            if (i_p >= 0)
              {
                //assert(i_p != j_p && i_p>=0)
                auto Mat_C = child.L(ii);
                if (i_p>j_p) {
                  for (int v=0; v<9; v++)
                    parent(i_p+v%3,j_p+v/3) += Mat_C[v];
                } else {
                  // add the transpose
                  for (int v=0; v<9; v++)
                    parent(j_p+v/3,i_p+v%3) += Mat_C[v];
                }}}}}
    } else {
	  Aggregate_Data<T>& child=agg_workspace[level-1];
      for (int j=0; j<child.n/d; j++) {
        const int j_p = mapping[j][child_index]*d;
        if (j_p>=0) {
          for (int w=0; w<d; w++)
            for (int v=w; v<d; v++)
              parent(j_p+v, j_p+w) += child(j*d+v, j*d+w);
          for (int i = j+1; i<child.n/d; i++) {
            const int i_p = mapping[i][child_index]*d;
            if (i_p>=0) {
              if (i_p>j_p) {
                for (int w=0; w<d; w++)
                  for (int v=0; v<d; v++)
                    parent(i_p+v, j_p+w) += child(i*d+v, j*d+w);
              } else {
                // add the transpose
                for (int w=0; w<d; w++)
                  for (int v=0; v<d; v++)
                    parent(j_p+w, i_p+v) += child(i*d+v, j*d+w);
              }}}}}
    }
  }

  //#####################################################################
  // Function Hierarchical_Cholesky::Allocate_Space_For_Schur_Complements
  //#####################################################################
  template<class T>
  void Hierarchical_Cholesky<T>::Allocate_Space_For_Schur_Complements()
  {
      long long int num=(1<<depth)*(6*sbd_csc_topology.n+9*sbd_csc_topology.offsets[sbd_csc_topology.n]);
    //std::cout<<"sbd schur comp = "<<(1<<dpeth)*(6*sbd_csc_topology.n+9*sbd_csc_topology.offsets[sbd_csc_topology.n])<<" float"<<std::endl;
    for (int l=1; l<=depth; l++)
        num += (size_t) (1<<(depth-l))*(agg_topologies[l-1]->Ninterface*3)*(agg_topologies[l-1]->Ninterface*3+1)/2;
    std::cout<<"schur comp = "<<num<<" float"<<std::endl;

    //PhysBAM::LOG::SCOPE scope("Allocate_Space_For_Schur_Complements");
    sbd_schur.resize(1<<depth);
    for (int r=0; r<1<<depth; r++) {
        sbd_schur[r].L_diagonal = new T[(size_t) 6*sbd_csc_topology.n];
        sbd_schur[r].L_lower = new T[(size_t)9*sbd_csc_topology.offsets[sbd_csc_topology.n]];
    }
    #if 1
    agg_schur.resize(depth);
    for (int l=1; l<=depth; l++) {
      agg_schur[l-1].resize(1<<(depth-l));
      for (int r=0; r<1<<(depth-l); r++) {
          agg_schur[l-1][r] = new T[(size_t)(agg_topologies[l-1]->Ninterface*3)*(agg_topologies[l-1]->Ninterface*3+1)/2];
      }
	}
    #endif
  }


  //#####################################################################
  // Function Hierarchical_Cholesky::Compute_Schur_Complements
  //#####################################################################
  template<class T>
  void Hierarchical_Cholesky<T>::Compute_Schur_Complements(int level, int rank) {
    if (level == 0) {
      BlockedCSCSymmetricMatrix3<T> child{sbd_topology.L.n, static_cast<void*>(sbd_workspace.L_diagonal), static_cast<void*>(sbd_workspace.L_lower), &sbd_topology.L.offsets[0], &sbd_topology.L.row[0]};
      Fill_Subdomain_CSC_Data(child, rank);
      {
        //PhysBAM::LOG::SCOPE scope("Subdomain_Cholesky");
        Sparse_Cholesky_Partial<T,T,int>(child, sbd_topology.Nsubdomain);
      }
      if (depth>0)
        Accumulate_to_Next_Level(level, rank%2);
      Extract_Schur(level, rank);
      // clear workspace
      std::fill(sbd_workspace.L_diagonal, sbd_workspace.L_diagonal+6*sbd_topology.L.n,T());
      std::fill(sbd_workspace.L_lower, sbd_workspace.L_lower+6*sbd_topology.L.offsets[sbd_topology.L.n],T());
    } else {
	  agg_workspace[level-1].A_ii = agg_schur[level-1][rank];
      std::fill(agg_workspace[level-1].A_ii, agg_workspace[level-1].A_ii+(agg_workspace[level-1].k)*(agg_workspace[level-1].k+1)/2,T());
      if (level<depth) {
        std::fill(agg_workspace[level-1].A_ir, agg_workspace[level-1].A_ir+(agg_workspace[level-1].k)*(agg_workspace[level-1].n-agg_workspace[level-1].k),T());
        std::fill(agg_workspace[level-1].A_rr, agg_workspace[level-1].A_rr+(agg_workspace[level-1].n-agg_workspace[level-1].k)*(agg_workspace[level-1].n-agg_workspace[level-1].k+1)/2,T());
      }

      Compute_Schur_Complements(level-1, rank*2);
      Compute_Schur_Complements(level-1, rank*2+1);
      // cholesky
      Cholesky_on_Aggregate(level);
      if (level<depth)
        Accumulate_to_Next_Level(level, rank%2);
      //  Extract_Schur(level, rank);
      //std::fill(agg_workspace[level-1].data, agg_workspace[level-1].data+(agg_topologies[level-1]->Nvariables*3)*(agg_topologies[level-1]->Nvariables*3+1)/2,T());

    }
  }

  template<class T>
  void Hierarchical_Cholesky<T>::Compute_Schur_Complements() {
      //std::cout<<"schur"<<std::endl;
      sbd_workspace.L_diagonal = new T[(size_t)6*sbd_topology.L.n]();
      sbd_workspace.L_lower = new T[(size_t)9*sbd_topology.L.offsets[sbd_topology.L.n]]();

      agg_workspace.resize(depth);
      for (int l=1; l<=depth; l++) {
          agg_workspace[l-1].n = agg_topologies[l-1]->Nvariables*3;
          agg_workspace[l-1].k = agg_topologies[l-1]->Ninterface*3;
          if (l<depth) {
              agg_workspace[l-1].A_ir= new T[(size_t)(agg_workspace[l-1].k)*(agg_workspace[l-1].n-agg_workspace[l-1].k)]();
              agg_workspace[l-1].A_rr= new T[(size_t)(agg_workspace[l-1].n-agg_workspace[l-1].k)*(agg_workspace[l-1].n-agg_workspace[l-1].k+1)/2]();
          }
      }
      Compute_Schur_Complements(depth, 0);
      long long int num = 6*sbd_topology.L.n+9*sbd_topology.L.offsets[sbd_topology.L.n];
      for (int l=1; l<depth; l++) {
          num+=(agg_workspace[l-1].k)*(agg_workspace[l-1].n-agg_workspace[l-1].k)+(agg_workspace[l-1].n-agg_workspace[l-1].k)*(agg_workspace[l-1].n-agg_workspace[l-1].k+1)/2;
      }
      std::cout<<"workspace = "<<num<<" float"<<std::endl;

      delete[] sbd_workspace.L_diagonal;
      delete[] sbd_workspace.L_lower;
      for (int l=1; l<depth; l++) {
          delete[] agg_workspace[l-1].A_ir;
          delete[] agg_workspace[l-1].A_rr;
      }
  }

  //#####################################################################
  // Function Hierarchical_Cholesky::Fill_Subdomain_CSC_Data
  //#####################################################################
  template<class T>
  void Hierarchical_Cholesky<T>::Fill_Subdomain_CSC_Data(BlockedCSCSymmetricMatrix3<T>& CSC, const int rank)
  {
    //PhysBAM::LOG::SCOPE scope("Fill_Subdomain_CSC_Data");
    // assert (CSC.n == coordinates.size())
    const Coord_Type& sbd_size = sbd_topology.subdomain_size;
    const Subdomain_Topology::Coord_Array_Type& coordinates=sbd_topology.coordinates;
    const Blocked_Array_3D<T> data_array{domain_size[0], domain_size[1], domain_size[2], data};

    Coord_Type min_corner;
    Coord_Type Sstep;
    Aggregate_Topology::Get_sbd_Mapping(min_corner, Sstep, sbd_size,  rank);

    using Matrix3_type = T (&)[9];
    using Matrix3_ptr = T (*)[9];
    using Matrix3_const_type =  const T (&)[9];
    using Lower_Triangular_Matrix3_type = T (&)[6];
    using Lower_Triangular_Matrix3_ptr = T (*)[6];
    using Lower_Triangular_Matrix3_const_type =  const T (&)[6];

    for (int j=0; j<CSC.n; j++) {
      // filling diagonal terms
      const Coord_Type& index_j = coordinates[j];
      Coord_Type coord_j{index_j[0]*Sstep[0]+min_corner[0], index_j[1]*Sstep[1]+min_corner[1], index_j[2]*Sstep[2]+min_corner[2]};

      if (coord_j[0]>=1 && coord_j[0]<=domain_size[0] && coord_j[1]>=1 && coord_j[1]<=domain_size[1] && coord_j[2]>=1 && coord_j[2]<=domain_size[2]) {
        const auto k_data = data_array(coord_j[0]-1,coord_j[1]-1,coord_j[2]-1,1,1,1);
        // why does this work? Reference to a const temporary object?
        Lower_Triangular_Matrix3_const_type local {
          k_data[0], k_data[3], k_data[6],
                     k_data[4], k_data[7],
                                k_data[8]
            };

        if (j<sbd_topology.Nsubdomain) {
          for (int v=0; v<6; v++ )
            CSC.D(j)[v] = local[v];
        } else {
          int nExtr = 0;
          for (int v=0; v<d; v++)
            if (index_j[v] == 0 || index_j[v] == sbd_size[v]+1)
              nExtr++;
          for (int v=0; v<6; v++ )
            CSC.D(j)[v] = local[v]/(1<<nExtr);
        }

        for (int ii=CSC.offsets[j]; ii<CSC.offsets[j+1]; ii++) {
          // filling off-diagonal terms
          const int i = CSC.row[ii];
          const Coord_Type& index_i = coordinates[i];
          const Coord_Type coord_i{index_i[0]*Sstep[0]+min_corner[0], index_i[1]*Sstep[1]+min_corner[1], index_i[2]*Sstep[2]+min_corner[2]};
          const Coord_Type delta{coord_j[0]-coord_i[0], coord_j[1]-coord_i[1], coord_j[2]-coord_i[2]};
          const Coord_Type delta_abs{abs(delta[0]), abs(delta[1]), abs(delta[2])};

          if (delta_abs[0]<=1 && delta_abs[1]<=1 && delta_abs[2]<=1 && coord_i[0]>=1 && coord_i[0]<=domain_size[0] && coord_i[1]>=1 && coord_i[1]<=domain_size[1] && coord_i[2]>=1 && coord_i[2]<=domain_size[2]) {
            const auto k_data = data_array(coord_i[0]-1,coord_i[1]-1,coord_i[2]-1,delta[0]+1,delta[1]+1,delta[2]+1);
            Matrix3_const_type local {
              k_data[0], k_data[3], k_data[6],
                k_data[1], k_data[4], k_data[7],
                k_data[2], k_data[5], k_data[8]
                };

            if (i<sbd_topology.Nsubdomain) {
              for (int v=0; v<9; v++)
                CSC.L(ii)[v] = local[v];
            } else {
              int nExtr = 0;
              for (int v=0; v<d; v++)
                if ((index_i[v] == 0 || index_i[v] == sbd_size[v]+1) && index_i[v] == index_j[v])
                  nExtr++;
              for (int v=0; v<9; v++ )
                CSC.L(ii)[v] = local[v]/(1<<nExtr);
            }
          } else {
            for (int v=0; v<9; v++ )
              CSC.L(ii)[v] = 0;
          }}
      } else {
        for (int v=0; v<6; v++ )
          CSC.D(j)[v] = T();
        for (int ii=CSC.offsets[j]; ii<CSC.offsets[j+1]; ii++) {
          for (int v=0; v<9; v++ )
            CSC.L(ii)[v] = T();
        }
      }}
  }

  //#####################################################################
  // Function Hierarchical_Cholesky::Cholesky_on_Aggregate
  //#####################################################################
  template<>
  void Hierarchical_Cholesky<float>::Cholesky_on_Aggregate(const int level) {
    //PhysBAM::LOG::SCOPE scope("Cholesky_on_Aggregate");
    lapack_int status = LAPACKE_spptrf (LAPACK_COL_MAJOR , 'L' , agg_workspace[level-1].k,  agg_workspace[level-1].A_ii);
    //std::cout<<std::endl<<status<<std::endl;
    if (level<depth) {
      status = LAPACKE_stptrs (LAPACK_COL_MAJOR , 'L' , 'N' , 'N' , agg_workspace[level-1].k , agg_workspace[level-1].n-agg_workspace[level-1].k, agg_workspace[level-1].A_ii, agg_workspace[level-1].A_ir, agg_workspace[level-1].k);
      //std::cout<<status<<std::endl;
      status = LAPACKE_ssfrk (LAPACK_COL_MAJOR , 'N', 'L', 'T', agg_workspace[level-1].n-agg_workspace[level-1].k, agg_workspace[level-1].k , -1. , agg_workspace[level-1].A_ir , agg_workspace[level-1].k, 1., agg_workspace[level-1].A_rr);
      //std::cout<<status<<std::endl;
    }
  }
  //#####################################################################
  template class Hierarchical_Cholesky<float>;
  //    template class Hierarchical_Cholesky<double>;
}
