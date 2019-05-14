//#####################################################################
// Copyright (c) 2018, Eftychios Sifakis, Qisi Wang
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include "Subdomain_Topology.h"
#include "Aggregate_Topology.h"
#include <iostream>
#include <iomanip>

namespace {
    void Aggregate_Morton(Macroblock_Utilities::Aggregate_Topology& agg, const int child_Nvariables, const int child_Ninterface, const std::array<int,3>& child_size, const typename Macroblock_Utilities::Aggregate_Topology::Hashtable_Type& indices)
  {
    using namespace Macroblock_Utilities;
    constexpr int d = 3;
    using Coord_Type = std::array<int,d>;

    // setup
    const int aggregation_axis=2-(agg.level-1)%3;
    agg.map.resize(child_Nvariables);

    int count_i=0;
    int count_v=agg.Ninterface;
    // First batch of subdomain indices correspond to interior nodes; those are not mapped to the aggregate
    for(int sbd_n=0;sbd_n<child_Ninterface;sbd_n++)
      agg.map[sbd_n]=std::array<int,2>{-1,-1};


    // create map (iterate through geometric domain)
    // iterating through top child (0 indexed)
    {
        // shared face first
      Coord_Type range[2] = {{0,0,0}, {child_size[0]-1, child_size[1]-1, child_size[2]-1}};
      range[0][aggregation_axis] = range[1][aggregation_axis];
      for (int i=range[0][0]; i<=range[1][0]; i++)
      for (int j=range[0][1]; j<=range[1][1]; j++)
      for (int k=range[0][2]; k<=range[1][2]; k++) {
        Coord_Type coord_upper = {i,j,k};
        Coord_Type coord_lower = coord_upper;
        coord_lower[aggregation_axis] = child_size[aggregation_axis]-1-coord_upper[aggregation_axis];
        // test if interior of the face
        bool interior = true;
        for (int v=(aggregation_axis+1)%d; v!=aggregation_axis; v=(v+1)%d)
          if (coord_upper[v] == range[0][v] || coord_upper[v] == range[1][v])
            {interior = false; break;}
        if (!interior) {
            agg.map[indices.at(coord_upper)][0] = agg.map[indices.at(coord_lower)][1] = count_v;
            agg.index[coord_upper] = count_v++;
        }
        else {
          agg.map[indices.at(coord_upper)][0] = agg.map[indices.at(coord_lower)][1] = count_i;
          agg.index[coord_upper] = count_i++;
        }
      }
    }


    {
      // 4 faces
      Coord_Type range[2] = {{0,0,0}, {child_size[0]-1, child_size[1]-1, child_size[2]-1}};
      range[1][aggregation_axis] -= 1;
      range[0][aggregation_axis] = 1;
      for (int i=range[0][0]; i<=range[1][0]; i++)
      for (int j=range[0][1]; j<=range[1][1]; j++)
      for (int k=range[0][2]; k<=range[1][2]; k++) {
        Coord_Type coord_upper = {i,j,k};
        Coord_Type coord_lower = coord_upper;
        coord_lower[aggregation_axis] = child_size[aggregation_axis]-1-coord_upper[aggregation_axis];
        // test if interior of the face
        bool interior = true;
        for (int v=(aggregation_axis+1)%d; v!=aggregation_axis; v=(v+1)%d)
          if (coord_upper[v] == range[0][v] || coord_upper[v] == range[1][v])
            {interior = false; break;}
        if (!interior) {
          agg.map[indices.at(coord_upper)][0] = count_v;
          agg.index[coord_upper] = count_v++;
          agg.map[indices.at(coord_lower)][1] = count_v;
          coord_lower[aggregation_axis]+=child_size[aggregation_axis]-1;
          agg.index[coord_lower] = count_v++;
        }
      }
    }


    {
      // last face
      Coord_Type range[2] = {{0,0,0}, {child_size[0]-1, child_size[1]-1, child_size[2]-1}};
      range[1][aggregation_axis] = range[0][aggregation_axis];
      for (int i=range[0][0]; i<=range[1][0]; i++)
      for (int j=range[0][1]; j<=range[1][1]; j++)
      for (int k=range[0][2]; k<=range[1][2]; k++) {
        Coord_Type coord_upper = {i,j,k};
        Coord_Type coord_lower = coord_upper;
        coord_lower[aggregation_axis] = child_size[aggregation_axis]-1-coord_upper[aggregation_axis];

        agg.map[indices.at(coord_upper)][0] = count_v;
        agg.index[coord_upper] = count_v++;
        agg.map[indices.at(coord_lower)][1] = count_v;
        coord_lower[aggregation_axis]+=child_size[aggregation_axis]-1;
        agg.index[coord_lower] = count_v++;
      }
    }

    // std::cout<<agg.Nvariables<<" "<<count_v<<std::endl;
    // std::cout<<agg.Ninterface<<" "<<count_i<<std::endl;
  }

    void Aggregate(Macroblock_Utilities::Aggregate_Topology& agg, const std::array<int,3>& bdr_idx, const int aggregation_axis, const int child_Nvariables, const int child_Ninterface, const std::vector<std::array<int,3>> & child_coordinates, int Nshare) {
        // std::cout<<"Nshare = "<<Nshare<<std::endl;
      /*
       * Interface
       * rank 0 children nodes
       * rank 1 children nodes
       * other shared nodes of the 2 children
       */

    using namespace Macroblock_Utilities;
    constexpr int d = 3;
    using Coord_Type = std::array<int,d>;
    agg.map.resize(child_Nvariables);
    int count_i=0;
    int count_l=agg.Ninterface;
    int count_u = (agg.Nvariables-((bdr_idx[aggregation_axis]>2?2:bdr_idx[aggregation_axis])+1)*(agg.Ninterface+Nshare))/2+(bdr_idx[aggregation_axis]>1?(agg.Ninterface+Nshare):0)+agg.Ninterface;
    int count_s=agg.Nvariables-Nshare;

    // std::cout<<"agg.Ninterface = "<<agg.Ninterface<<std::endl;
    // std::cout<<"agg.Nvariables = "<<agg.Nvariables<<std::endl;
    // std::cout<<"Nshare = "<<Nshare<<std::endl;
    // std::cout<<"bdr_index[aggregation_axis] = "<<bdr_idx[aggregation_axis]<<std::endl;
    // std::cout<<"non_share_offset = "<<non_share_offset<<std::endl;


    // First batch of subdomain indices correspond to interior nodes; those are not mapped to the aggregate
    for(int sbd_n=0;sbd_n<child_Ninterface;sbd_n++)
      agg.map[sbd_n]=std::array<int,2>{-1,-1};

    for(int sbd_n=child_Ninterface;sbd_n<child_Nvariables;sbd_n++){
      Coord_Type coord_lower(child_coordinates[sbd_n]);
      Coord_Type coord_upper(child_coordinates[sbd_n]);
      coord_upper[aggregation_axis]=agg.aggregate_size[aggregation_axis]-coord_upper[aggregation_axis]-1;

      if (coord_lower[aggregation_axis]==coord_upper[aggregation_axis]) {
        bool is_interior = true;
        bool on_domain_bdr = false;
        for (int v=(aggregation_axis+1)%d; v!=aggregation_axis&&!on_domain_bdr; v=(v+1)%d)
          {
            if (coord_lower[v]==0) {
              is_interior=false;
              if (bdr_idx[v] < 2)
                on_domain_bdr = true;
            }
            else if (coord_lower[v]==agg.aggregate_size[v]-1) {
              is_interior=false;
              if (bdr_idx[v] == 0)
                on_domain_bdr = true;
            }
          }
        if (on_domain_bdr) {
            agg.map[sbd_n] = std::array<int,2>{-2, -2};
            if (agg.level > 1)
              std::cerr<<"at agg.level "<<agg.level<<" mapped out of domain bdr, yet was indexed in child"<<std::endl;
          }
        else {
            if (is_interior) {
              agg.map[sbd_n] = std::array<int,2>{count_i,count_i};
              agg.index.insert({coord_lower,count_i++});
            }
            else {
              agg.map[sbd_n] = std::array<int,2>{count_s,count_s};
              agg.index.insert({coord_lower,count_s++});
            }
          }
      } else {
        // determine if on domain bdr, 0: none on interior of domain; 1: upper side interior of domain; 2: both side interiror of domain
        int n_domain_interior = 2;
        for (int v=(aggregation_axis+1)%d; v!=aggregation_axis&&n_domain_interior; v=(v+1)%d) {
            if (bdr_idx[v] == 1) {
                if (coord_lower[v] == 0)
                  n_domain_interior = 0;
              }
            else if (bdr_idx[v] == 0) {
                if (coord_lower[v] == 0 || coord_lower[v] == agg.aggregate_size[v]-1)
                  n_domain_interior = 0;
              }
        }
        if (n_domain_interior) {
            if (coord_lower[aggregation_axis] == 0)
              n_domain_interior = bdr_idx[aggregation_axis]<2?bdr_idx[aggregation_axis]:2;
        }

        switch (n_domain_interior)
          {
          case 0:
            agg.map[sbd_n] = std::array<int,2>{-2, -2}; // this case should never exist when level > 1, if everything is corret
            if (agg.level > 1)
              std::cerr<<"at agg.level "<<agg.level<<" mapped out of domain bdr, yet was indexed in child"<<std::endl;
            break;
          case 1:
            agg.map[sbd_n][0] = -2;
            agg.map[sbd_n][1] = count_u;
            agg.index.insert({coord_upper,count_u++});
            break;
          case 2:
            agg.map[sbd_n][0] = count_l;
            agg.index.insert({coord_lower,count_l++});
            agg.map[sbd_n][1] = count_u;
            agg.index.insert({coord_upper,count_u++});
            break;
          default:
            std::cerr<<"impossible n_domain_interior"<<std::endl;
            exit(1);
          }
      }
    }
    if(count_i!=agg.Ninterface) {
        std::cerr<<"Mismatch in number of interface variables at aggregate agg.level "<<agg.level<<": count_i = "<<count_i<<", agg.Ninterface = "<< agg.Ninterface<<std::endl;
        //        exit(1);
    }
    if(count_l!= (agg.Nvariables-((bdr_idx[aggregation_axis]>2?2:bdr_idx[aggregation_axis])+1)*(agg.Ninterface+Nshare))/2+(bdr_idx[aggregation_axis]>1?(agg.Ninterface+Nshare):0)+agg.Ninterface) {
        std::cerr<<"Mismatch in number of interface variables at aggregate agg.level "<<agg.level<<": count_l = "<<count_l<<", agg.Ninterface = "<< agg.Ninterface<<std::endl;
        //        exit(1);
    }
    if (count_u!=agg.Nvariables-Nshare) {
        std::cerr<<"Mismatch in number of interface variables at aggregate agg.level "<<agg.level<<": count_u = "<<count_u<<", agg.Ninterface = "<< agg.Ninterface<<std::endl;
      //        exit(1);
    }
    if(count_s!=agg.Nvariables) {
      std::cerr<<"Mismatch in number of total variables at aggregate agg.level "<<agg.level<<": count_s = "<<count_s<<", agg.Nvariables = "<< agg.Nvariables<<std::endl;
      // exit(1);
    }
  }
}


namespace Macroblock_Utilities {

//#####################################################################
// Function Subdomain::Initialize_Coordinates
//#####################################################################
void Aggregate_Topology::
Initialize_Coordinates()
{
    coordinates.resize(index.size());
    for (auto it : index)
        coordinates[it.second] = it.first;
}

//#####################################################################
// Function Aggregate_Topology::Aggregate (from Subdomain, with depth)
//#####################################################################
Aggregate_Topology* Aggregate_Topology::
Aggregate_Morton(const Subdomain_Topology& child, const int depth)
{
    int level=1;
    int aggregation_axis=2-(level-1)%3;
    // const std::array<int,d> bdr_idx = {
    //   (depth/3-level/3),
    //   (depth/3+(depth%3>1?1:0)-level/3-(level%3>1?1:0)),
    //   (depth/3+(depth%3>0?1:0)-level/3-(level%3>0?1:0)),
    // }; // there should exists simpler logic
    // std::cout<<"bdr_idx = ["<<bdr_idx[2]<<", "<<bdr_idx[1]<<", "<<bdr_idx[2]<<"]"<<std::endl;
    Coord_Type aggregate_size;
    for(int v=0;v<d;v++)
        if(v==aggregation_axis)
            aggregate_size[v]=2*child.subdomain_size[v]+3;
        else
            aggregate_size[v]=child.subdomain_size[v]+2;

    int Nvariables=1;
    {
      // for (int v = 0; v < d; ++v) {
      //   if (bdr_idx[v] == 0)
      //     Nvariables *= aggregate_size[v] - 2;
      //   else if (bdr_idx[v] == 1)
      //     Nvariables *= aggregate_size[v] - 1;
      //   else
      //     Nvariables *= aggregate_size[v];
      // }
      for(auto x:aggregate_size) Nvariables*=x;
      int Ninterior=1;
      for(auto x:aggregate_size) Ninterior*=(x-2);
      Nvariables-=Ninterior;
    }

      int Ninterface=1;
      for(int v=0;v<d;v++)
        if(v!=aggregation_axis)
          Ninterface*=(aggregate_size[v]-2);
      Nvariables+=Ninterface;
      Aggregate_Topology& agg = *(new Aggregate_Topology{level, aggregate_size, Nvariables, Ninterface});


      ::Aggregate_Morton(agg, child.Nvariables, child.Nsubdomain, {child.subdomain_size[0]+2,child.subdomain_size[1]+2,child.subdomain_size[2]+2},child.index);
      return &agg;
    }

Aggregate_Topology* Aggregate_Topology::
Aggregate(const Subdomain_Topology& child, const int depth)
{
    int level=1;
    const int aggregation_axis=Aggregate_Axis(level);
    const std::array<int,d> total_agg = Aggregation_In_Each_Dimension(depth);
    const std::array<int,d> cur_agg = Aggregation_In_Each_Dimension(level);
    const std::array<int,d> bdr_idx = {total_agg[0]-cur_agg[0], total_agg[1]-cur_agg[1], total_agg[2]-cur_agg[2]};
    // std::cout<<"bdr_idx = ["<<bdr_idx[2]<<", "<<bdr_idx[1]<<", "<<bdr_idx[2]<<"]"<<std::endl;
    Coord_Type aggregate_size;
    for(int v=0;v<d;v++)
        if(v==aggregation_axis)
            aggregate_size[v]=2*child.subdomain_size[v]+3;
        else
            aggregate_size[v]=child.subdomain_size[v]+2;
    int Nvariables=1;
    {
      for (int v = 0; v < d; ++v) {
        if (bdr_idx[v] == 0)
          Nvariables *= aggregate_size[v] - 2;
        else if (bdr_idx[v] == 1)
          Nvariables *= aggregate_size[v] - 1;
        else
          Nvariables *= aggregate_size[v];
      }
      int Ninterior=1;
      for(auto x:aggregate_size) Ninterior*=(x-2);
      Nvariables-=Ninterior;
    }
    int Nshare=1;
    {
        for (int v = (aggregation_axis+1)%d; v != aggregation_axis; v=(v+1)%d) {
            if (bdr_idx[v] == 0)
                Nshare *= aggregate_size[v] - 2;
            else if (bdr_idx[v] == 1)
                Nshare *= aggregate_size[v] - 1;
            else
                Nshare *= aggregate_size[v];
        }
    }


      int Ninterface=1;
      for(int v=0;v<d;v++)
        if(v!=aggregation_axis)
          Ninterface*=(aggregate_size[v]-2);
      Nvariables+=Ninterface;
      Nshare -= Ninterface;
      Aggregate_Topology& agg = *(new Aggregate_Topology{level, aggregate_size, Nvariables, Ninterface});

      ::Aggregate(agg, bdr_idx, aggregation_axis, child.Nvariables, child.Nsubdomain, child.coordinates, Nshare);

      return &agg;
    }

//#####################################################################
// Function Aggregate_Topology::Aggregate (from Aggregate, with depth)
//#####################################################################
  Aggregate_Topology* Aggregate_Topology::
Aggregate(const Aggregate_Topology& child, const int depth)
{
    const int level=child.level+1;
    const int aggregation_axis=Aggregate_Axis(level);
    Coord_Type aggregate_size;
    for(int v=0;v<d;v++)
        if(v==aggregation_axis)
            aggregate_size[v]=2*child.aggregate_size[v]-1;
        else
            aggregate_size[v]=child.aggregate_size[v];
    int Ninterface=1;
    for(int v=0;v<d;v++)
        if(v!=aggregation_axis)
            Ninterface*=(aggregate_size[v]-2);
    const std::array<int,d> total_agg = Aggregation_In_Each_Dimension(depth);
    const std::array<int,d> cur_agg = Aggregation_In_Each_Dimension(level);
    const std::array<int,d> bdr_idx = {total_agg[0]-cur_agg[0], total_agg[1]-cur_agg[1], total_agg[2]-cur_agg[2]};
    //std::cout<<"bdr_idx = ["<<bdr_idx[2]<<", "<<bdr_idx[1]<<", "<<bdr_idx[2]<<"]"<<std::endl;

    int Nvariables=1;
    {
      for (int v = 0; v < d; ++v) {
        if (bdr_idx[v] == 0)
          Nvariables *= aggregate_size[v] - 2;
        else if (bdr_idx[v] == 1)
          Nvariables *= aggregate_size[v] - 1;
        else
          Nvariables *= aggregate_size[v];
      }
        int Ninterior=1;
        for(auto x:aggregate_size) Ninterior*=(x-2);
        Nvariables-=Ninterior;
    }

    Nvariables+=Ninterface;

    int Nshare=1;
    {
        for (int v = (aggregation_axis+1)%d; v != aggregation_axis; v=(v+1)%d) {
            if (bdr_idx[v] == 0)
                Nshare *= aggregate_size[v] - 2;
            else if (bdr_idx[v] == 1)
                Nshare *= aggregate_size[v] - 1;
            else
                Nshare *= aggregate_size[v];
        }
    }
    Nshare-=Ninterface;

    Aggregate_Topology& agg = *(new Aggregate_Topology{level, aggregate_size, Nvariables, Ninterface});

    ::Aggregate(agg, bdr_idx, aggregation_axis, child.Nvariables, child.Ninterface, child.coordinates,Nshare);
    return &agg;
}


int Aggregate_Topology::
Aggregate_Axis(const int level)
{
  // the aggregation direction to form the level given
  return level==1?0:level==2?1:level==3?2:2-(level-1)%3; // xyzzyxzyx
}


std::array<int,3> Aggregate_Topology::
Aggregation_In_Each_Dimension(const int level)
{
  // number of aggregation in each direction all the way to the given level
  // Aggregation_In_Each_Dimension(level)[v] = Aggregation_In_Each_Dimension(level-1)[v]+1 (if v==Aggregate_Axis(level))
  //                                           Aggregation_In_Each_Dimension(level-1)[v]   (if v!=Aggregate_Axis(level))
  switch (level) {
  case 1:
    return {1,0,0};
  case 2:
    return {1,1,0};
  default:
    return {level/3, level/3+(level%3>1?1:0), level/3+(level%3>0?1:0)};
  }
}

Aggregate_Topology* Aggregate_Topology::
Aggregate_Morton(const Aggregate_Topology& child, const int depth)
{
    const int level=child.level+1;
    const int aggregation_axis=2-(level-1)%3;
    Coord_Type aggregate_size;
    for(int v=0;v<d;v++)
        if(v==aggregation_axis)
            aggregate_size[v]=2*child.aggregate_size[v]-1;
        else
            aggregate_size[v]=child.aggregate_size[v];
    int Ninterface=1;
    for(int v=0;v<d;v++)
        if(v!=aggregation_axis)
            Ninterface*=(aggregate_size[v]-2);

    // const std::array<int,d> bdr_idx = {
    //   (depth/3-level/3),
    //   (depth/3+(depth%3>1?1:0)-level/3-(level%3>1?1:0)),
    //   (depth/3+(depth%3>0?1:0)-level/3-(level%3>0?1:0)),
    // }; // there should exists simpler logic
    //std::cout<<"bdr_idx = ["<<bdr_idx[2]<<", "<<bdr_idx[1]<<", "<<bdr_idx[2]<<"]"<<std::endl;

    int Nvariables=1;
    {
      //for (int v = 0; v < d; ++v) {
        // if (bdr_idx[v] == 0)
        //   Nvariables *= aggregate_size[v] - 2;
        // else if (bdr_idx[v] == 1)
        //   Nvariables *= aggregate_size[v] - 1;
        // else
        //Nvariables *= aggregate_size[v];
          //}
      for(auto x:aggregate_size) Nvariables*=x;
      int Ninterior=1;
        for(auto x:aggregate_size) Ninterior*=(x-2);
        Nvariables-=Ninterior;
    }
    Nvariables+=Ninterface;

    Aggregate_Topology& agg = *(new Aggregate_Topology{level, aggregate_size, Nvariables, Ninterface});

    ::Aggregate_Morton(agg, child.Nvariables, child.Ninterface, child.aggregate_size,child.index);
    //::Aggregate(agg, bdr_idx, child.Nvariables, child.Ninterface, child.coordinates);
    return &agg;
}

//#####################################################################
// Function Get_sbd_Mapping
//#####################################################################
  void Aggregate_Topology::
  Get_sbd_Mapping(std::array<int, 3>& min_corner, std::array<int, 3>& Sstep, const std::array<int, 3>& sbd_size, const int rank) {
        constexpr int d = 3;
        using Coord_Type = std::array<int,d>;
        min_corner = Coord_Type{0};
        Sstep = Coord_Type{1,1,1};

        Coord_Type b_size{sbd_size[0]+1,sbd_size[1]+1,sbd_size[2]+1};

        for (int r=rank, s=0; r; r>>=1, s++) {
          const int axis=Aggregate_Axis(s+1);
          b_size[axis]<<=1;

            if (r&1) {
                min_corner[axis] = b_size[axis] - min_corner[axis];
                Sstep[axis] = -Sstep[axis];
            }
        }
    }

    void Aggregate_Topology::
  Get_sbd_Mapping_Morton(std::array<int, 3>& min_corner, std::array<int, 3>& Sstep, const std::array<int, 3>& sbd_size, const int rank) {
        constexpr int d = 3;
        using Coord_Type = std::array<int,d>;
        min_corner = Coord_Type{0};
        Sstep = Coord_Type{1,1,1};

        for (int b=0, r=rank; r; b++, r>>=1) {
          min_corner[2-b%d]|=((r&1)<<(b/3)); // need to fix, the first accum axis is not nessissary z with different sbd size
        }
      for (int v=0; v<d; v++)
        min_corner[v] *= (sbd_size[v]+1);

    }

//#####################################################################
// Function Aggregate_Mapping
//#####################################################################
    std::array<std::array<int,3>,2>Aggregate_Topology::
    Aggregate_Mapping(int level, int rank, const std::array<int, 3>& sbd_size) {
    // mapping from hierarchical rank and its level to the range of a aggregate (interfaces and subdomains and boundaries e.g. 5x5x5 for a subdomain)
    // agg_axes - aggregate axis from subdomains
    constexpr int d = 3;
    std::array<int,d> box_size{sbd_size[0]+1,sbd_size[1]+1,sbd_size[2]+1};
    //std::cout<<"inner level = "<<level<<std::endl;
    // std::cout<<"box_size = ["<<box_size[0]<<","<<box_size[1]<<","<<box_size[2]<<"]"<<std::endl;
    // std::cout<<level<<std::endl
    for (int l=0; l<level; l++) {
        box_size[Aggregate_Axis(l+1)] <<= 1;
        //box_size[agg_axes[l]] <<= 1;
    }
    //std::cout<<"box_size = ["<<box_size[0]<<","<<box_size[1]<<","<<box_size[2]<<"]"<<std::endl;


    std::array<int,d> min_corner{};
    std::array<int,d> max_corner(box_size);
    for (int r=rank,l=level; r; r>>=1, l++) {
//        const int axis = agg_axes[l];
        const int axis = Aggregate_Axis(l+1);
        box_size[axis]<<=1;
        //std::cout<<"box_size = ["<<box_size[0]<<","<<box_size[1]<<","<<box_size[2]<<"]"<<std::endl;
        if (r&1) {
            min_corner[axis] = box_size[axis] - min_corner[axis];
            max_corner[axis] = box_size[axis] - max_corner[axis];
        }
    }
    //std::cout<<"min_corner = ["<<min_corner[0]<<","<<min_corner[1]<<","<<min_corner[2]<<"]"<<std::endl;
    //std::cout<<"max_corner = ["<<max_corner[0]<<","<<max_corner[1]<<","<<max_corner[2]<<"]"<<std::endl;

    return {min_corner, max_corner};
}

//#####################################################################
// Function Level_Range
//#####################################################################
    std::array<std::array<int,3>,2>Aggregate_Topology::
    Level_Range(int level, int rank, int depth) {
        // it is assumed that the sbd_size = {3,3,3}
    // mapping from hierarchical rank and its level to the range of a level (interface or sbd)
    // agg_axes - aggregate axis from subdomains
    constexpr int d = 3;
    std::array<std::array<int, d>,2> range{0};
    if (level<=depth) {
        //const int axis = agg_axes[level-1];
        const int axis = Aggregate_Axis(level);
        std::array<int,d> bit_delta{1,1,1};
        std::array<int,d> bit_cnt{0,0,0};

        for (int l=level, r=rank; l<depth; l++, r>>=1) {
            //const int v = agg_axes[l];
            const int v = Aggregate_Axis(l+1);
            range[0][v]|=r&1?bit_delta[v]:0;
            bit_delta[v]<<=1;
            bit_cnt[v]++;
        }

// should reverse the bit traverse order
        bit_delta = {1,1,1};
        for (int v=0; v<d; v++) {
            for (;bit_cnt[v];bit_cnt[v]--) {
                if (range[0][v]&bit_delta[v])
                    range[0][v]^=(bit_delta[v]-1);
                bit_delta[v]<<=1;
            }
        }

        if (level) {
            for (int v=0; v<d ;v++) {
                if (v!=axis) {
                    range[1][v] = range[0][v]+1;
                    range[1][v] <<= 2;
                    range[0][v] <<= 2;
                } else {
                    range[0][v] <<= 1;
                    range[0][v] |= 1;
                    range[0][v] <<= 2;
                    range[1][v] = range[0][v];
                }
            }

            for (int l=1; l<level; l++) {
                for (int s=0; s<2; s++)
                    range[s][Aggregate_Axis(l)]<<=1;
                    //range[s][agg_axes[l-1]]<<=1;
            }
            for (int v=(axis+1)%d; v!=axis; v=(v+1)%d) {
                range[0][v] += 1;
                range[1][v] -= 1;
            }
        } else {
            for (int v=0; v<d ;v++) {
                range[1][v] = range[0][v]+1;
                range[1][v] <<= 2;
                range[0][v] <<= 2;
            }

            for (int l=1; l<=level; l++) {
                for (int s=0; s<2; s++)
                    range[s][Aggregate_Axis(l)]<<=1;
            }
            for (int v=0; v<d; v++) {
                range[0][v] += 1;
                range[1][v] -= 1;
            }
        }

    }
    return range;
}

//#####################################################################
// Function Level_Mapping
//#####################################################################
    std::array<std::array<int,3>,2>Aggregate_Topology::
    Level_Mapping(int level, int rank, const std::array<int, 3>& sbd_size) {
        //std::cout<<"level mapping level = "<<level<<std::endl;
        auto range = Aggregate_Mapping(level,rank,sbd_size);
        //std::cout<<"range = ["<<range[0][0]<<","<<range[0][1]<<","<<range[0][2]<<"] -> ["<<range[1][0]<<","<<range[1][1]<<","<<range[1][2]<<"]"<<std::endl;

        for (int v=0; v<d; v++) {
            // remove bdr dof's
            if (range[0][v]<range[1][v]) {
                range[0][v]++;
                range[1][v]--;
            } else {
                range[1][v]++;
                range[0][v]--;
            }
        }
        if (level) {
            const int axis = Aggregate_Axis(level);
            range[0][axis] = (range[0][axis]+range[1][axis])/2;
            range[1][axis] = range[0][axis];
        }
        return range;
    }
//#####################################################################

}
