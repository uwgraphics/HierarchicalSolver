//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T2> EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
EXTRAPOLATION_HIGHER_ORDER()
//    :grid(grid_input),ghost(ghost_input),phi(phi_input),dt((T).5,(T).33333,(T).25)
{
//    dt*=grid.dX.Min();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class T2> EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
~EXTRAPOLATION_HIGHER_ORDER()
{
}
//#####################################################################
// Function Fill_Level
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Fill_Level(const GRID<TV>& grid,const T_LEVELSET& phi,int ghost,MAPPING& m,ARRAY<TV>& normal,ARRAY<VECTOR<STENCIL,TV::m> >& stencil,int order,T distance)
{
    ARRAY<TV_INT> inside;
    int max_ghost=ghost+3;//2*order+2;
    m.node_to_index.Resize(grid.Domain_Indices(max_ghost)); // Need extra rings for the sentinals
    m.node_to_index.Fill(0);
    m.index_to_node.Append(TV_INT::All_Ones_Vector()*INT_MAX); // First index is the "outside" index.
    normal.Append(TV::All_Ones_Vector()*(typename TV::SCALAR)INT_MAX);

    // Cells that must be solved for normally.
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,max_ghost);it.Valid();it.Next()){
        T p=phi.Phi(it.Location());
        if(p<=0){
            if(p>-(T)2.5*grid.dX.Max()){inside.Append(it.index);m.node_to_index(it.index)=-1;}} // Register two levels to prevent the closure below from leaking inside.
        else if(p<=(distance+(T)1e-4)*grid.dX.Max()) m.node_to_index(it.index)=m.index_to_node.Append(it.index);}

    // Ensure that two upwind nodes are being solved for.
    for(int i=2;i<=m.index_to_node.m;i++){
        TV N=phi.Normal(grid.X(m.index_to_node(i)));
        normal.Append(N);
        for(int d=1;d<=TV::m;d++){
            int s=N(d)<0?1:-1;
            TV_INT ind=m.index_to_node(i);
            for(int j=1;j<=2;j++){
                ind(d)+=s;
                if(!grid.Inside_Domain(ind,max_ghost))continue;
                int& k=m.node_to_index(ind);
                if(!k) k=m.index_to_node.Append(ind);}}}
    m.max_solve_index(1)=m.index_to_node.m;

    // Register additional cells inside for derivatives.
    for(int o=1,i=1;o<=order*2+1;o++){
        int previous=m.index_to_node.m,mx=inside.m;
        for(;i<=mx;i++){
            bool added=false;
            for(int k=1;k<=TV::m*2;k++){
                TV_INT ind=grid.Node_Neighbor(inside(i),k);
                if(!grid.Inside_Domain(ind,max_ghost))continue;
                int& n=m.node_to_index(ind);
                if(!n){n=-1;inside.Append(ind);}
                else if(n>0 && n<=previous && !added){
                    m.node_to_index(inside(i))=m.index_to_node.Append(inside(i));
                    normal.Append(phi.Normal(grid.X(inside(i))));
                    added=true;}}
            if(!added) inside.Append(TV_INT(inside(i)));}  // Use a copy to avoid but on resize
        m.max_solve_index(o+1)=m.index_to_node.m;}

    m.ghost_node.Resize(m.index_to_node.m);
    m.ghost_node(1)=false;
    for(int i=2;i<=m.max_solve_index(2*order+2);i++) m.ghost_node(i)=!grid.Inside_Domain(m.index_to_node(i));

    // Register sentinal layer and precompute stencils.
    stencil.Resize(m.max_solve_index(order));
    for(int i=2;i<=m.max_solve_index(order);i++)if(!m.ghost_node(i)){
        TV N=normal(i);
        for(int d=1;d<=TV::m;d++){
            int s=N(d)<0?1:-1;
            STENCIL& st=stencil(i)(d);
            st.scale=s*N(d)*grid.one_over_dX(d);
            TV_INT ind=m.index_to_node(i);
            for(int j=1;j<=3;j++){
                st.nodes(j+1)=m.node_to_index(ind);
                ind(d)+=s;}
            ind(d)-=4*s;
            int& n=m.node_to_index(ind);
            if(!n) n=1;
            st.nodes(1)=n;}}
}
//#####################################################################
// Function Fill_un
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Fill_un(const GRID<TV>& grid,const MAPPING& m,const TV& one_over_dx,const ARRAY<TV>& normal,const ARRAY<T2>& x,ARRAY<T2>& xn,int o,int mo,BOUNDARY_UNIFORM<GRID<TV>,T2>& boundary)
{
    RANGE<TV_INT> safe_domain=m.node_to_index.domain.Thickened(-1);
    xn.Resize(m.max_solve_index(mo*2-o+2));
    for(int i=m.max_solve_index(o+1)+1;i<=m.max_solve_index(mo*2-o+2);i++){
        const TV_INT& index=m.index_to_node(i);
        if(!safe_domain.Lazy_Inside(index)) continue;
        T2 v=0;
        for(int d=1;d<=TV::m;d++){
            TV_INT a=index,b=index;a(d)--;b(d)++;
            v+=(x(m.node_to_index(b))-x(m.node_to_index(a)))*(T).5*one_over_dx(d)*normal(i)(d);}
        xn(i)=v;}
    Fill_Ghost_Region(grid,m,xn,boundary);
}
//#####################################################################
// Function Constant_Extrapolate_FE
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_FE(const GRID<TV>& grid,const MAPPING& m,const ARRAY<VECTOR<STENCIL,TV::m> >& stencil,const ARRAY<T2>& u,ARRAY<T2>& y,const ARRAY<T2>* z,int o,T dt,T alpha,BOUNDARY_UNIFORM<GRID<TV>,T2>& boundary)
{
    for(int i=2;i<=m.max_solve_index(o);i++)if(!m.ghost_node(i)){
        T2 dot=0,zi=z?(*z)(i):0,a=0;
        for(int d=1;d<=TV::m;d++){
            // Second order ENO.
            const STENCIL& s=stencil(i)(d);
            PHYSBAM_ASSERT(s.nodes(1)<=u.m && s.nodes(2)<=u.m && s.nodes(3)<=u.m && s.nodes(4)<=u.m);
            VECTOR<T2,4> f(u(s.nodes(1)),u(s.nodes(2)),u(s.nodes(3)),u(s.nodes(4)));
            VECTOR<T2,3> df(f(2)-f(1),f(3)-f(2),f(4)-f(3));
            VECTOR<T2,2> ddf(df(2)-df(1),df(3)-df(2));
            if(abs(ddf(1))<abs(ddf(2))) a=(T).5*(f(3)-f(1));
            else a=df(2)-(T).5*ddf(2);
            dot+=a*s.scale;}
        y(i)+=alpha*(u(i)-y(i)-dt*(dot-zi));}
    Fill_Ghost_Region(grid,m,y,boundary);
}
//#####################################################################
// Function Constant_Extrapolate_RK2
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_RK2(const GRID<TV>& grid,const MAPPING& m,const ARRAY<VECTOR<STENCIL,TV::m> >& stencil,ARRAY<T2>& u,const ARRAY<T2>* z,ARRAY<T2>& tmp,int o,T dt,BOUNDARY_UNIFORM<GRID<TV>,T2>& boundary)
{
    tmp.Resize(u.m);
    for(int i=m.max_solve_index(o)+1;i<=u.m;i++) tmp(i)=u(i);
    Extrapolate_FE(grid,m,stencil,u,tmp,z,o,dt,1,boundary);
    Extrapolate_FE(grid,m,stencil,tmp,u,z,o,dt,(T).5,boundary);
}
//#####################################################################
// Function Quadratic_Extrapolate
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_Node(const GRID<TV>& grid,const T_LEVELSET& phi,int ghost,ARRAYS_ND_BASE<VECTOR<T2,TV::m> >& u,int iterations,int order,T distance,BOUNDARY_UNIFORM<GRID<TV>,T2>& boundary)
{
    PHYSBAM_ASSERT(order>=1 && order<=3);
    PHYSBAM_ASSERT(!grid.Is_MAC_Grid());
    T dt=grid.dX.Max()/(TV::m+1);
    MAPPING m;
    ARRAY<TV> normal;
    ARRAY<VECTOR<STENCIL,TV::m> > stencil;
    Fill_Level(grid,phi,ghost,m,normal,stencil,order,distance);
    ARRAY<T2> du[3];
    du[0].Resize(m.max_solve_index(2*order+2));
    ARRAY<T2> tmp(m.max_solve_index(2*order+2));
    for(int i=m.max_solve_index(1)+1;i<=m.max_solve_index(2*order+2);i++) du[0](i)=u(m.index_to_node(i));
    for(int o=1;o<order;o++) Fill_un(grid,m,grid.one_over_dX,normal,du[o-1],du[o],o,order,boundary);
    for(int i=0;i<order;i++) du[i](1)=FLT_MAX/100; // Sentinal values for ENO.
    tmp(1)=FLT_MAX/100;
    for(int o=order;o>0;o--) for(int i=1;i<=iterations;i++) Extrapolate_RK2(grid,m,stencil,du[o-1],(o!=order?&du[o]:0),tmp,o,dt,boundary);
    for(int i=2;i<=m.max_solve_index(1);i++) u(m.index_to_node(i))=du[0](i);
}
//#####################################################################
// Function Fill_Ghost_Region
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Fill_Ghost_Region(const GRID<TV>& grid,const MAPPING& m,ARRAY<T2>& target,BOUNDARY_UNIFORM<GRID<TV>,T2>& boundary)
{
    static ARRAY<T2,TV_INT> tmp;
    tmp.Resize(m.node_to_index.domain);
    for(int i=2;i<=target.m;i++) tmp(m.index_to_node(i))=target(i);
    boundary.Fill_Ghost_Cells_Cell(grid,tmp,tmp,0,m.node_to_index.Number_Of_Ghost_Cells());
    for(int i=2;i<=target.m;i++) target(i)=tmp(m.index_to_node(i));
}
//#####################################################################
// Function Quadratic_Extrapolate
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_Cell(const GRID<TV>& grid,const T_LEVELSET& phi,int ghost,ARRAYS_ND_BASE<VECTOR<T2,TV::m> >& u,int iterations,int order,T distance,BOUNDARY_UNIFORM<GRID<TV>,T2>& boundary)
{
    GRID<TV> node_grid(grid.Get_Regular_Grid_At_MAC_Positions());
    Extrapolate_Node(node_grid,phi,ghost,u,iterations,order,distance,boundary);
}
//#####################################################################
// Function Quadratic_Extrapolate
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_Face(const GRID<TV>& grid,const T_LEVELSET& phi,int ghost,ARRAY<T2,FACE_INDEX<TV::m> >& u,int iterations,int order,T distance,BOUNDARY_UNIFORM<GRID<TV>,T2>& boundary)
{
    for(int i=1;i<=TV::m;i++){
        GRID<TV> node_grid(grid.Get_Axis_X_Face_Grid(i));
        Extrapolate_Node(node_grid,phi,ghost,u.Component(i),iterations,order,distance,boundary);}
}
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<float,1>,float>;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<float,2>,float>;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<float,3>,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<double,1>,double>;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<double,2>,double>;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<double,3>,double>;
#endif
