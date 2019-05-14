//#####################################################################
// Copyright 2009, Eftychios Sifakis, Yongning Zhu, Nathan Mitchell
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARBITRARY_RANGE_ITERATOR
//#####################################################################
#ifndef __ARBITRARY_RANGE_ITERATOR__
#define __ARBITRARY_RANGE_ITERATOR__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>

namespace PhysBAM{

template<int d,int stride=1> class ARBITRARY_RANGE_ITERATOR;

template<int d> 
class ARBITRARY_RANGE_ITERATOR<d,1>
{
    typedef VECTOR<int,d> TV_INT;
    typedef RANGE<TV_INT> T_RANGE;
    typedef VECTOR<int,d> TV_DIR;

    const T_RANGE& range;
    TV_INT index;

public:
    TV_DIR direction;
    


    ARBITRARY_RANGE_ITERATOR(const T_RANGE& range_input)
        :range(range_input)
    {
        for( int i=1; i<=d; i++)
            direction(i) = range.min_corner(i) < range.max_corner(i) ? 1 : -1;
        Reset();
    }

    void Reset()
    {index=range.min_corner;}

    bool Valid() const
    {
        return (direction.x*index.x)<=(direction.x*range.max_corner.x);
    }

    void Next()
    {for(int i=d;i>=1;i--) if((direction(i)*index(i))<(direction(i)*range.max_corner(i)) || i==1){index(i)+=direction(i);return;} else index(i)=range.min_corner(i);}

    const TV_INT& Index()
    {return index;}
};

template<int d,int stride>
class ARBITRARY_RANGE_ITERATOR
{
    STATIC_ASSERT((stride!=1));
    typedef VECTOR<int,d> TV_INT;
    typedef RANGE<TV_INT> T_RANGE;
    typedef VECTOR<int,d> TV_DIR;

    const T_RANGE& range;
    TV_INT index;
public:

    TV_DIR direction;

    ARBITRARY_RANGE_ITERATOR(const T_RANGE& range_input)
        :range(range_input)
    {
        for( int i=1; i<=d; i++)
            direction(i) = range.min_corner(i) < range.max_corner(i) ? 1 : -1;
        Reset();
    }

    void Reset()
    {index=range.min_corner;}

    bool Valid() const          
    {
        return (direction.x*index.x)<=(direction.x*range.max_corner.x);
    }
    
    void Next()
    {for(int i=d;i>=1;i--) if((direction(i)*(index(i)+(direction(i)*stride)))<=(direction(i)*range.max_corner(i)) || i==1){index(i)+=direction(i)*stride;return;} else index(i)=range.min_corner(i);}

    const TV_INT& Index()
    {return index;}
};

//#####################################################################
}
#endif
