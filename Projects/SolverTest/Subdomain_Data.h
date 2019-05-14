//#####################################################################
// Copyright (c) 2018, Eftychios Sifakis, Qisi Wang
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Subdomain_Data_h__
#define __Subdomain_Data_h__

namespace Macroblock_Utilities{
//#####################################################################
// Class Subdomain_Data
//#####################################################################
    template <class T_DATA>
    struct Subdomain_Data
    {
        /* T_DATA* A_diagonal; */
        /* T_DATA* A_lower; */
        T_DATA* L_diagonal;
        T_DATA* L_lower;
    };
}
#endif//__Subdomain_Data_h__
