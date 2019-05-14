//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __VOLUME_PRESERVATION_DEVIATION_H__
#define __VOLUME_PRESERVATION_DEVIATION_H__
template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Volume_Preservation_Deviation
{
static void Run(T_DATA (&M), const T_DATA (&Sigma)[3]);
};
#endif
