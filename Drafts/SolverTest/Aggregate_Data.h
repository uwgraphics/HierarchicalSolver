//#####################################################################
// Copyright (c) 2018, Eftychios Sifakis, Qisi Wang
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Aggregate_Data_h__
#define __Aggregate_Data_h__

#include <cassert>

namespace Macroblock_Utilities{
//#####################################################################
// Class Aggregate_Data
//#####################################################################
    template <class T_DATA>
    struct Aggregate_Data
    {
      size_t n;
      size_t k;
      T_DATA* A_ii; // Packed
      T_DATA* A_ir; // Dense column major 
      T_DATA* A_rr; // RFP

      inline const T_DATA& operator()(size_t i, size_t j) const {
        assert((0<=j) && (j<=i) && (i<=n));
	if (i<k && j<k) {
	  return A_ii[j*k+i-(j*(j+1))/2];
	} else if (i>=k && j>=k) {
	  i -= k;
	  j -= k;
	  const size_t m = n-k;
	  size_t l=m/2;
	  if(m%2) // N is odd
            if(j<=l)
	      return A_rr[j*m+i];
            else
	      return A_rr[(i-l)*m+j-l-1];
	  else // M is even
            if(j<l)
	      return A_rr[j*(m+1)+i+1];
            else
	      return A_rr[(i-l)*(m+1)+j-l];
	} else {
	  i-=k;
	  //	  const size_t m = n-k;
	  return A_ir[i*k+j];
	}

      }

      inline T_DATA& operator()(size_t i, size_t j) {
        assert((0<=j) && (j<=i) && (i<=n));
	if (i<k && j<k) {
	  return A_ii[j*k+i-(j*(j+1))/2];
	} else if (i>=k && j>=k) {
	  i -= k;
	  j -= k;
	  const size_t m = n-k;
	  size_t l=m/2;
	  if(m%2) // N is odd
            if(j<=l)
	      return A_rr[j*m+i];
            else
	      return A_rr[(i-l)*m+j-l-1];
	  else // M is even
            if(j<l)
	      return A_rr[j*(m+1)+i+1];
            else
	      return A_rr[(i-l)*(m+1)+j-l];
	} else {
	  i-=k;
	  //	  const size_t m = n-k;
	  return A_ir[i*k+j];
	}

      }


    };
}
#endif//__Aggregate_Data_h__
