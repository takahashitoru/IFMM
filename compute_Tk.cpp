#include "compute_Tk.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

void compute_Tk(dvec & nodes, dmat & T) {
  const int r = nodes.m;
  

  	if (r>1){
	  for (int i = 0; i < r; ++i) {
	    const double x = nodes(i);
	    T(0, i) = 1;
	    T(1, i) = x;
	    for (int k = 2; k < r; k++)
	      T(k, i) = 2.0 * x * T(k - 1, i) - T(k - 2, i);
	  }
	}
	else{
		T(0,0)=1;
	}


}