#ifndef _INIT_M2L_ENTRIES_
#define _INIT_M2L_ENTRIES_

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <cstring>
#include <list>
#include <cassert>
#include <sys/time.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <vector>


#include "container.h"
#include "ifmm.h"

using namespace std;
using namespace Eigen;

class vec3;
struct oct_node;

void init_M2L_entries(int n, int l, double S, ivec & Kidx, kernel * kfun, dvec & M2LOp);

inline int M2L_index_map(vec3 & i);

#endif