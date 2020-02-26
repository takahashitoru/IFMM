#ifndef _INIT_M2L_BLOCK_
#define _INIT_M2L_BLOCK_

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
class kernel;

void init_M2L_block(double b_size, dvec & nodes, vec3 & c_f, vec3 & c_s, dmat & Kmat, kernel * kfun);

#endif
