#ifndef _TRANSFER_L2L_TO_L2P_
#define _TRANSFER_L2L_TO_L2P_

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

// class vec3;
struct oct_node;

void transfer_L2L_to_L2P(oct_node * node);

#endif