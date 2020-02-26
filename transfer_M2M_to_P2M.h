#ifndef _TRANSFER_M2M_TO_P2M_
#define _TRANSFER_M2M_TO_P2M_

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

void transfer_M2M_to_P2M(oct_node * node);

#endif