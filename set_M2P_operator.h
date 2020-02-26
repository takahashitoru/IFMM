#ifndef _SET_M2P_OPERATOR_
#define _SET_M2P_OPERATOR_

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

struct oct_node;

void set_M2P_operator(oct_node * node, int &n, bool &IsLeaf);

#endif