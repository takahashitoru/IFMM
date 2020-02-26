#ifndef _SET_P2P_OPERATOR_
#define _SET_P2P_OPERATOR_

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
class kernel;
// class diag_plus_lr;


void set_P2P_operator(oct_node * node,kernel * kfun, bool &IsLeaf);
// void set_P2P_operator(oct_node * node,diag_plus_lr * dlr_fun, bool &IsLeaf);


#endif