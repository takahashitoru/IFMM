#ifndef _SET_L2P_OPERATOR_
#define _SET_L2P_OPERATOR_

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

using namespace std;
using namespace Eigen;

void set_L2P_operator(MatrixXd & L2P_operator,dmat & Sr_x, dmat & Sr_y, dmat & Sr_z);

#endif