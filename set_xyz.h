#ifndef _SET_XYZ_
#define _SET_XYZ_

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

void set_xyz(oct_node * node);

#endif