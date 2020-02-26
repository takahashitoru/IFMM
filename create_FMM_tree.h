#ifndef _CREATE_FMM_TREE_
#define _CREATE_FMM_TREE_

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

void create_FMM_tree(int levels, int idx, vec3 & point, oct_node * node);

#endif