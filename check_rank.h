#ifndef _CHECK_RANK_
#define _CHECK_RANK_

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

void check_rank(oct_node * node,int &curr_level, double &epsilon);

#if defined(PARA_CHECK_RANK)
void check_rank_self_setup(oct_node *node, double &epsilon);
void check_rank_self(oct_node * node, int &curr_level);
#endif

#endif
