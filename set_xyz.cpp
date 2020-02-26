#include "set_xyz.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

void set_xyz(oct_node * node){

  node->x=VectorXd::Zero(node->nDOF);
  node->y=VectorXd::Zero(node->rank);
  node->z=VectorXd::Zero(node->rank);

}