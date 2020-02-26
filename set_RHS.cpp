#include "set_RHS.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

void set_RHS(oct_node * node, bool &IsLeaf){

  if (IsLeaf){
    node->RHS_leaf=Map<VectorXd>(&node->sigma(0),node->pnt.m);
    node->RHS=VectorXd::Zero(node->rank);
  }
  else{
    node->RHS=VectorXd::Zero(node->rank);
  }

}