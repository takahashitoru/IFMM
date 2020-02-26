#include "set_P2L_operator.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

void set_P2L_operator(oct_node * node, int &n, bool &IsLeaf){
  int n3=n*n*n;

  if (IsLeaf){ // leaf level
    for (int iNgbr = 0; iNgbr < node->ngbr.m; ++iNgbr) {
      // set_P2L_operator
      node->P2L_operator.push_back(MatrixXd::Zero(n3,node->pnt.m));
    }
  }
  else{
     for (int iNgbr = 0; iNgbr < node->ngbr.m; ++iNgbr) {
      // set_P2L_operator
      node->P2L_operator.push_back(MatrixXd::Zero(n3,node->nDOF));
    }
  }
}