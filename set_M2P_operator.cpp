#include "set_M2P_operator.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Eigen;

void set_M2P_operator(oct_node * node, int &n, bool &IsLeaf){

  int n3=n*n*n;

  if (IsLeaf){ // leaf level
    for (int iNgbr = 0; iNgbr < node->ngbr.m; ++iNgbr) {
    // set_M2P_operator
    node->M2P_operator.push_back(MatrixXd::Zero(node->ngbr(iNgbr)->pnt.m,n3));
    }
  }
  else{
    for (int iNgbr = 0; iNgbr < node->ngbr.m; ++iNgbr) {
    // set_M2P_operator
    node->M2P_operator.push_back(MatrixXd::Zero(node->ngbr(iNgbr)->nDOF,n3));
    }
  }
};