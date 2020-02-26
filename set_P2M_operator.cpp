#include "set_P2M_operator.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Eigen;

void set_P2M_operator(MatrixXd & P2M_operator,dmat & Sr_x, dmat & Sr_y, dmat & Sr_z){
  const int nx = Sr_x.m;
  const int r = Sr_x.n;

  int k = 0;
  for (int l1 = 0; l1 < r; l1++) {
    for (int l2 = 0; l2 < r; l2++) {
      for (int l3 = 0; l3 < r; l3++, k++) {
        for (int j = 0; j < nx; j++) {         
          P2M_operator(k,j)=Sr_x(j, l1) * Sr_y(j, l2) * Sr_z(j, l3);
        }
      }
    }
  }
};