#include "set_L2P_operator.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

void set_L2P_operator(MatrixXd & L2P_operator,dmat & Sr_x, dmat & Sr_y, dmat & Sr_z){

  const int nx = Sr_x.m;
  const int r = Sr_x.n;

  for (int i = 0; i < nx; i++) {
    int k = 0;
    for (int l1 = 0; l1 < r; l1++) {
      for (int l2 = 0; l2 < r; l2++) {
        double tmp = Sr_x(i, l1) * Sr_y(i, l2);
        for (int l3 = 0; l3 < r; l3++, k++)
        {
          L2P_operator(i,k)=tmp * Sr_z(i, l3);
        }
      }
    }
  }
}