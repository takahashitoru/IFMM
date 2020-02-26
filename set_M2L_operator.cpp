#include "set_M2L_operator.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Eigen;

void set_M2L_operator(oct_node * node, int &n){
  int n3=n*n*n;

  for (int iInteract = 0; iInteract < node->Ilist.m; ++iInteract) {      
        node->M2L_operator.push_back(MatrixXd::Zero(n3,n3));
        MatrixXd temp = Map<MatrixXd>(node->K(iInteract),n3,n3);
        node->M2L_operator[iInteract] = (temp).transpose();

        // cout << "node->M2L_operator[iInteract]" << endl << node->M2L_operator[iInteract] << endl;
        
  }


  for (int iNgbr = 0; iNgbr < node->ngbr.m; ++iNgbr) {
        node->M2L_operator.push_back(MatrixXd::Zero(n3,n3));
  }

}