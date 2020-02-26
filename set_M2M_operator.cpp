#include "set_M2M_operator.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Eigen;

void set_M2M_operator(oct_node * node,int &n,dvec & Sr){

  // cout << "test_M2M" << endl;

  int n2 = n*n;
  int n3 = n*n*n;

  node->M2M_operator.resize(8);

  for (int iChild=0;iChild<8;iChild++){
    node->M2M_operator[iChild]=MatrixXd::Zero(n3,n3);
    // node->M2L_operator.push_back(MatrixXd::Zero(n3,n3));

    int id_node_local=iChild;

    int mirr_x = (id_node_local < 4) ? 0:1;
    int mirr_y = (id_node_local <2 || (id_node_local>3 && id_node_local<6)) ? 0:1;
    int mirr_z = (id_node_local % 2 == 0) ? 0:1;


    MatrixXd Sz_temp=MatrixXd::Zero(n3,n3);
    MatrixXd Sy_temp=MatrixXd::Zero(n3,n3);
    MatrixXd Sx_temp=MatrixXd::Zero(n3,n3);

    // Transform z-coordinate
    if (mirr_z==0){
    for (int index=0;index<n3;index+=n){
      for (int l1 = 0; l1 < n; l1++){
        for (int l2 = 0; l2 < n; l2++){
          Sz_temp(index+l1,index+l2)=Sr(l1 * n + l2);
        }
      }
    }
    }
    else{
    for (int index=0;index<n3;index+=n){
      for (int l1 = 0; l1 < n; l1++){
        for (int l2 = 0; l2 < n; l2++){
          int temp_1=(n-1)-l1;
          int temp_2=(n-1)-l2;
          Sz_temp(index+l1,index+l2)=Sr(temp_1 * n + temp_2);
        }
      }
    }
    }

    // Transform y-coordinate
    if (mirr_y==0){
    for (int l1 = 0; l1 < n3; l1 += n2) {
      for (int l3 = l1; l3 < l1 + n; l3++){
        for (int ind_i = 0; ind_i < n ; ind_i++){
          for (int ind_j =0; ind_j < n ; ind_j++){
            Sy_temp(l3 + ind_i*n,l3 + ind_j*n)=Sr(ind_i * n + ind_j);
          }
        }
      }
    }
    }
    else{
    for (int l1 = 0; l1 < n3; l1 += n2) {
      for (int l3 = l1; l3 < l1 + n; l3++){
        for (int ind_i = 0; ind_i < n ; ind_i++){
          for (int ind_j =0; ind_j < n ; ind_j++){
              int temp_1=(n-1)-ind_i;
              int temp_2=(n-1)-ind_j;
              Sy_temp(l3 + ind_i*n,l3 + ind_j*n)=Sr(temp_1 * n + temp_2);
          }
        }
      }
    }
    }

    // Transform x-coordinate
    if (mirr_x==0){
    for (int l2 = 0; l2 < n2; l2 += n) {
      for (int l3 = l2; l3 < l2 + n; l3++) {
        for (int ind_i = 0; ind_i < n ; ind_i++){
          for (int ind_j =0; ind_j < n ; ind_j++){
            Sx_temp(l3 + ind_i*n2,l3 + ind_j*n2)=Sr(ind_i * n + ind_j);
          }
        }
      }
    }
    }
    else{
    for (int l2 = 0; l2 < n2; l2 += n) {
      for (int l3 = l2; l3 < l2 + n; l3++) {
        for (int ind_i = 0; ind_i < n ; ind_i++){
          for (int ind_j =0; ind_j < n ; ind_j++){
              int temp_1=(n-1)-ind_i;
              int temp_2=(n-1)-ind_j;
              Sx_temp(l3 + ind_i*n2,l3 + ind_j*n2)=Sr(temp_1 * n + temp_2);
          }
        }
      }
    }
    }
      
    // Form M2M-operator (denoted as P2M here, although this is a non-leaf level)
    node->M2M_operator[iChild]=Sx_temp*(Sy_temp*Sz_temp);

    // MatrixXd temp=Sz_temp.transpose()*(Sy_temp.transpose()*Sx_temp.transpose());
    // node->M2M_operator[iChild]=temp;


}




//   int n2 = n*n;
//   int n3 = n*n*n;

//   int id_node_local;

//   for (int iChild=0;iChild<8;iChild++){
//    if (node->parent->child[iChild] != NULL){
//      if (node->parent->child[iChild]->id == node->id){
//        id_node_local=iChild;
//       }
//     }
//   }

// int mirr_x = (id_node_local < 4) ? 0:1;
// int mirr_y = (id_node_local <2 || (id_node_local>3 && id_node_local<6)) ? 0:1;
// int mirr_z = (id_node_local % 2 == 0) ? 0:1;


// MatrixXd Sz_temp=MatrixXd::Zero(n3,n3);
// MatrixXd Sy_temp=MatrixXd::Zero(n3,n3);
// MatrixXd Sx_temp=MatrixXd::Zero(n3,n3);

// // Transform z-coordinate
// if (mirr_z==0){
//   for (int index=0;index<n3;index+=n){
//     for (int l1 = 0; l1 < n; l1++){
//       for (int l2 = 0; l2 < n; l2++){
//         Sz_temp(index+l1,index+l2)=Sr(l1 * n + l2);
//       }
//     }
//   }
// }
// else{
//   for (int index=0;index<n3;index+=n){
//     for (int l1 = 0; l1 < n; l1++){
//       for (int l2 = 0; l2 < n; l2++){
//         int temp_1=(n-1)-l1;
//         int temp_2=(n-1)-l2;
//         Sz_temp(index+l1,index+l2)=Sr(temp_1 * n + temp_2);
//       }
//     }
//   }
// }

// // Transform y-coordinate
// if (mirr_y==0){
//   for (int l1 = 0; l1 < n3; l1 += n2) {
//     for (int l3 = l1; l3 < l1 + n; l3++){
//       for (int ind_i = 0; ind_i < n ; ind_i++){
//         for (int ind_j =0; ind_j < n ; ind_j++){
//           Sy_temp(l3 + ind_i*n,l3 + ind_j*n)=Sr(ind_i * n + ind_j);
//         }
//       }
//     }
//   }
// }
// else{
//   for (int l1 = 0; l1 < n3; l1 += n2) {
//     for (int l3 = l1; l3 < l1 + n; l3++){
//       for (int ind_i = 0; ind_i < n ; ind_i++){
//         for (int ind_j =0; ind_j < n ; ind_j++){
//             int temp_1=(n-1)-ind_i;
//             int temp_2=(n-1)-ind_j;
//             Sy_temp(l3 + ind_i*n,l3 + ind_j*n)=Sr(temp_1 * n + temp_2);
//         }
//       }
//     }
//   }
// }

// // Transform x-coordinate
// if (mirr_x==0){
//   for (int l2 = 0; l2 < n2; l2 += n) {
//     for (int l3 = l2; l3 < l2 + n; l3++) {
//       for (int ind_i = 0; ind_i < n ; ind_i++){
//         for (int ind_j =0; ind_j < n ; ind_j++){
//           Sx_temp(l3 + ind_i*n2,l3 + ind_j*n2)=Sr(ind_i * n + ind_j);
//         }
//       }
//     }
//   }
// }
// else{
//   for (int l2 = 0; l2 < n2; l2 += n) {
//     for (int l3 = l2; l3 < l2 + n; l3++) {
//       for (int ind_i = 0; ind_i < n ; ind_i++){
//         for (int ind_j =0; ind_j < n ; ind_j++){
//             int temp_1=(n-1)-ind_i;
//             int temp_2=(n-1)-ind_j;
//             Sx_temp(l3 + ind_i*n2,l3 + ind_j*n2)=Sr(temp_1 * n + temp_2);
//         }
//       }
//     }
//   }
// }
    
// // Form M2M-operator (denoted as P2M here, although this is a non-leaf level)
// P2M_operator=Sx_temp*(Sy_temp*Sz_temp);

};