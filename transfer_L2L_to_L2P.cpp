#include "transfer_L2L_to_L2P.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Eigen;

void transfer_L2L_to_L2P(oct_node * node){

  	node->L2P_operator = MatrixXd::Zero(node->nDOF,node->rank);

  	for (int iChild=0; iChild<8; iChild++){
		if (node->child[iChild] != NULL){
			node->L2P_operator.block(node->rank_children_cum(iChild),0,node->child[iChild]->rank,node->rank)=node->L2L_operator[iChild];


			// if (node->id==17){
			// 	cout << "node->L2L_operator[iChild].rows()" << node->L2L_operator[iChild].rows() << endl;
			// 	cout << "node->L2L_operator[iChild].cols()" << node->L2L_operator[iChild].cols() << endl;
			// }
		}
	}
	// cout << "node->L2P_operator" << node->L2P_operator  << endl;

	// if (node->id==17){
	// 	cout << "node->L2P_operator.rows()" << node->L2P_operator.rows() << endl;
	// 	cout << "node->L2P_operator.cols()" << node->L2P_operator.cols() << endl;
	// }
}