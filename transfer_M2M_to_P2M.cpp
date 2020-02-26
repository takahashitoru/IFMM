#include "transfer_M2M_to_P2M.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Eigen;

void transfer_M2M_to_P2M(oct_node * node){

  	node->P2M_operator = MatrixXd::Zero(node->rank,node->nDOF);

	for (int iChild=0; iChild<8; iChild++){
		if (node->child[iChild] != NULL){
			node->P2M_operator.block(0,node->rank_children_cum(iChild),node->rank,node->child[iChild]->rank) = node->M2M_operator[iChild];


			// if (node->id==17){
				// cout << "node->M2M_operator[iChild].rows()" << node->M2M_operator[iChild].rows() << endl;
				// cout << "node->M2M_operator[iChild].cols()" << node->M2M_operator[iChild].cols() << endl;
			// }

		}
	}
	// cout << "node->P2M_operator" << node->P2M_operator  << endl;

	// if (node->id==17){
	// 	cout << "node->P2M_operator.rows()" << node->P2M_operator.rows() << endl;
	// 	cout << "node->P2M_operator.cols()" << node->P2M_operator.cols() << endl;
	// }

}