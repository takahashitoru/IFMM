#include "transfer_RHS.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Eigen;

void transfer_RHS(oct_node * node){

	// cout << "Transfer RHS" << endl;

	node->RHS_leaf=VectorXd::Zero(node->nDOF);

	for (int iChild=0; iChild<8; iChild++){
		if (node->child[iChild] != NULL){
			node->RHS_leaf.block(node->rank_children_cum(iChild),0,node->child[iChild]->rank,1) = node->child[iChild]->RHS;
		}
	}
}