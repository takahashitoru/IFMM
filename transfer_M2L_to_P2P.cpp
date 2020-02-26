#include "transfer_M2L_to_P2P.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Eigen;

void transfer_M2L_to_P2P(oct_node * node){

	// printf("node->id%d\n",node->id);

	for (int iNgbr = 0; iNgbr < node->ngbr.m; ++iNgbr) {
  		node->P2P_operator.push_back(MatrixXd::Zero(node->ngbr(iNgbr)->nDOF,node->nDOF));

		for (int iChild=0; iChild<8; iChild++){
			if (node->child[iChild] != NULL){
				for (int iNgbr_child =0; iNgbr_child < 8; iNgbr_child++){
					if (node->ngbr(iNgbr)->child[iNgbr_child] != NULL){

						int index_iChild_to_iNgbr_child=node->child[iChild]->ngbr.m;
						for (int i=0; i<node->child[iChild]->ngbr.m;i++){
							if(node->child[iChild]->ngbr(i)->id == node->ngbr(iNgbr)->child[iNgbr_child]->id){
								index_iChild_to_iNgbr_child=i;
								break;
							}
						}

				        bool WellSeparated=false;
						if (index_iChild_to_iNgbr_child==node->child[iChild]->ngbr.m){
				          WellSeparated=true;
						  for (int i=0; i<node->child[iChild]->Ilist.m;i++){
							if(node->child[iChild]->Ilist(i)->id == node->ngbr(iNgbr)->child[iNgbr_child]->id){
								index_iChild_to_iNgbr_child=i;
								break;
				            }
				          }
				      	}

				      	if (WellSeparated){
							node->P2P_operator[iNgbr].block(node->ngbr(iNgbr)->rank_children_cum(iNgbr_child),node->rank_children_cum(iChild),node->ngbr(iNgbr)->child[iNgbr_child]->rank,node->child[iChild]->rank)=node->child[iChild]->M2L_operator[index_iChild_to_iNgbr_child];
				      	}
				      	else{
							node->P2P_operator[iNgbr].block(node->ngbr(iNgbr)->rank_children_cum(iNgbr_child),node->rank_children_cum(iChild),node->ngbr(iNgbr)->child[iNgbr_child]->rank,node->child[iChild]->rank)=node->child[iChild]->M2L_operator[node->child[iChild]->Ilist.m+index_iChild_to_iNgbr_child];
				      	}
					}
				}
			}
		}
		// cout << "node->P2P_operator[iNgbr]" << node->P2P_operator[iNgbr] << endl;
	}
}