#include "set_nDOF.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

void set_nDOF(oct_node * node, int &n, bool &IsLeaf){
   int n3=n*n*n;

  if (IsLeaf){ // leaf level
    node->rank=n3;
    node->nDOF = node->pnt.m;
  }
  else{
    node->rank=n3;
    
    node->nDOF=0;
    node->rank_children_cum.resize(9);
    node->rank_children_cum(0)=0;

    for (int iChild=0; iChild<8; iChild++){
      if (node->child[iChild] != NULL){
        node->nDOF+=node->child[iChild]->rank;
        node->rank_children_cum(iChild+1)=node->rank_children_cum(iChild)+node->child[iChild]->rank;
      }
      else{
        node->rank_children_cum(iChild+1)=node->rank_children_cum(iChild);
      }
    }
  }
    // printf("node->id, node->nDOF %d %d\n",node->id,node->nDOF);
}