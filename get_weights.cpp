#include "get_weights.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

#ifdef DISABLE_JACOBISVD
#include "lapacksvd.hpp"
#endif

using namespace std;
using namespace Eigen;

void get_weights(oct_node * node, int &curr_level, double &epsilon){

    const double machine_eps = std::numeric_limits<double>::epsilon();


  	int nRows=node->nDOF;
  	int nColumns=0;
  	// INTERACTION LIST
    for (int lInteract=0; lInteract < node->Ilist.m; lInteract++){
    	nColumns+=node->Ilist(lInteract)->rank;
    }
    if(curr_level>2){
    	nColumns+=node->parent->rank;
    }



    // INCOMING EDGES
    MatrixXd Incoming_concat=MatrixXd::Zero(nRows,nColumns);
    int pos_block=0;
    // INTERACTION LIST
    for (int lInteract=0; lInteract < node->Ilist.m; lInteract++){

      // also reverse
      int index_lNgbr_to_iSelf=0;
      for (int i=0; i < node->Ilist(lInteract)->Ilist.m;i++){
        if (node->Ilist(lInteract)->Ilist(i)->id == node->id)    
        {
          index_lNgbr_to_iSelf=i; // index l to i
          break;
        }
      }
      Incoming_concat.block(0,pos_block,nRows,node->Ilist(lInteract)->rank) = node->L2P_operator*node->Ilist(lInteract)->M2L_operator[index_lNgbr_to_iSelf];
      pos_block+=node->Ilist(lInteract)->rank;  


    }


    if(curr_level>2){
    	int id_node_local=0;
      	for (int iChild=0;iChild<8;iChild++){
       		if (node->parent->child[iChild] != NULL){
        		if (node->parent->child[iChild]->id == node->id){
           			id_node_local=iChild;
          		}
        	}
      	}
    	Incoming_concat.block(0,pos_block,nRows,node->parent->rank) = node->L2P_operator*node->parent->L2L_operator[id_node_local];
    }

    
    // SVD OF Incoming_concat
#ifdef DISABLE_JACOBISVD
    LapackSVD<MatrixXd> svd_Incoming_concat(Incoming_concat, ComputeThinU | ComputeThinV);
#else
    JacobiSVD<MatrixXd> svd_Incoming_concat(Incoming_concat, ComputeThinU | ComputeThinV);
#endif
    VectorXd Sigma_Incoming_concat = svd_Incoming_concat.singularValues();


    int iSV_Incoming_concat=0;
    while (Sigma_Incoming_concat(iSV_Incoming_concat) > epsilon*Sigma_Incoming_concat(0) && Sigma_Incoming_concat(iSV_Incoming_concat) > machine_eps && iSV_Incoming_concat < (Sigma_Incoming_concat.rows()-1) ){
    	iSV_Incoming_concat++; 
    }
    if (Sigma_Incoming_concat(iSV_Incoming_concat) > epsilon*Sigma_Incoming_concat(0) && Sigma_Incoming_concat(iSV_Incoming_concat) > machine_eps){
      	iSV_Incoming_concat=Sigma_Incoming_concat.rows();
    }





    // OUTGOING EDGES
    MatrixXd Outgoing_concat=MatrixXd::Zero(nColumns,nRows);
    pos_block=0;
    // INTERACTION LIST
    for (int lInteract=0; lInteract < node->Ilist.m; lInteract++){
      Outgoing_concat.block(pos_block,0,node->Ilist(lInteract)->rank,nRows) = node->M2L_operator[lInteract]*node->P2M_operator;

      pos_block+=node->Ilist(lInteract)->rank;
    }
    if(curr_level>2){

    	int id_node_local=0;
      	for (int iChild=0;iChild<8;iChild++){
       		if (node->parent->child[iChild] != NULL){
        		if (node->parent->child[iChild]->id == node->id){
           			id_node_local=iChild;
          		}
        	}
      	}
      	Outgoing_concat.block(pos_block,0,node->parent->rank,nRows) = node->parent->M2M_operator[id_node_local]*node->P2M_operator;
    }



     // SVD OF Outgoing_concat
#ifdef DISABLE_JACOBISVD
    LapackSVD<MatrixXd> svd_Outgoing_concat(Outgoing_concat, ComputeThinU | ComputeThinV);
#else
    JacobiSVD<MatrixXd> svd_Outgoing_concat(Outgoing_concat, ComputeThinU | ComputeThinV);
#endif
    VectorXd Sigma_Outgoing_concat = svd_Outgoing_concat.singularValues();

    int iSV_Outgoing_concat=0;
    while (Sigma_Outgoing_concat(iSV_Outgoing_concat) > epsilon*Sigma_Outgoing_concat(0) && Sigma_Outgoing_concat(iSV_Outgoing_concat) > machine_eps && iSV_Outgoing_concat < (Sigma_Outgoing_concat.rows()-1) ){
    
    	iSV_Outgoing_concat++; 
    }
    if (Sigma_Outgoing_concat(iSV_Outgoing_concat) > epsilon*Sigma_Outgoing_concat(0) && Sigma_Outgoing_concat(iSV_Outgoing_concat) > machine_eps){
      	iSV_Outgoing_concat=Sigma_Outgoing_concat.rows();
    }



    // MAKE SURE Uj AND Vj HAVE THE SAME RANK
    if (iSV_Incoming_concat != iSV_Outgoing_concat){
      iSV_Incoming_concat=min(node->rank,max(iSV_Outgoing_concat,iSV_Incoming_concat));
      iSV_Outgoing_concat = iSV_Incoming_concat;
    }


    MatrixXd Sigma_diag_Incoming_concat=MatrixXd::Zero(iSV_Incoming_concat,iSV_Incoming_concat);
    for (int index=0;index<iSV_Incoming_concat;index++){
      Sigma_diag_Incoming_concat(index,index)=Sigma_Incoming_concat(index);
    }

    MatrixXd Sigma_diag_Outgoing_concat=MatrixXd::Zero(iSV_Outgoing_concat,iSV_Outgoing_concat);
    for (int index=0;index<iSV_Outgoing_concat;index++){
      Sigma_diag_Outgoing_concat(index,index)=Sigma_Outgoing_concat(index);
    }


    MatrixXd R_L2P_operator = Sigma_diag_Incoming_concat*svd_Incoming_concat.matrixV().block(0,0,nColumns,iSV_Incoming_concat).transpose();
    MatrixXd R_P2M_operator = svd_Outgoing_concat.matrixU().block(0,0,nColumns,iSV_Outgoing_concat)*Sigma_diag_Outgoing_concat;



    // UPDATE ALL OPERATORS
    // NEIGHBOUR LIST
    for (int lNgbr=0; lNgbr < node->ngbr.m; lNgbr++){

      int nRows_temp=node->M2L_operator[node->Ilist.m+lNgbr].rows();
      node->M2L_operator[node->Ilist.m+lNgbr].resize(nRows_temp,iSV_Outgoing_concat);
      node->M2L_operator[node->Ilist.m+lNgbr]=MatrixXd::Zero(nRows_temp,iSV_Outgoing_concat);


      // also reverse
      int index_lNgbr_to_iSelf=0;
      for (int i=0; i < node->ngbr(lNgbr)->ngbr.m;i++){
        if (node->ngbr(lNgbr)->ngbr(i)->id == node->id)    
        {
          index_lNgbr_to_iSelf=i; // index l to i
          break;
        }
      }
      int nCols_temp=node->ngbr(lNgbr)->M2L_operator[node->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_iSelf].cols();
      node->ngbr(lNgbr)->M2L_operator[node->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_iSelf].resize(iSV_Incoming_concat,nCols_temp);
      node->ngbr(lNgbr)->M2L_operator[node->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_iSelf]=MatrixXd::Zero(iSV_Incoming_concat,nCols_temp);

    }

    int pos_block_incoming=0;
    int pos_block_outgoing=0;
    // INTERACTION LIST
    for (int lInteract=0; lInteract < node->Ilist.m; lInteract++){
      node->M2L_operator[lInteract].resize(node->Ilist(lInteract)->rank,iSV_Outgoing_concat);
      node->M2L_operator[lInteract] =  R_P2M_operator.block(pos_block_outgoing,0,node->Ilist(lInteract)->rank,iSV_Outgoing_concat);
      pos_block_outgoing+=node->Ilist(lInteract)->rank;

      // also reverse
      int index_lNgbr_to_iSelf=0;
      for (int i=0; i < node->Ilist(lInteract)->Ilist.m;i++){
        if (node->Ilist(lInteract)->Ilist(i)->id == node->id)    
        {
          index_lNgbr_to_iSelf=i; // index l to i
          break;
        }
      }
      node->Ilist(lInteract)->M2L_operator[index_lNgbr_to_iSelf].resize(iSV_Incoming_concat,node->Ilist(lInteract)->rank);
      node->Ilist(lInteract)->M2L_operator[index_lNgbr_to_iSelf] = R_L2P_operator.block(0,pos_block_incoming,iSV_Incoming_concat,node->Ilist(lInteract)->rank);
      pos_block_incoming+=node->Ilist(lInteract)->rank;
    }

    
    if(curr_level>2){

      int id_node_local=0;
      for (int iChild=0;iChild<8;iChild++){
       if (node->parent->child[iChild] != NULL){
         if (node->parent->child[iChild]->id == node->id){
           id_node_local=iChild;
          }
        }
      }


      node->parent->M2M_operator[id_node_local]=R_P2M_operator.block(pos_block_outgoing,0,node->parent->rank,iSV_Outgoing_concat);
      node->parent->L2L_operator[id_node_local]=R_L2P_operator.block(0,pos_block_incoming,iSV_Incoming_concat,node->parent->rank);


    }




    // UPDATE Uj and Vj
    node->L2P_operator=svd_Incoming_concat.matrixU().block(0,0,nRows,iSV_Incoming_concat);
    node->P2M_operator=svd_Outgoing_concat.matrixV().block(0,0,nRows,iSV_Outgoing_concat).transpose();

    
    

    // ASSIGN WEIGHTS TO Uj and Vj
    node->Sigma_L2P=Sigma_diag_Incoming_concat;
    node->Sigma_P2M=Sigma_diag_Outgoing_concat;
	

    // UPDATE RANK 
    node->rank=iSV_Incoming_concat;




	
}
