#include "check_rank.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

#ifdef DISABLE_JACOBISVD
#include "lapacksvd.hpp"
#endif

using namespace std;
using namespace Eigen;

void check_rank(oct_node * node, int &curr_level, double &epsilon){


    // printf("curr_level: %d\n",curr_level);
    // printf("node->id, node->rank: %d %d \n",node->id,node->rank);
    // printf("node->id,rank, nDOF: %d %d %d \n",node->id, node->rank, node->nDOF);
    const double machine_eps = std::numeric_limits<double>::epsilon();


    // SVD OF P2M
#ifdef DISABLE_JACOBISVD
    MatrixXd P2M_operator_transpose = node->P2M_operator.transpose();
    LapackSVD<MatrixXd> svd_P2M_operator(P2M_operator_transpose, ComputeThinU | ComputeThinV);
#else
    JacobiSVD<MatrixXd> svd_P2M_operator(node->P2M_operator.transpose(), ComputeThinU | ComputeThinV);
#endif
    VectorXd Sigma_i_P2M = svd_P2M_operator.singularValues();

    // cout << "Sigma_i_P2M" << endl << Sigma_i_P2M << endl;

    int iSV_P2M=0;
    while (Sigma_i_P2M(iSV_P2M) > epsilon*Sigma_i_P2M(0) && Sigma_i_P2M(iSV_P2M) > machine_eps && iSV_P2M < (Sigma_i_P2M.rows()-1) ){
    // while (Sigma_i_P2M(iSV_P2M) > epsilon && Sigma_i_P2M(iSV_P2M) > machine_eps && iSV_P2M < (Sigma_i_P2M.rows()-1) ){
     iSV_P2M++; 
    }
   if (Sigma_i_P2M(iSV_P2M) > epsilon*Sigma_i_P2M(0) && Sigma_i_P2M(iSV_P2M) > machine_eps){
   // if (Sigma_i_P2M(iSV_P2M) > epsilon && Sigma_i_P2M(iSV_P2M) > machine_eps){
      iSV_P2M=Sigma_i_P2M.rows();
    }

    // iSV_P2M=Sigma_i_P2M.rows();

    // SVD OF L2P
#ifdef DISABLE_JACOBISVD
    LapackSVD<MatrixXd> svd_L2P_operator(node->L2P_operator, ComputeThinU | ComputeThinV);
#else
    JacobiSVD<MatrixXd> svd_L2P_operator(node->L2P_operator, ComputeThinU | ComputeThinV);
#endif
    VectorXd Sigma_i_L2P = svd_L2P_operator.singularValues();

    int iSV_L2P=0;
    while (Sigma_i_L2P(iSV_L2P) > epsilon*Sigma_i_L2P(0) && Sigma_i_L2P(iSV_L2P) > machine_eps && iSV_L2P < (Sigma_i_L2P.rows()-1) ){
    // while (Sigma_i_L2P(iSV_L2P) > epsilon && Sigma_i_L2P(iSV_L2P) > machine_eps && iSV_L2P < (Sigma_i_L2P.rows()-1) ){
     iSV_L2P++; 
    }
    if (Sigma_i_L2P(iSV_L2P) > epsilon*Sigma_i_L2P(0) && Sigma_i_L2P(iSV_L2P) > machine_eps){
    // if (Sigma_i_L2P(iSV_L2P) > epsilon && Sigma_i_L2P(iSV_L2P) > machine_eps){
      iSV_L2P=Sigma_i_L2P.rows();
    }

    // iSV_L2P=Sigma_i_L2P.rows();

    // cout << "Sigma_i_L2P" << endl << Sigma_i_L2P << endl;
    
    if (curr_level==2){
      // printf("iSV_P2M, iSV_L2P: %d %d\n",iSV_P2M, iSV_L2P);
      // cout << "Sigma_i_P2M" << endl << Sigma_i_P2M << endl;

    }

    // MAKE SURE Uj AND Vj HAVE THE SAME RANK
    if (iSV_P2M != iSV_L2P){
      iSV_L2P=min(node->rank,max(iSV_P2M,iSV_L2P));
      iSV_P2M = iSV_L2P;

    }

    // printf("iSV_P2M, iSV_L2P: %d %d\n",iSV_P2M, iSV_L2P);


    MatrixXd Sigma_P2M_operator=MatrixXd::Zero(iSV_P2M,iSV_P2M);
    for (int index=0;index<iSV_P2M;index++){
      Sigma_P2M_operator(index,index)=Sigma_i_P2M(index);
    }

    MatrixXd Sigma_L2P_operator=MatrixXd::Zero(iSV_L2P,iSV_L2P);
    for (int index=0;index<iSV_L2P;index++){
      Sigma_L2P_operator(index,index)=Sigma_i_L2P(index);
    }

    MatrixXd U_P2M_operator = svd_P2M_operator.matrixU().block(0,0,node->nDOF,iSV_P2M);
    MatrixXd R_P2M_operator = Sigma_P2M_operator*svd_P2M_operator.matrixV().block(0,0,node->rank,iSV_P2M).transpose();


    MatrixXd U_L2P_operator = svd_L2P_operator.matrixU().block(0,0,node->nDOF,iSV_L2P);
    MatrixXd R_L2P_operator = Sigma_L2P_operator*svd_L2P_operator.matrixV().block(0,0,node->rank,iSV_L2P).transpose();

    // cout << "node->P2M_operator.transpose()" << endl << node->P2M_operator.transpose() << endl;
    // cout << "test" << endl << U_P2M_operator*R_P2M_operator << endl;


    // NEIGHBOUR LIST
    for (int lNgbr=0; lNgbr < node->ngbr.m; lNgbr++){
      node->M2L_operator[node->Ilist.m+lNgbr] = node->M2L_operator[node->Ilist.m+lNgbr]*R_P2M_operator.transpose();

      // also reverse
      int index_lNgbr_to_iSelf=0;
      for (int i=0; i < node->ngbr(lNgbr)->ngbr.m;i++){
        if (node->ngbr(lNgbr)->ngbr(i)->id == node->id)    
        {
          index_lNgbr_to_iSelf=i; // index l to i
          break;
        }
      }
      node->ngbr(lNgbr)->M2L_operator[node->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_iSelf] = R_L2P_operator*node->ngbr(lNgbr)->M2L_operator[node->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_iSelf];
    }

    // INTERACTION LIST
    for (int lInteract=0; lInteract < node->Ilist.m; lInteract++){
      node->M2L_operator[lInteract] =  node->M2L_operator[lInteract]*R_P2M_operator.transpose();

      // also reverse
      int index_lNgbr_to_iSelf=0;
      for (int i=0; i < node->Ilist(lInteract)->Ilist.m;i++){
        if (node->Ilist(lInteract)->Ilist(i)->id == node->id)    
        {
          index_lNgbr_to_iSelf=i; // index l to i
          break;
        }
      }
      node->Ilist(lInteract)->M2L_operator[index_lNgbr_to_iSelf] = R_L2P_operator*node->Ilist(lInteract)->M2L_operator[index_lNgbr_to_iSelf];
    }

    
    if(curr_level>2){
      // UPDATE PARENT INTERACTIONS
      int id_node_local=0;
      for (int iChild=0;iChild<8;iChild++){
       if (node->parent->child[iChild] != NULL){
         if (node->parent->child[iChild]->id == node->id){
           id_node_local=iChild;
          }
        }
      }
      node->parent->M2M_operator[id_node_local]=node->parent->M2M_operator[id_node_local]*R_P2M_operator.transpose();
      node->parent->L2L_operator[id_node_local]=R_L2P_operator*node->parent->L2L_operator[id_node_local];

      // printf("node->parent->id%d\n",node->parent->id);
    }


    node->P2M_operator=U_P2M_operator.transpose();
    node->L2P_operator=U_L2P_operator;

    // UPDATE RANK 
    node->rank=iSV_P2M;


    // printf("node->id,rank, nDOF: %d %d %d \n",node->id, node->rank, node->nDOF);

    if (node->nDOF < node->rank){
      cout << "Check node !" << node->id << endl;

      // NEIGHBOUR LIST
      for (int lNgbr=0; lNgbr < node->ngbr.m; lNgbr++){
        node->M2L_operator[node->Ilist.m+lNgbr] = node->M2L_operator[node->Ilist.m+lNgbr]*node->P2M_operator;

        // also reverse
        int index_lNgbr_to_iSelf=0;
        for (int i=0; i < node->ngbr(lNgbr)->ngbr.m;i++){
          if (node->ngbr(lNgbr)->ngbr(i)->id == node->id)    
          {
            index_lNgbr_to_iSelf=i; // index l to i
            break;
          }
        }
        node->ngbr(lNgbr)->M2L_operator[node->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_iSelf] = node->L2P_operator*node->ngbr(lNgbr)->M2L_operator[node->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_iSelf] ;
      }

      // INTERACTION LIST
      for (int lInteract=0; lInteract < node->Ilist.m; lInteract++){
        node->M2L_operator[lInteract] =  node->M2L_operator[lInteract]*node->P2M_operator;

        // also reverse
        int index_lNgbr_to_iSelf=0;
        for (int i=0; i < node->Ilist(lInteract)->Ilist.m;i++){
          if (node->Ilist(lInteract)->Ilist(i)->id == node->id)    
          {
            index_lNgbr_to_iSelf=i; // index l to i
            break;
          }
        }
        node->Ilist(lInteract)->M2L_operator[index_lNgbr_to_iSelf] = node->L2P_operator*node->Ilist(lInteract)->M2L_operator[index_lNgbr_to_iSelf];
      }

      
      if(curr_level>2){
        // UPDATE PARENT INTERACTIONS
        int id_node_local=0;
        for (int iChild=0;iChild<8;iChild++){
         if (node->parent->child[iChild] != NULL){
           if (node->parent->child[iChild]->id == node->id){
             id_node_local=iChild;
            }
          }
        }
        node->parent->M2M_operator[id_node_local]=node->parent->M2M_operator[id_node_local]*node->P2M_operator;
        node->parent->L2L_operator[id_node_local]=node->L2P_operator*node->parent->L2L_operator[id_node_local];

        // // // printf("node->parent->id%d\n",node->parent->id);
      }


      node->P2M_operator=MatrixXd::Identity(node->nDOF,node->nDOF);
      node->L2P_operator=MatrixXd::Identity(node->nDOF,node->nDOF);

      // UPDATE RANK 
      node->rank=node->nDOF;

    }
    else{

     //  // SVD OF P2M
     //  JacobiSVD<MatrixXd> svd_P2M_operator(node->P2M_operator.transpose(), ComputeThinU | ComputeThinV);
     //  VectorXd Sigma_i_P2M = svd_P2M_operator.singularValues();

     //  // cout << "Sigma_i_P2M" << endl << Sigma_i_P2M << endl;

     //  int iSV_P2M=0;
     //  while (Sigma_i_P2M(iSV_P2M) > epsilon*Sigma_i_P2M(0) && iSV_P2M < (Sigma_i_P2M.rows()-1) ){
     //   iSV_P2M++; 
     //  }
     // if (Sigma_i_P2M(iSV_P2M) > epsilon*Sigma_i_P2M(0)){
     //    iSV_P2M=Sigma_i_P2M.rows();
     //  }

     //  // SVD OF L2P
     //  JacobiSVD<MatrixXd> svd_L2P_operator(node->L2P_operator, ComputeThinU | ComputeThinV);
     //  VectorXd Sigma_i_L2P = svd_L2P_operator.singularValues();

     //  int iSV_L2P=0;
     //  while (Sigma_i_L2P(iSV_L2P) > epsilon*Sigma_i_L2P(0) && iSV_L2P < (Sigma_i_L2P.rows()-1) ){
     //   iSV_L2P++; 
     //  }
     //  if (Sigma_i_L2P(iSV_L2P) > epsilon*Sigma_i_L2P(0)){
     //    iSV_L2P=Sigma_i_L2P.rows();
     //  }

     //  // MAKE SURE Uj AND Vj HAVE THE SAME RANK
     //  if (iSV_P2M != iSV_P2M){
     //    iSV_L2P=min(node->rank,max(iSV_P2M,iSV_L2P));
     //    iSV_P2M = iSV_L2P;
     //  }


     //  MatrixXd Sigma_P2M_operator=MatrixXd::Zero(iSV_P2M,iSV_P2M);
     //  for (int index=0;index<iSV_P2M;index++){
     //    Sigma_P2M_operator(index,index)=Sigma_i_P2M(index);
     //  }

     //  MatrixXd Sigma_L2P_operator=MatrixXd::Zero(iSV_L2P,iSV_L2P);
     //  for (int index=0;index<iSV_L2P;index++){
     //    Sigma_L2P_operator(index,index)=Sigma_i_L2P(index);
     //  }

     //  MatrixXd U_P2M_operator = svd_P2M_operator.matrixU().block(0,0,node->nDOF,iSV_P2M);
     //  MatrixXd R_P2M_operator = Sigma_P2M_operator*svd_P2M_operator.matrixV().block(0,0,node->rank,iSV_P2M).transpose();


     //  MatrixXd U_L2P_operator = svd_L2P_operator.matrixU().block(0,0,node->nDOF,iSV_L2P);
     //  MatrixXd R_L2P_operator = Sigma_L2P_operator*svd_L2P_operator.matrixV().block(0,0,node->rank,iSV_L2P).transpose();

     //  // cout << "node->P2M_operator.transpose()" << endl << node->P2M_operator.transpose() << endl;
     //  // cout << "test" << endl << U_P2M_operator*R_P2M_operator << endl;


     //  // NEIGHBOUR LIST
     //  for (int lNgbr=0; lNgbr < node->ngbr.m; lNgbr++){
     //    node->M2L_operator[node->Ilist.m+lNgbr] = node->M2L_operator[node->Ilist.m+lNgbr]*R_P2M_operator.transpose();

     //    // also reverse
     //    int index_lNgbr_to_iSelf=0;
     //    for (int i=0; i < node->ngbr(lNgbr)->ngbr.m;i++){
     //      if (node->ngbr(lNgbr)->ngbr(i)->id == node->id)    
     //      {
     //        index_lNgbr_to_iSelf=i; // index l to i
     //        break;
     //      }
     //    }
     //    node->ngbr(lNgbr)->M2L_operator[node->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_iSelf] = R_L2P_operator*node->ngbr(lNgbr)->M2L_operator[node->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_iSelf];
     //  }

     //  // INTERACTION LIST
     //  for (int lInteract=0; lInteract < node->Ilist.m; lInteract++){
     //    node->M2L_operator[lInteract] =  node->M2L_operator[lInteract]*R_P2M_operator.transpose();

     //    // also reverse
     //    int index_lNgbr_to_iSelf=0;
     //    for (int i=0; i < node->Ilist(lInteract)->Ilist.m;i++){
     //      if (node->Ilist(lInteract)->Ilist(i)->id == node->id)    
     //      {
     //        index_lNgbr_to_iSelf=i; // index l to i
     //        break;
     //      }
     //    }
     //    node->Ilist(lInteract)->M2L_operator[index_lNgbr_to_iSelf] = R_L2P_operator*node->Ilist(lInteract)->M2L_operator[index_lNgbr_to_iSelf];
     //  }

      
     //  if(curr_level>2){
     //    // UPDATE PARENT INTERACTIONS
     //    int id_node_local;
     //    for (int iChild=0;iChild<8;iChild++){
     //     if (node->parent->child[iChild] != NULL){
     //       if (node->parent->child[iChild]->id == node->id){
     //         id_node_local=iChild;
     //        }
     //      }
     //    }
     //    node->parent->M2M_operator[id_node_local]=node->parent->M2M_operator[id_node_local]*R_P2M_operator.transpose();
     //    node->parent->L2L_operator[id_node_local]=R_L2P_operator*node->parent->L2L_operator[id_node_local];

     //    printf("node->parent->id%d\n",node->parent->id);
     //  }


     //  node->P2M_operator=U_P2M_operator.transpose();
     //  node->L2P_operator=U_L2P_operator;

     //  // UPDATE RANK 
     //  node->rank=iSV_P2M;

    }

    // // // printf("rank: %d \n",node->rank);
    // printf("node->id, node->rank: %d %d \n",node->id,node->rank);


}



#if defined(PARA_CHECK_RANK)

void check_rank_self_setup(oct_node *node, double &epsilon)
{
  const double machine_eps = std::numeric_limits<double>::epsilon();
  
  // SVD OF P2M
  
  //20200108#ifdef DISABLE_JACOBISVD
  //20200108  MatrixXcd P2M_operator_adjoint = node->P2M_operator.adjoint();
  //20200108  LapackSVD<MatrixXcd> svd_P2M_operator(P2M_operator_adjoint, ComputeThinU | ComputeThinV);
  //20200108#else
  //20200108  JacobiSVD<MatrixXcd> svd_P2M_operator(node->P2M_operator.adjoint(), ComputeThinU | ComputeThinV);
  //20200108#endif
#ifdef DISABLE_JACOBISVD
  MatrixXd P2M_operator_transpose = node->P2M_operator.transpose();
  LapackSVD<MatrixXd> svd_P2M_operator(P2M_operator_transpose, ComputeThinU | ComputeThinV);
#else
  JacobiSVD<MatrixXd> svd_P2M_operator(node->P2M_operator.transpose(), ComputeThinU | ComputeThinV);
#endif
  VectorXd Sigma_i_P2M = svd_P2M_operator.singularValues();
  
  int iSV_P2M=0;
  while (Sigma_i_P2M(iSV_P2M) > epsilon*Sigma_i_P2M(0) && Sigma_i_P2M(iSV_P2M) > machine_eps && iSV_P2M < (Sigma_i_P2M.rows()-1) ){
    iSV_P2M++; 
  }
  if (Sigma_i_P2M(iSV_P2M) > epsilon*Sigma_i_P2M(0) && Sigma_i_P2M(iSV_P2M) > machine_eps){
    iSV_P2M=Sigma_i_P2M.rows();
  }
  
  // SVD OF L2P

  //20200108#ifdef DISABLE_JACOBISVD
  //20200108  LapackSVD<MatrixXcd> svd_L2P_operator(node->L2P_operator, ComputeThinU | ComputeThinV);
  //20200108#else
  //20200108  JacobiSVD<MatrixXcd> svd_L2P_operator(node->L2P_operator, ComputeThinU | ComputeThinV);
  //20200108#endif
#ifdef DISABLE_JACOBISVD
  LapackSVD<MatrixXd> svd_L2P_operator(node->L2P_operator, ComputeThinU | ComputeThinV);
#else
  JacobiSVD<MatrixXd> svd_L2P_operator(node->L2P_operator, ComputeThinU | ComputeThinV);
#endif
  VectorXd Sigma_i_L2P = svd_L2P_operator.singularValues();

  int iSV_L2P=0;
  while (Sigma_i_L2P(iSV_L2P) > epsilon*Sigma_i_L2P(0) && Sigma_i_L2P(iSV_L2P) > machine_eps && iSV_L2P < (Sigma_i_L2P.rows()-1) ){
    iSV_L2P++; 
  }
  if (Sigma_i_L2P(iSV_L2P) > epsilon*Sigma_i_L2P(0) && Sigma_i_L2P(iSV_L2P) > machine_eps){
    iSV_L2P=Sigma_i_L2P.rows();
  }

  // MAKE SURE Uj AND Vj HAVE THE SAME RANK

  if (iSV_P2M != iSV_L2P){
    iSV_L2P=min(node->rank,max(iSV_P2M,iSV_L2P));
    iSV_P2M = iSV_L2P;
  }
  
  //20200108  MatrixXcd Sigma_P2M_operator=MatrixXcd::Zero(iSV_P2M, iSV_P2M);
  //20200108  for (int index=0;index<iSV_P2M;index++){
  //20200108    Sigma_P2M_operator(index,index)=Sigma_i_P2M(index);
  //20200108  }
  MatrixXd Sigma_P2M_operator=MatrixXd::Zero(iSV_P2M,iSV_P2M);
  for (int index=0;index<iSV_P2M;index++){
    Sigma_P2M_operator(index,index)=Sigma_i_P2M(index);
  }

  //20200108  MatrixXcd Sigma_L2P_operator=MatrixXcd::Zero(iSV_L2P, iSV_L2P);
  //20200108  for (int index=0;index<iSV_L2P;index++){
  //20200108    Sigma_L2P_operator(index,index)=Sigma_i_L2P(index);
  //20200108  }
  MatrixXd Sigma_L2P_operator=MatrixXd::Zero(iSV_L2P,iSV_L2P);
  for (int index=0;index<iSV_L2P;index++){
    Sigma_L2P_operator(index,index)=Sigma_i_L2P(index);
  }
  

  //20200108  MatrixXcd U_P2M_operator = svd_P2M_operator.matrixU().block(0, 0, node->nDOF, iSV_P2M);
  //20200108  node->R_P2M_operator = Sigma_P2M_operator * svd_P2M_operator.matrixV().block(0, 0, node->rank, iSV_P2M).adjoint();
  //20200108  MatrixXcd U_L2P_operator = svd_L2P_operator.matrixU().block(0, 0, node->nDOF, iSV_L2P);
  //20200108  node->R_L2P_operator = Sigma_L2P_operator * svd_L2P_operator.matrixV().block(0, 0, node->rank, iSV_L2P).adjoint();
  //20200108  node->P2M_operator = U_P2M_operator.adjoint();
  //20200108  node->L2P_operator = U_L2P_operator;
  MatrixXd U_P2M_operator = svd_P2M_operator.matrixU().block(0, 0, node->nDOF, iSV_P2M);
  node->R_P2M_operator = Sigma_P2M_operator * svd_P2M_operator.matrixV().block(0, 0, node->rank, iSV_P2M).transpose();
  MatrixXd U_L2P_operator = svd_L2P_operator.matrixU().block(0, 0, node->nDOF, iSV_L2P);
  node->R_L2P_operator = Sigma_L2P_operator * svd_L2P_operator.matrixV().block(0, 0, node->rank, iSV_L2P).transpose();
  node->P2M_operator = U_P2M_operator.transpose();
  node->L2P_operator = U_L2P_operator;

  // UPDATE RANK 

  node->rank = iSV_P2M;

}


//20200108void check_rank_self(oct_node * node, int &curr_level, ofstream &outfile)
void check_rank_self(oct_node * node, int &curr_level)
{
  // NEIGHBOUR LIST

  for (int lNgbr = 0; lNgbr < node->ngbr.m; lNgbr ++) {
    //20200108    node->M2L_operator[node->Ilist.m + lNgbr] = node->ngbr(lNgbr)->R_L2P_operator * node->M2L_operator[node->Ilist.m + lNgbr] * node->R_P2M_operator.adjoint();
    node->M2L_operator[node->Ilist.m + lNgbr] = node->ngbr(lNgbr)->R_L2P_operator * node->M2L_operator[node->Ilist.m + lNgbr] * node->R_P2M_operator.transpose();
  }

  // INTERACTION LIST
  
  for (int lInteract = 0; lInteract < node->Ilist.m; lInteract ++) {
    //20200108    MatrixXcd temp = node->Ilist(lInteract)->R_L2P_operator * node->M2L_operator[lInteract];
    //20200108    node->M2L_operator[lInteract] = temp * node->R_P2M_operator.adjoint();
    MatrixXd temp = node->Ilist(lInteract)->R_L2P_operator * node->M2L_operator[lInteract];
    node->M2L_operator[lInteract] = temp * node->R_P2M_operator.transpose();
  }

  if (curr_level > 2) {
    // UPDATE PARENT INTERACTIONS
    int id_node_local = 0;
    for (int iChild = 0; iChild < 8; iChild ++) {
      if (node->parent->child[iChild] != NULL) {
	if (node->parent->child[iChild]->id == node->id) {
	  id_node_local=iChild;
	}
      }
    }
    //20200108    node->parent->M2M_operator[id_node_local] = node->parent->M2M_operator[id_node_local] * node->R_P2M_operator.adjoint();
    //20200108    node->parent->L2L_operator[id_node_local] = node->R_L2P_operator * node->parent->L2L_operator[id_node_local];
    node->parent->M2M_operator[id_node_local] = node->parent->M2M_operator[id_node_local] * node->R_P2M_operator.transpose();
    node->parent->L2L_operator[id_node_local] = node->R_L2P_operator * node->parent->L2L_operator[id_node_local];
  }



  if (node->nDOF < node->rank) {

    for (int lNgbr = 0; lNgbr < node->ngbr.m; lNgbr ++) {
      node->M2L_operator[node->Ilist.m + lNgbr] = node->ngbr(lNgbr)->L2P_operator * node->M2L_operator[node->Ilist.m + lNgbr] * node->P2M_operator; // slightly different from loop1
    }

    for (int lInteract = 0; lInteract < node->Ilist.m; lInteract ++) {
      node->M2L_operator[lInteract] = node->Ilist(lInteract)->L2P_operator * node->M2L_operator[lInteract] * node->P2M_operator; // slightly different from loop1
    }

    if (curr_level > 2) {
      // UPDATE PARENT INTERACTIONS
      int id_node_local = 0;
      for (int iChild = 0; iChild < 8; iChild ++) {
	if (node->parent->child[iChild] != NULL) {
	  if (node->parent->child[iChild]->id == node->id) {
	    id_node_local = iChild;
	  }
	}
      }
      node->parent->M2M_operator[id_node_local]=node->parent->M2M_operator[id_node_local]*node->P2M_operator; // same as original
      node->parent->L2L_operator[id_node_local]=node->L2P_operator*node->parent->L2L_operator[id_node_local]; // same as original
    }

    //20200108    node->P2M_operator = MatrixXcd::Identity(node->nDOF,node->nDOF);
    //20200108    node->L2P_operator = MatrixXcd::Identity(node->nDOF,node->nDOF);
    node->P2M_operator = MatrixXd::Identity(node->nDOF,node->nDOF);
    node->L2P_operator = MatrixXd::Identity(node->nDOF,node->nDOF);

    // UPDATE RANK 
    node->rank=node->nDOF;

  } else {

    /* this case is possible */

  }

}

#endif // PARA_CHECK_RANK
