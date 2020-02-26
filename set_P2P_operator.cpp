#include "set_P2P_operator.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

// void set_P2P_operator(oct_node * node,Stokes_mobility * mfun, bool &IsLeaf, double &DX){
void set_P2P_operator(oct_node * node,kernel * kfun, bool &IsLeaf){

  // printf("node->id: %d\n",node->id);

  if (IsLeaf){ // leaf level


#if defined(_OPENMP)

    // In this case, the loop including set_P2P_operator() is assumed
    // to be parallelized. Then, the symmetry of P2P operators in the
    // original code must not be used. Toru

    for (int iNgbr = 0; iNgbr < node->ngbr.m; ++iNgbr) {

      node->P2P_operator.push_back(MatrixXd::Zero(node->ngbr(iNgbr)->pnt.m,node->pnt.m));

      for(int i=0;i<node->ngbr(iNgbr)->pnt.m;++i){
	for(int j=0;j<node->pnt.m;++j){
	  vec3 r = node->pnt.p[j] - node->ngbr(iNgbr)->pnt.p[i];
	  node->P2P_operator[iNgbr](i,j) = (*kfun)(r);
	}
      }

    } // iNgbr


#else // original

    int iSelf=0; // index to current node (in the neighbour list of the current node)
    for (int iNgbr=0; iNgbr < node-> ngbr.m; iNgbr++){
      if (node->ngbr(iNgbr)->id==node->id){
        iSelf=iNgbr;
        break;
      }
    }


    for (int iNgbr = 0; iNgbr < node->ngbr.m; ++iNgbr) {
       node->P2P_operator.push_back(MatrixXd::Zero(node->ngbr(iNgbr)->pnt.m,node->pnt.m));

       // MatrixXd test = MatrixXd:: Zero(node->ngbr(iNgbr)->pnt.m,node->pnt.m);

       if (iNgbr==iSelf || node->ngbr(iNgbr)->P2P_operator.size()==0){

         for(int i=0;i<node->ngbr(iNgbr)->pnt.m;++i){
            for(int j=0;j<node->pnt.m;++j){
              vec3 r = node->pnt.p[j] - node->ngbr(iNgbr)->pnt.p[i];
              // for (int iDim=0; iDim<NDIM; iDim++){
                // for (int jDim=0; jDim<NDIM; jDim++){
                  // node->P2P_operator[iNgbr](i*NDIM+iDim,j*NDIM+jDim) = (*mfun)(r, DX, i*NDIM+iDim,j*NDIM+jDim);
                  node->P2P_operator[iNgbr](i,j) = (*kfun)(r);
                // }
              // }
            }
          }

        }
        else{
          int index_iNgbr_to_iSelf=0;
          for (int i=0; i < node->ngbr(iNgbr)->ngbr.m;i++){
            if (node->ngbr(iNgbr)->ngbr(i)->id == node->id)
            {
              index_iNgbr_to_iSelf=i; // index to current node (in the neighbour list of i)
              break;
            }
          }
          
          node->P2P_operator[iNgbr] = node->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf].transpose();
        
        }


        // printf("node->ngbr(iNgbr)->pnt.m,node->pnt.m,: %d %d\n",node->ngbr(iNgbr)->pnt.m,node->pnt.m);
        // printf("node->P2P_operator[iNgbr].rows(),node->P2P_operator[iNgbr].cols(): %d %d\n",node->P2P_operator[iNgbr].rows(),node->P2P_operator[iNgbr].cols());
        // cout << "node->P2P_operator[iNgbr]" << endl << node->P2P_operator[iNgbr] << endl;
    }

#endif

  }
  else{
    // P2P AT NON-LEAF LEVEL IS CONSTRUCTED BASED ON M2L-OPERATORS OF THE CHILD LEVEL
    // SEE TRANSFER_M2L_TO_P2P
  }

}
