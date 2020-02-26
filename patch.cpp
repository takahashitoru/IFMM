//160331// PRECOMPUTE CERTAIN MATRICES HERE
//160331on->P2P_operator_self_LU = on->P2P_operator[iSelf].partialPivLu();
//160331MatrixXd rho_self=on->P2P_operator_self_LU.solve(on->L2P_operator);
//160331
//160331MatrixXd alpha_self=MatrixXd::Zero(on->rank,on->rank);
//160331alpha_self = - on->P2M_operator*rho_self;
//160331
//160331// LU DECOMPOSITION OF alpha_self
//160331on->alpha_self_LU = alpha_self.partialPivLu();
//160331MatrixXd alpha_self_inv=on->alpha_self_LU.solve(MatrixXd::Identity(on->rank,on->rank));

/*
  Once again, it seems that the pivot (V^T P2P^-1 U)^-1 of a few nodes
  becomes nearly singular, leading to an erroneous solution. Is the
  solution mainly incorrect for only a few elements, but correct for
  the remaining elements (as in the previous case, where only ~100 of
  the ~65000 elements were incorrect)?

  If that's indeed the case, can you try to add the code in attachment
  to ifmm.cpp? It should be included just after you check the norm of
  alpha_self_inv. If the norm is too large (e.g., >1e10), this piece
  of code will set U and V equal to the identity matrix, basically
  pushing all the variables to the parent level (resulting in
  alpha_self = P2P_self, which should solve the problem).

  I hope this will help, at least for these particular cases where a
  problem with one or two nodes in the tree leads to a large error in
  the solution. We might need to look further into this to understand
  what's actually going wrong with the method.
*/

#ifndef THRESHOLD_NEARLY_SINGULAR
//160427#define THRESHOLD_NEARLY_SINGULAR 1e+10
#define THRESHOLD_NEARLY_SINGULAR 1e+5
#endif

// CHECK IF alpha_self IS NEARLY SINGULAR OR NOT
//160331if (alpha_self_inv.norm() > 1e10){
if (alpha_self_inv.norm() > THRESHOLD_NEARLY_SINGULAR){

	// NEIGHBOUR LIST
	for (int lNgbr=0; lNgbr < on->ngbr.m; lNgbr++){
        on->M2L_operator[on->Ilist.m+lNgbr] = on->M2L_operator[on->Ilist.m+lNgbr]*on->P2M_operator;

        // also reverse
        int index_lNgbr_to_iSelf=0;
        for (int i=0; i < on->ngbr(lNgbr)->ngbr.m;i++){
          if (on->ngbr(lNgbr)->ngbr(i)->id == on->id)    
          {
            index_lNgbr_to_iSelf=i; // index l to i
            break;
          }
        }
        on->ngbr(lNgbr)->M2L_operator[on->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_iSelf] = on->L2P_operator*on->ngbr(lNgbr)->M2L_operator[on->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_iSelf] ;
    }

    // INTERACTION LIST
  	for (int lInteract=0; lInteract < on->Ilist.m; lInteract++){
        on->M2L_operator[lInteract] =  on->M2L_operator[lInteract]*on->P2M_operator;

        // also reverse
        int index_lNgbr_to_iSelf=0;
        for (int i=0; i < on->Ilist(lInteract)->Ilist.m;i++){
          if (on->Ilist(lInteract)->Ilist(i)->id == on->id)    
          {
            index_lNgbr_to_iSelf=i; // index l to i
            break;
          }
        }
        on->Ilist(lInteract)->M2L_operator[index_lNgbr_to_iSelf] = on->L2P_operator*on->Ilist(lInteract)->M2L_operator[index_lNgbr_to_iSelf];
     }

      
 	if(curr_level>2){
        // UPDATE PARENT INTERACTIONS
        int id_on_local=0;
        for (int iChild=0;iChild<8;iChild++){
         if (on->parent->child[iChild] != NULL){
           if (on->parent->child[iChild]->id == on->id){
             id_on_local=iChild;
            }
          }
        }
        on->parent->M2M_operator[id_on_local]=on->parent->M2M_operator[id_on_local]*on->P2M_operator;
        on->parent->L2L_operator[id_on_local]=on->L2P_operator*on->parent->L2L_operator[id_on_local];

  }

	//160331  on->P2M_operator=MatrixXd::Identity(on->nDOF,on->nDOF);
	//160331  on->L2P_operator=MatrixXd::Identity(on->nDOF,on->nDOF);
  on->P2M_operator=MatrixXcd::Identity(on->nDOF,on->nDOF);
  on->L2P_operator=MatrixXcd::Identity(on->nDOF,on->nDOF);

  // UPDATE RANK 
  on->rank=on->nDOF;


  // RECOMPUTE CERTAIN MATRICES HERE
	rho_self=on->P2P_operator_self_LU.solve(on->L2P_operator);
	//160331	alpha_self = MatrixXd::Zero(on->rank,on->rank);
	alpha_self = MatrixXcd::Zero(on->rank,on->rank);
	alpha_self = - on->P2M_operator*rho_self;

	// LU DECOMPOSITION OF alpha_self
	on->alpha_self_LU = alpha_self.partialPivLu();
	//160331	alpha_self_inv=on->alpha_self_LU.solve(MatrixXd::Identity(on->rank,on->rank));
	alpha_self_inv=on->alpha_self_LU.solve(MatrixXcd::Identity(on->rank,on->rank));

}

//160331VectorXd kappa_s = on->P2P_operator_self_LU.solve(on->RHS_leaf);
//160331VectorXd iota_s = on->alpha_self_LU.solve(on->P2M_operator*(kappa_s));
//160331
//160331// vector<MatrixXd> tau_si(on->ngbr.m);
//160331MatrixXd tau_si;
//160331vector<MatrixXd> beta_is(on->ngbr.m);
//160331// vector<MatrixXd> gamma_si(on->ngbr.m);
//160331MatrixXd gamma_si;
//160331// vector<MatrixXd> delta_ii(on->ngbr.m);
//160331vector<MatrixXd> xi_si(on->ngbr.m);
//160331vector<MatrixXd> psi_si(on->ngbr.m);
//160331
//160331// vector<MatrixXd> omicron_si(on->ngbr.m);
//160331MatrixXd omicron_si;
//160331// vector<MatrixXd> nu_si(on->ngbr.m);
//160331MatrixXd nu_si;
//160331vector<MatrixXd> mu_is(on->ngbr.m);
//160331// vector<MatrixXd> chi_ii(on->ngbr.m);
//160331vector<MatrixXd> zeta_si(on->ngbr.m);
//160331vector<MatrixXd> eta_si(on->ngbr.m);
//160331
//160331vector<MatrixXd> p2lto(on->ngbr.m);
//160331vector<MatrixXd> p2pto(on->ngbr.m);
