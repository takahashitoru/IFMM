
		MatrixXd temp_K_ji_prime;
		jSV_M2P=0;
		// ACA_FullyPivoted(M2P_fillin_ji,U_j_prime,temp_K_ji_prime,epsilon,jSV_M2P,jSV_M2P);
		ACA_FullyPivoted(M2P_fillin_ji,U_j_prime,temp_K_ji_prime,epsilon_rel_fillin,jSV_M2P,0);
		K_ji_prime=temp_K_ji_prime.transpose();

		// if (jSV_M2P==21){
		//   cout << "M2P_fillin_ji" << endl << M2P_fillin_ji << endl;
		//   cout << "U_j_prime" << endl << U_j_prime << endl;
		//   cout << "K_ji_prime" << endl << K_ji_prime << endl;
		// }
		// printf("jSV_M2P: %d\n",jSV_M2P);




		jSV_P2L=0;
		// ACA_FullyPivoted(P2L_fillin_ij,K_ij_prime,V_j_prime,epsilon,jSV_P2L,jSV_P2L);
		ACA_FullyPivoted(P2L_fillin_ij,K_ij_prime,V_j_prime,epsilon_rel_fillin,jSV_P2L,0);

		NeedUpdate = (jSV_M2P==0 && jSV_P2L==0) ? false:true;




		// if (jSV_M2P==21){
		//   cout << "P2L_fillin_ij" << endl << P2L_fillin_ij << endl;
		//   cout << "K_ij_prime" << endl << K_ij_prime << endl;
		//   cout << "V_j_prime" << endl << V_j_prime << endl;
		// }

		// printf("M2P_fillin_ji.rows(),M2P_fillin_ji.cols(): %d %d\n",M2P_fillin_ji.rows(),M2P_fillin_ji.cols());
		// printf("P2L_fillin_ij.rows(),P2L_fillin_ij.cols(): %d %d\n",P2L_fillin_ij.rows(),P2L_fillin_ij.cols());
		// printf("jSV_M2P,jSV_P2L: %d %d\n",jSV_M2P,jSV_P2L);
                    
		// printf("on->ngbr(jNgbr)->rank: %d\n",on->ngbr(jNgbr)->rank);
		// printf("on->ngbr(jNgbr)->nDOF: %d\n",on->ngbr(jNgbr)->nDOF);


		outfile << "Rank M2P_ji (" << on->id << " " << on->ngbr(iNgbr)->id << " " << on->ngbr(jNgbr)->id << ")" << " " << jSV_M2P << endl;
		outfile << "Rank P2L_ij (" << on->id << " " << on->ngbr(iNgbr)->id << " " << on->ngbr(jNgbr)->id << ")" << " " << jSV_P2L << endl;



		if (NeedUpdate==true){




		  // RECOMPRESS THE INTERPOLATION OPERATOR Uj
		  U_comb_j = MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank+jSV_M2P);
		  U_comb_j.block(0,0,on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank) = on->ngbr(jNgbr)->L2P_operator;
		  U_comb_j.block(0,on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->nDOF,jSV_M2P) = U_j_prime;

		  // RECOMPRESS THE ANTERPOLATION OPERATOR Vj
		  V_comb_j = MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank+jSV_P2L);
		  V_comb_j.block(0,0,on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank) = on->ngbr(jNgbr)->P2M_operator.transpose();
		  V_comb_j.block(0,on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->nDOF,jSV_P2L) = V_j_prime;





		  MatrixXd temp_R_recomp_U_j;
		  SV_U_comb_j=0;
		  // SV_U_comb_j=on->ngbr(jNgbr)->rank; // FIX THE RANK !!!
		  // ACA_FullyPivoted(U_comb_j,U_recomp_j,temp_R_recomp_U_j,epsilon,SV_U_comb_j,SV_U_comb_j);
		  ACA_FullyPivoted(U_comb_j,U_recomp_j,temp_R_recomp_U_j,epsilon_rel_basis,SV_U_comb_j,0);
		  R_recomp_U_j=temp_R_recomp_U_j.transpose();

		  // printf("SV_U_comb_j: %d\n",SV_U_comb_j);

		  MatrixXd temp_R_recomp_V_j;
		  SV_V_comb_j=0;
		  // SV_V_comb_j=on->ngbr(jNgbr)->rank; // FIX THE RANK !!!
		  ACA_FullyPivoted(V_comb_j,V_recomp_j,temp_R_recomp_V_j,epsilon_rel_basis,SV_V_comb_j,SV_U_comb_j);
		  R_recomp_V_j=temp_R_recomp_V_j.transpose();

		  // printf("SV_V_comb_j: %d\n",SV_V_comb_j);



		  // cout << "V_comb_j" << endl << V_comb_j << endl;
		  // cout << "V_recomp_j" << endl << V_recomp_j << endl;

		  // cout << "on->ngbr(jNgbr)->P2M_operator.transpose()" << endl << on->ngbr(jNgbr)->P2M_operator.transpose() << endl;
		  // cout << "V_j_prime" << endl << V_j_prime << endl;



                    
		  // printf("pre-resize %d %d \n",SV_U_comb_j,SV_V_comb_j);
		  // MAKE SURE Uj AND Vj HAVE THE SAME RANK
		  if (SV_U_comb_j != SV_V_comb_j){
		    // printf("resize %d %d \n",SV_U_comb_j,SV_V_comb_j);
		    int temprank=min(min(on->ngbr(jNgbr)->rank+jSV_M2P,on->ngbr(jNgbr)->rank+jSV_P2L),max(SV_U_comb_j,SV_V_comb_j));
		    // cout << "temprank" << temprank << endl;

		    if (temprank==SV_U_comb_j){
		      SV_V_comb_j=temprank;

		      ACA_FullyPivoted(V_comb_j,V_recomp_j,temp_R_recomp_V_j,epsilon_rel_basis,SV_V_comb_j,SV_V_comb_j);
		      R_recomp_V_j=temp_R_recomp_V_j.transpose();
		    }
		    else if(temprank==SV_V_comb_j){
		      SV_U_comb_j=temprank;

		      ACA_FullyPivoted(U_comb_j,U_recomp_j,temp_R_recomp_U_j,epsilon_rel_basis,SV_U_comb_j,SV_U_comb_j);
		      R_recomp_U_j=temp_R_recomp_U_j.transpose();
		    }
		    else{
		      SV_U_comb_j=temprank;
		      SV_V_comb_j=temprank;
                        
		      ACA_FullyPivoted(U_comb_j,U_recomp_j,temp_R_recomp_U_j,epsilon_rel_basis,SV_U_comb_j,SV_U_comb_j);
		      R_recomp_U_j=temp_R_recomp_U_j.transpose();
		      ACA_FullyPivoted(V_comb_j,V_recomp_j,temp_R_recomp_V_j,epsilon_rel_basis,SV_V_comb_j,SV_V_comb_j);
		      R_recomp_V_j=temp_R_recomp_V_j.transpose();
		    }
		  }

		  // printf("check resize %d %d \n",SV_U_comb_j,SV_V_comb_j);
		  // printf("U_recomp_j.cols(): %d\n",U_recomp_j.cols());
		  // printf("V_recomp_j.cols(): %d\n",V_recomp_j.cols());

		}


