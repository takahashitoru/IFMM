
		// cout << "M2P_rsvd" << endl;

		int rank_RSVD = max(rank_RSVD_min,(int) min(M2P_fillin_ji.rows(),M2P_fillin_ji.cols())/10);
		if (rank_RSVD > min(M2P_fillin_ji.rows(),M2P_fillin_ji.cols())){
		  rank_RSVD = min(M2P_fillin_ji.rows(),M2P_fillin_ji.cols());
		}

		// SVD OF M2P_ji
		// printf("rank_RSVD: %d\n",rank_RSVD);
		RedSVD<MatrixXd> rsvd_M2P_fillin_ji(M2P_fillin_ji,rank_RSVD);
		VectorXd Sigma_j_M2P = rsvd_M2P_fillin_ji.singularValues();


		jSV_M2P=0;
		while (Sigma_j_M2P(jSV_M2P) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j_M2P(jSV_M2P)>machine_eps && jSV_M2P < (rank_RSVD-1) ){
		  jSV_M2P++; 
		}
		if (Sigma_j_M2P(jSV_M2P) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j_M2P(jSV_M2P)>machine_eps){
		  // rank_RSVD=min(M2P_fillin_ji.rows(),M2P_fillin_ji.cols());

		  // SVD OF M2P_ji
		  rsvd_M2P_fillin_ji.compute(M2P_fillin_ji,min(M2P_fillin_ji.rows(),M2P_fillin_ji.cols()));
		  Sigma_j_M2P = rsvd_M2P_fillin_ji.singularValues();



		  jSV_M2P=0;
		  while (Sigma_j_M2P(jSV_M2P) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j_M2P(jSV_M2P)>machine_eps && jSV_M2P < (min(M2P_fillin_ji.rows(),M2P_fillin_ji.cols())-1)){
		    jSV_M2P++; 
		  }
		  if (Sigma_j_M2P(jSV_M2P) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j_M2P(jSV_M2P)>machine_eps){
		    jSV_M2P=min(M2P_fillin_ji.rows(),M2P_fillin_ji.cols());
		  }
		}

		// printf("jSV_M2P: %d\n",jSV_M2P);




		// SVD OF P2L_ij
		// printf("rank_RSVD: %d\n",rank_RSVD);
		RedSVD<MatrixXd> rsvd_P2L_fillin_ij(P2L_fillin_ij,rank_RSVD);
		VectorXd Sigma_j_P2L = rsvd_P2L_fillin_ij.singularValues();


		jSV_P2L=0;
		while (Sigma_j_P2L(jSV_P2L) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j_P2L(jSV_P2L)>machine_eps && jSV_P2L < (rank_RSVD-1) ){
		  jSV_P2L++; 
		}
		// printf("jSV_P2L: %d\n",jSV_P2L);

		if (Sigma_j_P2L(jSV_P2L) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j_P2L(jSV_P2L)>machine_eps){
		  // rank_RSVD=min(M2P_fillin_ji.rows(),M2P_fillin_ji.cols());

		  // SVD OF P2L_ij
		  rsvd_P2L_fillin_ij.compute(P2L_fillin_ij,min(P2L_fillin_ij.rows(),P2L_fillin_ij.cols()));
		  Sigma_j_P2L = rsvd_P2L_fillin_ij.singularValues();


		  jSV_P2L=0;
		  while (Sigma_j_P2L(jSV_P2L) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j_P2L(jSV_P2L)>machine_eps && jSV_P2L < (min(P2L_fillin_ij.rows(),P2L_fillin_ij.cols())-1)){
		    jSV_P2L++; 
		  }
		  if (Sigma_j_P2L(jSV_P2L) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j_P2L(jSV_P2L)>machine_eps){
		    jSV_P2L=min(P2L_fillin_ij.rows(),P2L_fillin_ij.cols());
		  }
		}

		// printf("jSV_M2P,jSV_P2L: %d %d\n",jSV_M2P,jSV_P2L);


		if (jSV_P2L != jSV_M2P){
		  int jSV_M2P_temp=jSV_M2P;
		  int jSV_P2L_temp=jSV_P2L;

		  jSV_M2P=min(on->ngbr(jNgbr)->rank,max(jSV_M2P,jSV_P2L));
		  jSV_P2L = jSV_M2P;

                      
		  if (jSV_M2P>rank_RSVD && jSV_M2P > jSV_M2P_temp){
		    // RedSVD<MatrixXd> rsvd_M2P_fillin_ji_new(M2P_fillin_ji,jSV_M2P);
		    rsvd_M2P_fillin_ji.compute(M2P_fillin_ji,jSV_M2P);
		    Sigma_j_M2P = rsvd_M2P_fillin_ji.singularValues();
		  }
		  if (jSV_P2L>rank_RSVD && jSV_P2L > jSV_P2L_temp){
		    // RedSVD<MatrixXd> rsvd_P2L_fillin_ij_new(P2L_fillin_ij,jSV_P2L);
		    rsvd_P2L_fillin_ij.compute(P2L_fillin_ij,jSV_P2L);
		    Sigma_j_P2L = rsvd_P2L_fillin_ij.singularValues();
		  }

                        
		  // printf("jSV_M2P,jSV_P2L: %d %d\n",jSV_M2P,jSV_P2L);

                      
		  // // RedSVD<MatrixXd> rsvd_M2P_fillin_ji(M2P_fillin_ji,jSV_M2P);
		  // rsvd_M2P_fillin_ji.compute(M2P_fillin_ji,max(jSV_M2P,rank_RSVD));
		  // Sigma_j_M2P = rsvd_M2P_fillin_ji.singularValues();

		  // rsvd_P2L_fillin_ij.compute(P2L_fillin_ij,max(jSV_P2L,rank_RSVD));
		  // Sigma_j_P2L = rsvd_P2L_fillin_ij.singularValues();

		  // while ((Sigma_j_P2L(jSV_M2P-1) < 1000*machine_eps || Sigma_j_M2P(jSV_P2L-1)< 1000*machine_eps) && jSV_M2P>0){
		  //   cout << "t0" << endl;
		  //   jSV_M2P--;
		  //   jSV_P2L--;
		  // }

		  while (jSV_M2P>0){
		    if (Sigma_j_P2L(jSV_M2P-1) < 1000*machine_eps || Sigma_j_M2P(jSV_P2L-1)< 1000*machine_eps){
		      jSV_M2P--;
		      jSV_P2L--;
		    }
		    else{
		      break;
		    }
		  }
                      
		  // if (Sigma_j_P2L(jSV_M2P-1) < 1000*machine_eps || Sigma_j_M2P(jSV_P2L-1)< 1000*machine_eps){
		  //   jSV_P2L=min(jSV_M2P_temp,jSV_P2L_temp);
		  //   jSV_M2P=min(jSV_M2P_temp,jSV_P2L_temp);
		  // }

		}



		// ANALYSE outfile << "Rank M2P_ji (" << on->id << " " << on->ngbr(iNgbr)->id << " " << on->ngbr(jNgbr)->id << ")" << " " << jSV_M2P << endl;
		// ANALYSE outfile << "Rank P2L_ij (" << on->id << " " << on->ngbr(iNgbr)->id << " " << on->ngbr(jNgbr)->id << ")" << " " << jSV_P2L << endl;


		NeedUpdate = (jSV_M2P==0) ? false:true;

		if (NeedUpdate==true){     

		  // cout << "NeedUpdate" << endl;                 

		  MatrixXd Sigma_ji_prime=MatrixXd::Zero(jSV_M2P,jSV_M2P);
		  for (int index=0;index<jSV_M2P;index++){
		    Sigma_ji_prime(index,index)=Sigma_j_M2P(index);
		  }

                      

		  U_j_prime = rsvd_M2P_fillin_ji.matrixU().block(0,0,on->ngbr(jNgbr)->nDOF,jSV_M2P);
		  K_ji_prime = Sigma_ji_prime*rsvd_M2P_fillin_ji.matrixV().block(0,0,on->ngbr(iNgbr)->rank,jSV_M2P).transpose();
                      

		  MatrixXd Sigma_ij_prime=MatrixXd::Zero(jSV_P2L,jSV_P2L);
		  for (int index=0;index<jSV_P2L;index++){
		    Sigma_ij_prime(index,index)=Sigma_j_P2L(index);
		  }


		  K_ij_prime = rsvd_P2L_fillin_ij.matrixU().block(0,0,on->ngbr(iNgbr)->rank,jSV_P2L)*Sigma_ij_prime;
		  V_j_prime = rsvd_P2L_fillin_ij.matrixV().block(0,0,on->ngbr(jNgbr)->nDOF,jSV_P2L);
                      

		  // // RECOMPRESS THE INTERPOLATION OPERATOR Uj
		  // U_comb_j = MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank+jSV_M2P);
		  // U_comb_j.block(0,0,on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank) = on->ngbr(jNgbr)->L2P_operator;
		  // U_comb_j.block(0,on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->nDOF,jSV_M2P) = U_j_prime;


		  // JacobiSVD<MatrixXd> svd_U_comb_j(U_comb_j, ComputeThinU | ComputeThinV);
		  // VectorXd Sigma_U_comb_j = svd_U_comb_j.singularValues();

		  // while (Sigma_U_comb_j(SV_U_comb_j) > epsilon_rel*Sigma_U_comb_j(0) && Sigma_U_comb_j(SV_U_comb_j)>machine_eps && SV_U_comb_j < (min(U_comb_j.rows(),U_comb_j.cols())-1) ){
		  // // while (Sigma_U_comb_j(SV_U_comb_j) > epsilon && Sigma_U_comb_j(SV_U_comb_j)>machine_eps && SV_U_comb_j < (min(U_comb_j.rows(),U_comb_j.cols())-1) ){
		  //  SV_U_comb_j++; 
		  // }
		  // if (Sigma_U_comb_j(SV_U_comb_j) > epsilon_rel*Sigma_U_comb_j(0)  && Sigma_U_comb_j(SV_U_comb_j)>machine_eps){
		  // // if (Sigma_U_comb_j(SV_U_comb_j) > epsilon && Sigma_U_comb_j(SV_U_comb_j)>machine_eps){
		  //   SV_U_comb_j=min(U_comb_j.rows(),U_comb_j.cols());
		  // }

		  // // RECOMPRESS THE ANTERPOLATION OPERATOR Vj
		  // V_comb_j = MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank+jSV_P2L);
		  // V_comb_j.block(0,0,on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank) = on->ngbr(jNgbr)->P2M_operator.transpose();
		  // V_comb_j.block(0,on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->nDOF,jSV_P2L) = V_j_prime;



		  // JacobiSVD<MatrixXd> svd_V_comb_j(V_comb_j, ComputeThinU | ComputeThinV);
		  // VectorXd Sigma_V_comb_j = svd_V_comb_j.singularValues();

		  // while (Sigma_V_comb_j(SV_V_comb_j) > epsilon_rel*Sigma_V_comb_j(0) && Sigma_V_comb_j(SV_V_comb_j)>machine_eps && SV_V_comb_j < (min(V_comb_j.rows(),V_comb_j.cols())-1) ){
		  // // while (Sigma_V_comb_j(SV_V_comb_j) > epsilon && Sigma_V_comb_j(SV_V_comb_j)>machine_eps && SV_V_comb_j < (min(V_comb_j.rows(),V_comb_j.cols())-1) ){
		  //  SV_V_comb_j++; 
		  // }
		  // if (Sigma_V_comb_j(SV_V_comb_j) > epsilon_rel*Sigma_V_comb_j(0) && Sigma_V_comb_j(SV_V_comb_j)>machine_eps){
		  // // if (Sigma_V_comb_j(SV_V_comb_j) > epsilon && Sigma_V_comb_j(SV_V_comb_j)>machine_eps){
		  //   SV_V_comb_j=min(V_comb_j.rows(),V_comb_j.cols());
		  // }


                      
		  // // MAKE SURE Uj AND Vj HAVE THE SAME RANK
		  // if (SV_U_comb_j != SV_V_comb_j){
		  //   SV_U_comb_j=min(min(on->ngbr(jNgbr)->rank+jSV_M2P,on->ngbr(jNgbr)->rank+jSV_P2L),max(SV_U_comb_j,SV_V_comb_j));
		  //   SV_V_comb_j = SV_U_comb_j;
		  // }


		  // MatrixXd Sigma_diag_U_comb_j=MatrixXd::Zero(SV_U_comb_j,SV_U_comb_j);
		  // for (int index=0;index<SV_U_comb_j;index++){
		  //   Sigma_diag_U_comb_j(index,index)=Sigma_U_comb_j(index);
		  // }

		  // U_recomp_j = svd_U_comb_j.matrixU().block(0,0,on->ngbr(jNgbr)->nDOF,SV_U_comb_j);
		  // R_recomp_U_j = Sigma_diag_U_comb_j*(svd_U_comb_j.matrixV().block(0,0,on->ngbr(jNgbr)->rank+jSV_M2P,SV_U_comb_j).transpose());

		  // MatrixXd Sigma_diag_V_comb_j=MatrixXd::Zero(SV_V_comb_j,SV_V_comb_j);
		  // for (int index=0;index<SV_V_comb_j;index++){
		  //   Sigma_diag_V_comb_j(index,index)=Sigma_V_comb_j(index);
		  // }
		  // V_recomp_j = svd_V_comb_j.matrixU().block(0,0,on->ngbr(jNgbr)->nDOF,SV_V_comb_j);
		  // R_recomp_V_j = Sigma_diag_V_comb_j*(svd_V_comb_j.matrixV().block(0,0,on->ngbr(jNgbr)->rank+jSV_P2L,SV_V_comb_j).transpose());

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


		  // MAKE SURE Uj AND Vj HAVE THE SAME RANK
		  if (SV_U_comb_j != SV_V_comb_j){
		    // printf("resize %d %d \n",SV_U_comb_j,SV_V_comb_j);
		    int temprank=min(min(on->ngbr(jNgbr)->rank+jSV_M2P,on->ngbr(jNgbr)->rank+jSV_P2L),max(SV_U_comb_j,SV_V_comb_j));
		    // int temprank=min(min(on->ngbr(jNgbr)->rank+jSV_M2P,on->ngbr(jNgbr)->rank+jSV_P2L),min(SV_U_comb_j,SV_V_comb_j));
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

                      

		}

