
		// SVD OF M2P_ji
#ifdef DISABLE_JACOBISVD
		LapackSVD<MatrixXd> svd_M2P_fillin_ji(M2P_fillin_ji, ComputeThinU | ComputeThinV);
#else
		JacobiSVD<MatrixXd> svd_M2P_fillin_ji(M2P_fillin_ji, ComputeThinU | ComputeThinV);
#endif
		// BDCSVD<MatrixXd> svd_M2P_fillin_ji(M2P_fillin_ji, ComputeThinU | ComputeThinV);
		VectorXd Sigma_j_M2P = svd_M2P_fillin_ji.singularValues();


		jSV_M2P=0;
		// while (Sigma_j_M2P(jSV_M2P) > epsilon*Sigma_j_M2P(0) && Sigma_j_M2P(jSV_M2P)>machine_eps && jSV_M2P < (min(M2P_fillin_ji.rows(),M2P_fillin_ji.cols())-1) ){
		while (Sigma_j_M2P(jSV_M2P) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j_M2P(jSV_M2P)>machine_eps && jSV_M2P < (min(M2P_fillin_ji.rows(),M2P_fillin_ji.cols())-1) ){
		  jSV_M2P++; 
		}
		// if (Sigma_j_M2P(jSV_M2P) > epsilon*Sigma_j_M2P(0) && Sigma_j_M2P(jSV_M2P)>machine_eps){
		if (Sigma_j_M2P(jSV_M2P) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j_M2P(jSV_M2P)>machine_eps){
		  jSV_M2P=min(M2P_fillin_ji.rows(),M2P_fillin_ji.cols());
		}
		// cout << "epsilon*Sigma_j_M2P(0)" << endl<< epsilon*Sigma_j_M2P(0) << endl;
		// cout << "Sigma_j_M2P(jSV_M2P-1)" << endl<< Sigma_j_M2P(jSV_M2P-1) << endl;
		// cout << "Sigma_j_M2P" << endl<< Sigma_j_M2P << endl;



		// SVD OF P2L_ij
#ifdef DISABLE_JACOBISVD
		LapackSVD<MatrixXd> svd_P2L_fillin_ij(P2L_fillin_ij, ComputeThinU | ComputeThinV);
#else
		JacobiSVD<MatrixXd> svd_P2L_fillin_ij(P2L_fillin_ij, ComputeThinU | ComputeThinV);
#endif
		// BDCSVD<MatrixXd> svd_P2L_fillin_ij(P2L_fillin_ij, ComputeThinU | ComputeThinV);
		VectorXd Sigma_j_P2L = svd_P2L_fillin_ij.singularValues();


		jSV_P2L=0;
		// while (Sigma_j_P2L(jSV_P2L) > epsilon*Sigma_j_P2L(0) && Sigma_j_P2L(jSV_P2L)>machine_eps && jSV_P2L < (min(P2L_fillin_ij.rows(),P2L_fillin_ij.cols())-1)){
		while (Sigma_j_P2L(jSV_P2L) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j_P2L(jSV_P2L)>machine_eps && jSV_P2L < (min(P2L_fillin_ij.rows(),P2L_fillin_ij.cols())-1)){
		  jSV_P2L++; 
		}
		// if (Sigma_j_P2L(jSV_P2L) > epsilon*Sigma_j_P2L(0) && Sigma_j_P2L(jSV_P2L)>machine_eps){
		if (Sigma_j_P2L(jSV_P2L) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j_P2L(jSV_P2L)>machine_eps){
		  jSV_P2L=min(P2L_fillin_ij.rows(),P2L_fillin_ij.cols());
		}

		if (jSV_P2L != jSV_M2P){
		  jSV_M2P=min(on->ngbr(jNgbr)->rank,max(jSV_M2P,jSV_P2L));
		  jSV_P2L = jSV_M2P;
		}

		// printf("M2P_fillin_ji.rows(), M2P_fillin_ji.cols(): %d %d \n",M2P_fillin_ji.rows(), M2P_fillin_ji.cols());
		// printf("P2L_fillin_ij.rows(), P2L_fillin_ij.cols(): %d %d \n",P2L_fillin_ij.rows(), P2L_fillin_ij.cols());

		outfile << "Rank M2P_ji (" << on->id << " " << on->ngbr(iNgbr)->id << " " << on->ngbr(jNgbr)->id << ")" << " " << jSV_M2P << endl;
		outfile << "Rank P2L_ij (" << on->id << " " << on->ngbr(iNgbr)->id << " " << on->ngbr(jNgbr)->id << ")" << " " << jSV_P2L << endl;

		NeedUpdate = (jSV_M2P==0) ? false:true;


		if (NeedUpdate==true){



		  Sigma_ji_prime=MatrixXd::Zero(jSV_M2P,jSV_M2P);
		  for (int index=0;index<jSV_M2P;index++){
		    Sigma_ji_prime(index,index)=Sigma_j_M2P(index);
		  }

		  U_j_prime = svd_M2P_fillin_ji.matrixU().block(0,0,on->ngbr(jNgbr)->nDOF,jSV_M2P);
		  K_ji_prime = Sigma_ji_prime*svd_M2P_fillin_ji.matrixV().block(0,0,on->ngbr(iNgbr)->rank,jSV_M2P).transpose();

		  // MatrixXd t1 = U_j_prime*K_ji_prime;
		  // printf("M2P_fillin_ji.norm(),t1.norm():%20.5e %20.5e\n",M2P_fillin_ji.norm(),t1.norm());

		  Sigma_ij_prime=MatrixXd::Zero(jSV_P2L,jSV_P2L);
		  for (int index=0;index<jSV_P2L;index++){
		    Sigma_ij_prime(index,index)=Sigma_j_P2L(index);
		  }

		  K_ij_prime = svd_P2L_fillin_ij.matrixU().block(0,0,on->ngbr(iNgbr)->rank,jSV_P2L)*Sigma_ij_prime;
		  V_j_prime = svd_P2L_fillin_ij.matrixV().block(0,0,on->ngbr(jNgbr)->nDOF,jSV_P2L);

		  // MatrixXd t2 = K_ij_prime*V_j_prime.transpose();
		  // printf("P2L_fillin_ij.norm(),t2.norm():%20.5e %20.5e\n",P2L_fillin_ij.norm(),t2.norm());


		  // OLD VERSION (WITHOUT SCALING) //
		  //-----------------------------------------------------------------------------------------------------------------------------//

		  // // printf("jSV_M2P,jSV_P2L: %d %d\n",jSV_M2P,jSV_P2L);
		  // // RECOMPRESS THE INTERPOLATION OPERATOR Uj
		  // U_comb_j = MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank+jSV_M2P);
		  // U_comb_j.block(0,0,on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank) = on->ngbr(jNgbr)->L2P_operator;
		  // U_comb_j.block(0,on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->nDOF,jSV_M2P) = U_j_prime;


		  // JacobiSVD<MatrixXd> svd_U_comb_j(U_comb_j, ComputeThinU | ComputeThinV);
		  // // BDCSVD<MatrixXd> svd_U_comb_j(U_comb_j, ComputeThinU | ComputeThinV);
		  // VectorXd Sigma_U_comb_j = svd_U_comb_j.singularValues();


		  // // cout << "on->ngbr(jNgbr)->L2P_operator" << endl << on->ngbr(jNgbr)->L2P_operator << endl;
		  // // cout << "U_j_prime" << endl << U_j_prime << endl;

		  // // // // printf("Rank M2P_ji (jSV_M2P): %d\n",jSV_M2P);
		  // // // // printf("Rank P2L_ij (jSV_P2L): %d\n",jSV_P2L);


                    

                    
		  // // int SV_U_comb_j=0;
		  // while (Sigma_U_comb_j(SV_U_comb_j) > epsilon_rel_basis*Sigma_U_comb_j(0) && Sigma_U_comb_j(SV_U_comb_j)>machine_eps && SV_U_comb_j < (min(U_comb_j.rows(),U_comb_j.cols())-1) ){
		  //  SV_U_comb_j++; 
		  // }
		  // if (Sigma_U_comb_j(SV_U_comb_j) > epsilon_rel_basis*Sigma_U_comb_j(0) && Sigma_U_comb_j(SV_U_comb_j)>machine_eps){
		  //   SV_U_comb_j=min(U_comb_j.rows(),U_comb_j.cols());
		  // }
                    

		  // // printf("jSV_M2P, jSV_P2L: %d %d \n",jSV_M2P, jSV_P2L);


                   
		  //  // RECOMPRESS THE ANTERPOLATION OPERATOR Vj
		  //  V_comb_j = MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank+jSV_P2L);
		  //  V_comb_j.block(0,0,on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank) = on->ngbr(jNgbr)->P2M_operator.transpose();
		  //  V_comb_j.block(0,on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->nDOF,jSV_P2L) = V_j_prime;


		  //  JacobiSVD<MatrixXd> svd_V_comb_j(V_comb_j, ComputeThinU | ComputeThinV);
		  //  // BDCSVD<MatrixXd> svd_V_comb_j(V_comb_j, ComputeThinU | ComputeThinV);
		  //  VectorXd Sigma_V_comb_j = svd_V_comb_j.singularValues();


		  //  // int SV_V_comb_j = 0;
		  //  // while (Sigma_V_comb_j(SV_V_comb_j) > epsilon*Sigma_V_comb_j(0) && Sigma_V_comb_j(SV_V_comb_j)>machine_eps && SV_V_comb_j < (min(V_comb_j.rows(),V_comb_j.cols())-1) ){
		  // // int SV_V_comb_j = 0;
		  //  while (Sigma_V_comb_j(SV_V_comb_j) > epsilon_rel_basis*Sigma_V_comb_j(0) && Sigma_V_comb_j(SV_V_comb_j)>machine_eps && SV_V_comb_j < (min(V_comb_j.rows(),V_comb_j.cols())-1) ){
		  //  // while (Sigma_V_comb_j(SV_V_comb_j) > epsilon_rel*sigma0_Aextended && Sigma_V_comb_j(SV_V_comb_j)>machine_eps && SV_V_comb_j < (min(V_comb_j.rows(),V_comb_j.cols())-1) ){
		  //   SV_V_comb_j++; 
		  //  }
		  //  if (Sigma_V_comb_j(SV_V_comb_j) > epsilon_rel_basis*Sigma_V_comb_j(0) && Sigma_V_comb_j(SV_V_comb_j)>machine_eps){
		  //  // if (Sigma_V_comb_j(SV_V_comb_j) > epsilon_rel*sigma0_Aextended && Sigma_V_comb_j(SV_V_comb_j)>machine_eps){
		  //    SV_V_comb_j=min(V_comb_j.rows(),V_comb_j.cols());
		  //  }


                   
                    
		  //  // MAKE SURE Uj AND Vj HAVE THE SAME RANK
		  //  if (SV_U_comb_j != SV_V_comb_j){
		  //    SV_U_comb_j=min(min(on->ngbr(jNgbr)->rank+jSV_M2P,on->ngbr(jNgbr)->rank+jSV_P2L),max(SV_U_comb_j,SV_V_comb_j));
		  //    // SV_U_comb_j=min(SV_U_comb_j,SV_V_comb_j);
		  //    SV_V_comb_j = SV_U_comb_j;
		  //  }

		  //  // cout << "Sigma_U_comb_j" << endl << Sigma_U_comb_j << endl;
		  //  // cout << "svd_U_comb_j.matrixU()" << endl << svd_U_comb_j.matrixU() << endl;

		  //  MatrixXd Sigma_diag_U_comb_j=MatrixXd::Zero(SV_U_comb_j,SV_U_comb_j);
		  //  for (int index=0;index<SV_U_comb_j;index++){
		  //    Sigma_diag_U_comb_j(index,index)=Sigma_U_comb_j(index);
		  //  }

		  //  U_recomp_j = svd_U_comb_j.matrixU().block(0,0,on->ngbr(jNgbr)->nDOF,SV_U_comb_j);
		  //  R_recomp_U_j = Sigma_diag_U_comb_j*(svd_U_comb_j.matrixV().block(0,0,on->ngbr(jNgbr)->rank+jSV_M2P,SV_U_comb_j).transpose());

		  //  MatrixXd Sigma_diag_V_comb_j=MatrixXd::Zero(SV_V_comb_j,SV_V_comb_j);
		  //  for (int index=0;index<SV_V_comb_j;index++){
		  //    Sigma_diag_V_comb_j(index,index)=Sigma_V_comb_j(index);
		  //  }
		  //  V_recomp_j = svd_V_comb_j.matrixU().block(0,0,on->ngbr(jNgbr)->nDOF,SV_V_comb_j);
		  //  R_recomp_V_j = Sigma_diag_V_comb_j*(svd_V_comb_j.matrixV().block(0,0,on->ngbr(jNgbr)->rank+jSV_P2L,SV_V_comb_j).transpose());

		  //-----------------------------------------------------------------------------------------------------------------------------//

		  // NEW VERSION (WITH SCALING)

		  // RECOMPRESS THE SCALED (!) INTERPOLATION OPERATOR Uj
		  U_comb_j = MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank+jSV_M2P);
		  U_comb_j.block(0,0,on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank) = on->ngbr(jNgbr)->L2P_operator*on->ngbr(jNgbr)->Sigma_L2P;
		  U_comb_j.block(0,on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->nDOF,jSV_M2P) = U_j_prime * Sigma_ji_prime;

#ifdef DISABLE_JACOBISVD
		  LapackSVD<MatrixXd> svd_U_comb_j(U_comb_j, ComputeThinU | ComputeThinV);
#else
		  JacobiSVD<MatrixXd> svd_U_comb_j(U_comb_j, ComputeThinU | ComputeThinV);
#endif
		  VectorXd Sigma_U_comb_j = svd_U_comb_j.singularValues();

		  SV_U_comb_j=0;
		  while (Sigma_U_comb_j(SV_U_comb_j) > epsilon_rel_basis*Sigma_U_comb_j(0) && Sigma_U_comb_j(SV_U_comb_j)>machine_eps && SV_U_comb_j < (min(U_comb_j.rows(),U_comb_j.cols())-1) ){
		    SV_U_comb_j++; 
		  }
		  if (Sigma_U_comb_j(SV_U_comb_j) > epsilon_rel_basis*Sigma_U_comb_j(0) && Sigma_U_comb_j(SV_U_comb_j)>machine_eps){
		    SV_U_comb_j=min(U_comb_j.rows(),U_comb_j.cols());
		  }

		  // RECOMPRESS THE SCALED (!) ANTERPOLATION OPERATOR Vj
		  V_comb_j = MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank+jSV_P2L);
		  V_comb_j.block(0,0,on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank) = on->ngbr(jNgbr)->P2M_operator.transpose()*on->ngbr(jNgbr)->Sigma_P2M;
		  V_comb_j.block(0,on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->nDOF,jSV_P2L) = V_j_prime*Sigma_ij_prime;

#ifdef DISABLE_JACOBISVD
		  LapackSVD<MatrixXd> svd_V_comb_j(V_comb_j, ComputeThinU | ComputeThinV);
#else
		  JacobiSVD<MatrixXd> svd_V_comb_j(V_comb_j, ComputeThinU | ComputeThinV);
#endif
		  VectorXd Sigma_V_comb_j = svd_V_comb_j.singularValues();

		  SV_V_comb_j = 0;
		  while (Sigma_V_comb_j(SV_V_comb_j) > epsilon_rel_basis*Sigma_V_comb_j(0) && Sigma_V_comb_j(SV_V_comb_j)>machine_eps && SV_V_comb_j < (min(V_comb_j.rows(),V_comb_j.cols())-1) ){
		    SV_V_comb_j++; 
		  }
		  if (Sigma_V_comb_j(SV_V_comb_j) > epsilon_rel_basis*Sigma_V_comb_j(0) && Sigma_V_comb_j(SV_V_comb_j)>machine_eps){
		    SV_V_comb_j=min(V_comb_j.rows(),V_comb_j.cols());
		  }
                   
                    
		  // MAKE SURE Uj AND Vj HAVE THE SAME RANK
		  if (SV_U_comb_j != SV_V_comb_j){
		    SV_U_comb_j=min(min(on->ngbr(jNgbr)->rank+jSV_M2P,on->ngbr(jNgbr)->rank+jSV_P2L),max(SV_U_comb_j,SV_V_comb_j));
		    // SV_U_comb_j=min(SV_U_comb_j,SV_V_comb_j);
		    SV_V_comb_j = SV_U_comb_j;
		  }


		  // RECOMPRESS Uj
		  Sigma_diag_U_comb_j=MatrixXd::Zero(SV_U_comb_j,SV_U_comb_j);
		  for (int index=0;index<SV_U_comb_j;index++){
		    Sigma_diag_U_comb_j(index,index)=Sigma_U_comb_j(index);
		  }

		  U_recomp_j = svd_U_comb_j.matrixU().block(0,0,on->ngbr(jNgbr)->nDOF,SV_U_comb_j);
		  R_recomp_U_j = Sigma_diag_U_comb_j*(svd_U_comb_j.matrixV().block(0,0,on->ngbr(jNgbr)->rank+jSV_M2P,SV_U_comb_j).transpose());

		  // RECOMPRESS Vj
		  Sigma_diag_V_comb_j=MatrixXd::Zero(SV_V_comb_j,SV_V_comb_j);
		  for (int index=0;index<SV_V_comb_j;index++){
		    Sigma_diag_V_comb_j(index,index)=Sigma_V_comb_j(index);
		  }
		  V_recomp_j = svd_V_comb_j.matrixU().block(0,0,on->ngbr(jNgbr)->nDOF,SV_V_comb_j);
		  R_recomp_V_j = Sigma_diag_V_comb_j*(svd_V_comb_j.matrixV().block(0,0,on->ngbr(jNgbr)->rank+jSV_P2L,SV_V_comb_j).transpose());



                    

		}

