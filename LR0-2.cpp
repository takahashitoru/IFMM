		  // cout << "P2P_fillin_ij.rows(): " << P2P_fillin_ij.rows() << endl;
		  // cout << "P2P_fillin_ij.cols(): " << P2P_fillin_ij.cols() << endl;

		  // SVD OF P2P_ij
#ifdef DISABLE_JACOBISVD
		  LapackSVD<MatrixXd> svd_P2P_fillin_ij(P2P_fillin_ij, ComputeThinU | ComputeThinV);
#else
		  JacobiSVD<MatrixXd> svd_P2P_fillin_ij(P2P_fillin_ij, ComputeThinU | ComputeThinV);
#endif
		  // BDCSVD<MatrixXd> svd_P2P_fillin_ij(P2P_fillin_ij, ComputeThinU | ComputeThinV);
		  VectorXd Sigma_i = svd_P2P_fillin_ij.singularValues();
                    
                    
		  iSV=0;
		  // while (Sigma_i(iSV) > epsilon*Sigma_i(0) && Sigma_i(iSV) > machine_eps && iSV< (min(P2P_fillin_ij.rows(),P2P_fillin_ij.cols())-1)){
		  while (Sigma_i(iSV) > epsilon_rel_fillin*sigma0_Aextended && Sigma_i(iSV) > machine_eps && iSV< (min(P2P_fillin_ij.rows(),P2P_fillin_ij.cols())-1)){
		    iSV++; 
		  }
		  // if (Sigma_i(iSV) > epsilon*Sigma_i(0) && Sigma_i(iSV) > machine_eps){
		  if (Sigma_i(iSV) > epsilon_rel_fillin*sigma0_Aextended && Sigma_i(iSV) > machine_eps){
		    iSV=min(P2P_fillin_ij.rows(),P2P_fillin_ij.cols());
		  }

		  // SVD OF P2P_ji
#ifdef DISABLE_JACOBISVD
		  LapackSVD<MatrixXd> svd_P2P_fillin_ji(P2P_fillin_ji, ComputeThinU | ComputeThinV);
#else
		  JacobiSVD<MatrixXd> svd_P2P_fillin_ji(P2P_fillin_ji, ComputeThinU | ComputeThinV);
#endif
		  // BDCSVD<MatrixXd> svd_P2P_fillin_ji(P2P_fillin_ji, ComputeThinU | ComputeThinV);
		  VectorXd Sigma_j = svd_P2P_fillin_ji.singularValues();

                    
		  jSV=0;
		  // while (Sigma_j(jSV) > epsilon*Sigma_j(0) && Sigma_j(jSV) > machine_eps && jSV < (min(P2P_fillin_ji.rows(),P2P_fillin_ji.cols())-1)){
		  while (Sigma_j(jSV) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j(jSV) > machine_eps && jSV < (min(P2P_fillin_ji.rows(),P2P_fillin_ji.cols())-1)){
		    jSV++; 
		  }
		  // if (Sigma_j(jSV) > epsilon*Sigma_j(0) && Sigma_j(jSV) > machine_eps){
		  if (Sigma_j(jSV) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j(jSV) > machine_eps){
		    jSV=min(P2P_fillin_ji.rows(),P2P_fillin_ji.cols());
		  }

		  //   printf("iSV,jSV,%d %d \n",iSV,jSV);
		  //   cout << "Sigma_i" << endl << Sigma_i << endl;
		  //   cout << "Sigma_j" << endl << Sigma_j << endl;

		  if (iSV != jSV){
		    iSV=min(min(on->ngbr(iNgbr)->nDOF,on->ngbr(jNgbr)->nDOF),max(jSV,jSV));
		    jSV = iSV;
		  }


		  outfile << "Rank P2P_ij (" << on->id << " " << on->ngbr(iNgbr)->id << " " << on->ngbr(jNgbr)->id << ")" << " " << iSV << endl;
		  outfile << "Rank P2P_ji (" << on->id << " " << on->ngbr(iNgbr)->id << " " << on->ngbr(jNgbr)->id << ")" << " " << jSV << endl;


		  NeedUpdate = (iSV==0) ? false:true;


		  if (NeedUpdate==true){




                    K_ij_prime=MatrixXd::Zero(iSV,iSV);
                    for (int index=0;index<iSV;index++){
                      K_ij_prime(index,index)=Sigma_i(index);
                    }

                    U_i_prime = svd_P2P_fillin_ij.matrixU().block(0,0,on->ngbr(iNgbr)->nDOF,iSV);
                    V_j_prime = svd_P2P_fillin_ij.matrixV().block(0,0,on->ngbr(jNgbr)->nDOF,iSV);


                    // OLD VERSION (WITHOUT SCALING) //
                    //----------------------------------------------------------------------------------------------------------------------------//
                    // // RECOMPRESS THE INTERPOLATION OPERATOR Ui
                    // U_comb_i = MatrixXd::Zero(on->ngbr(iNgbr)->nDOF,on->ngbr(iNgbr)->rank+iSV);
                    // U_comb_i.block(0,0,on->ngbr(iNgbr)->nDOF,on->ngbr(iNgbr)->rank) = on->ngbr(iNgbr)->L2P_operator;
                    // U_comb_i.block(0,on->ngbr(iNgbr)->rank,on->ngbr(iNgbr)->nDOF,iSV) = U_i_prime;

                    
                    // JacobiSVD<MatrixXd> svd_U_comb_i(U_comb_i, ComputeThinU | ComputeThinV);
                    // // BDCSVD<MatrixXd> svd_U_comb_i(U_comb_i, ComputeThinU | ComputeThinV);
                    // VectorXd Sigma_U_comb_i = svd_U_comb_i.singularValues();


                    // SV_U_comb_i=0;
                    // while (Sigma_U_comb_i(SV_U_comb_i) > epsilon_rel_basis*Sigma_U_comb_i(0) && Sigma_U_comb_i(SV_U_comb_i) > machine_eps && SV_U_comb_i<(Sigma_U_comb_i.rows()-1) ){
                    // // while (Sigma_U_comb_i(SV_U_comb_i) > epsilon_rel_fillin*sigma0_Aextended && Sigma_U_comb_i(SV_U_comb_i) > machine_eps && SV_U_comb_i<(Sigma_U_comb_i.rows()-1) ){
                    //  SV_U_comb_i++; 
                    // }
                    // if (Sigma_U_comb_i(SV_U_comb_i) > epsilon_rel_basis*Sigma_U_comb_i(0) && Sigma_U_comb_i(SV_U_comb_i) > machine_eps){
                    // // if (Sigma_U_comb_i(SV_U_comb_i) > epsilon_rel*sigma0_Aextended && Sigma_U_comb_i(SV_U_comb_i) > machine_eps){
                    //   SV_U_comb_i=Sigma_U_comb_i.rows();
                    // }

                    // // RECOMPRESS THE ANTERPOLATION OPERATOR Vj
                    // V_comb_j =  MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank+iSV);
                    // V_comb_j.block(0,0,on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank) = on->ngbr(jNgbr)->P2M_operator.transpose();
                    // V_comb_j.block(0,on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->nDOF,iSV) = V_j_prime;

                    // JacobiSVD<MatrixXd> svd_V_comb_j(V_comb_j, ComputeThinU | ComputeThinV);
                    // // BDCSVD<MatrixXd> svd_V_comb_j(V_comb_j, ComputeThinU | ComputeThinV);
                    // VectorXd Sigma_V_comb_j = svd_V_comb_j.singularValues();


                    // // int SV_V_comb_j=0;
                    // // while (Sigma_V_comb_j(SV_V_comb_j) > epsilon*Sigma_V_comb_j(0) && Sigma_V_comb_j(SV_V_comb_j) > machine_eps && SV_V_comb_j<(Sigma_V_comb_j.rows()-1) ){
                    // while (Sigma_V_comb_j(SV_V_comb_j) > epsilon_rel_basis*Sigma_V_comb_j(0) && Sigma_V_comb_j(SV_V_comb_j) > machine_eps && SV_V_comb_j<(Sigma_V_comb_j.rows()-1) ){
                    // // while (Sigma_V_comb_j(SV_V_comb_j) > epsilon_rel*sigma0_Aextended && Sigma_V_comb_j(SV_V_comb_j) > machine_eps && SV_V_comb_j<(Sigma_V_comb_j.rows()-1) ){
                    //  SV_V_comb_j++; 
                    // }
                    // if (Sigma_V_comb_j(SV_V_comb_j) > epsilon_rel_basis*Sigma_V_comb_j(0) && Sigma_V_comb_j(SV_V_comb_j) > machine_eps){
                    // // if (Sigma_V_comb_j(SV_V_comb_j) > epsilon_rel*sigma0_Aextended && Sigma_V_comb_j(SV_V_comb_j) > machine_eps){
                    //   SV_V_comb_j=Sigma_V_comb_j.rows();
                    // }
                    //----------------------------------------------------------------------------------------------------------------------------//

                   
                    K_ji_prime=MatrixXd::Zero(jSV,jSV);
                    for (int index=0;index<jSV;index++){
                      K_ji_prime(index,index)=Sigma_j(index);
                    }

                    U_j_prime = svd_P2P_fillin_ji.matrixU().block(0,0,on->ngbr(jNgbr)->nDOF,jSV);
                    V_i_prime = svd_P2P_fillin_ji.matrixV().block(0,0,on->ngbr(iNgbr)->nDOF,jSV);

                    // OLD VERSION (WITHOUT SCALING) //
                    //----------------------------------------------------------------------------------------------------------------------------//
                    // RECOMPRESS THE INTERPOLATION OPERATOR Uj
                    // U_comb_j = MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank+jSV);
                    // U_comb_j.block(0,0,on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank) = on->ngbr(jNgbr)->L2P_operator;
                    // U_comb_j.block(0,on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->nDOF,jSV) = U_j_prime;

                    // JacobiSVD<MatrixXd> svd_U_comb_j(U_comb_j, ComputeThinU | ComputeThinV);
                    // // BDCSVD<MatrixXd> svd_U_comb_j(U_comb_j, ComputeThinU | ComputeThinV);
                    // VectorXd Sigma_U_comb_j = svd_U_comb_j.singularValues();

                    // SV_U_comb_j=0;
                    // while (Sigma_U_comb_j(SV_U_comb_j) > epsilon_rel_basis*Sigma_U_comb_j(0) && Sigma_U_comb_j(SV_U_comb_j) > machine_eps && SV_U_comb_j < (Sigma_U_comb_j.rows()-1)){
                    // // while (Sigma_U_comb_j(SV_U_comb_j) > epsilon_rel*sigma0_Aextended && Sigma_U_comb_j(SV_U_comb_j) > machine_eps && SV_U_comb_j < (Sigma_U_comb_j.rows()-1)){
                    //  SV_U_comb_j++; 
                    // }
                    // if (Sigma_U_comb_j(SV_U_comb_j) > epsilon_rel_basis*Sigma_U_comb_j(0) && Sigma_U_comb_j(SV_U_comb_j) > machine_eps){
                    // // if (Sigma_U_comb_j(SV_U_comb_j) > epsilon_rel*sigma0_Aextended && Sigma_U_comb_j(SV_U_comb_j) > machine_eps){
                    //   SV_U_comb_j=Sigma_U_comb_j.rows();
                    // }

                    // // RECOMPRESS THE ANTERPOLATION OPERATOR Vi
                    // V_comb_i =  MatrixXd::Zero(on->ngbr(iNgbr)->nDOF,on->ngbr(iNgbr)->rank+jSV);
                    // V_comb_i.block(0,0,on->ngbr(iNgbr)->nDOF,on->ngbr(iNgbr)->rank) = on->ngbr(iNgbr)->P2M_operator.transpose();
                    // V_comb_i.block(0,on->ngbr(iNgbr)->rank,on->ngbr(iNgbr)->nDOF,jSV) = V_i_prime;


                    // JacobiSVD<MatrixXd> svd_V_comb_i(V_comb_i, ComputeThinU | ComputeThinV);
                    // // BDCSVD<MatrixXd> svd_V_comb_i(V_comb_i, ComputeThinU | ComputeThinV);
                    // VectorXd Sigma_V_comb_i = svd_V_comb_i.singularValues();


                    // SV_V_comb_i=0;
                    // while (Sigma_V_comb_i(SV_V_comb_i) > epsilon_rel_basis*Sigma_V_comb_i(0) && Sigma_V_comb_i(SV_V_comb_i) > machine_eps && SV_V_comb_i < (Sigma_V_comb_i.rows()-1)){
                    // // while (Sigma_V_comb_i(SV_V_comb_i) > epsilon_rel*sigma0_Aextended && Sigma_V_comb_i(SV_V_comb_i) > machine_eps && SV_V_comb_i < (Sigma_V_comb_i.rows()-1)){
                    //  SV_V_comb_i++; 
                    // }
                    // if (Sigma_V_comb_i(SV_V_comb_i) > epsilon_rel_basis*Sigma_V_comb_i(0) && Sigma_V_comb_i(SV_V_comb_i) > machine_eps){
                    // // if (Sigma_V_comb_i(SV_V_comb_i) > epsilon_rel*sigma0_Aextended && Sigma_V_comb_i(SV_V_comb_i) > machine_eps){
                    //   SV_V_comb_i=Sigma_V_comb_i.rows();
                    // }


                    // // printf("SV_U_comb_i,SV_V_comb_i: %d %d\n",SV_U_comb_i,SV_V_comb_i);
                    // // printf("SV_U_comb_j,SV_V_comb_j: %d %d\n",SV_U_comb_j,SV_V_comb_j);

                    // // MAKE SURE Ui AND Vi HAVE THE SAME RANK
                    // if (SV_U_comb_i != SV_V_comb_i){
                    //   SV_U_comb_i=min(min(on->ngbr(iNgbr)->rank+iSV,on->ngbr(iNgbr)->rank+jSV),max(SV_U_comb_i,SV_V_comb_i));
                    //   // SV_U_comb_i=min(SV_U_comb_i,SV_V_comb_i);
                    //   SV_V_comb_i = SV_U_comb_i;
                    // }

                    // // MAKE SURE Uj AND Vj HAVE THE SAME RANK
                    // if (SV_U_comb_j != SV_V_comb_j){
                    //   SV_U_comb_j=min(min(on->ngbr(jNgbr)->rank+jSV,on->ngbr(jNgbr)->rank+iSV),max(SV_U_comb_j,SV_V_comb_j));
                    //   // SV_U_comb_j=min(SV_U_comb_j,SV_V_comb_j);
                    //   SV_V_comb_j = SV_U_comb_j;
                    // }



                    // // RECOMPRESS Ui
                    // MatrixXd Sigma_diag_U_comb_i=MatrixXd::Zero(SV_U_comb_i,SV_U_comb_i);
                    // for (int index=0;index<SV_U_comb_i;index++){
                    //   Sigma_diag_U_comb_i(index,index)=Sigma_U_comb_i(index);
                    // }
                    // U_recomp_i = svd_U_comb_i.matrixU().block(0,0,on->ngbr(iNgbr)->nDOF,SV_U_comb_i);
                    // R_recomp_U_i = Sigma_diag_U_comb_i*(svd_U_comb_i.matrixV().block(0,0,on->ngbr(iNgbr)->rank+iSV,SV_U_comb_i).transpose());

                    // // RECOMPRESS Uj
                    // MatrixXd Sigma_diag_U_comb_j=MatrixXd::Zero(SV_U_comb_j,SV_U_comb_j);
                    // for (int index=0;index<SV_U_comb_j;index++){
                    //   Sigma_diag_U_comb_j(index,index)=Sigma_U_comb_j(index);
                    // }
                    // U_recomp_j = svd_U_comb_j.matrixU().block(0,0,on->ngbr(jNgbr)->nDOF,SV_U_comb_j);
                    // R_recomp_U_j = Sigma_diag_U_comb_j*(svd_U_comb_j.matrixV().block(0,0,on->ngbr(jNgbr)->rank+jSV,SV_U_comb_j).transpose());
                    
                    // // RECOMPRESS Vi
                    // MatrixXd Sigma_diag_V_comb_i=MatrixXd::Zero(SV_V_comb_i,SV_V_comb_i);
                    // for (int index=0;index<SV_V_comb_i;index++){
                    //   Sigma_diag_V_comb_i(index,index)=Sigma_V_comb_i(index);
                    // }
                    // V_recomp_i = svd_V_comb_i.matrixU().block(0,0,on->ngbr(iNgbr)->nDOF,SV_V_comb_i);
                    // R_recomp_V_i = Sigma_diag_V_comb_i*(svd_V_comb_i.matrixV().block(0,0,on->ngbr(iNgbr)->rank+jSV,SV_V_comb_i).transpose());

                    // // RECOMPRESS Vj 
                    // MatrixXd Sigma_diag_V_comb_j=MatrixXd::Zero(SV_V_comb_j,SV_V_comb_j);
                    // for (int index=0;index<SV_V_comb_j;index++){
                    //   Sigma_diag_V_comb_j(index,index)=Sigma_V_comb_j(index);
                    // }
                    // V_recomp_j = svd_V_comb_j.matrixU().block(0,0,on->ngbr(jNgbr)->nDOF,SV_V_comb_j);
                    // R_recomp_V_j = Sigma_diag_V_comb_j*(svd_V_comb_j.matrixV().block(0,0,on->ngbr(jNgbr)->rank+iSV,SV_V_comb_j).transpose());
                    //----------------------------------------------------------------------------------------------------------------------------//

                    // NEW VERSION (WITH SCALING) //
                    // RECOMPRESS THE SCALED (!) INTERPOLATION OPERATOR Ui
                    U_comb_i = MatrixXd::Zero(on->ngbr(iNgbr)->nDOF,on->ngbr(iNgbr)->rank+iSV);
                    U_comb_i.block(0,0,on->ngbr(iNgbr)->nDOF,on->ngbr(iNgbr)->rank) = on->ngbr(iNgbr)->L2P_operator*on->ngbr(iNgbr)->Sigma_L2P;
                    U_comb_i.block(0,on->ngbr(iNgbr)->rank,on->ngbr(iNgbr)->nDOF,iSV) = U_i_prime*K_ij_prime;

#ifdef DISABLE_JACOBISVD
                    LapackSVD<MatrixXd> svd_U_comb_i(U_comb_i, ComputeThinU | ComputeThinV);
#else
                    JacobiSVD<MatrixXd> svd_U_comb_i(U_comb_i, ComputeThinU | ComputeThinV);
#endif
                    VectorXd Sigma_U_comb_i = svd_U_comb_i.singularValues();

                    SV_U_comb_i=0;
                    while (Sigma_U_comb_i(SV_U_comb_i) > epsilon_rel_basis*Sigma_U_comb_i(0) && Sigma_U_comb_i(SV_U_comb_i) > machine_eps && SV_U_comb_i<(Sigma_U_comb_i.rows()-1) ){
		      SV_U_comb_i++; 
                    }
                    if (Sigma_U_comb_i(SV_U_comb_i) > epsilon_rel_basis*Sigma_U_comb_i(0) && Sigma_U_comb_i(SV_U_comb_i) > machine_eps){
                      SV_U_comb_i=Sigma_U_comb_i.rows();
                    }



                    // RECOMPRESS THE SCALED (!) ANTERPOLATION OPERATOR Vj
                    V_comb_j =  MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank+iSV);
                    V_comb_j.block(0,0,on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank) = on->ngbr(jNgbr)->P2M_operator.transpose()*on->ngbr(jNgbr)->Sigma_P2M;
                    V_comb_j.block(0,on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->nDOF,iSV) = V_j_prime*K_ij_prime;

#ifdef DISABLE_JACOBISVD
                    LapackSVD<MatrixXd> svd_V_comb_j(V_comb_j, ComputeThinU | ComputeThinV);
#else
                    JacobiSVD<MatrixXd> svd_V_comb_j(V_comb_j, ComputeThinU | ComputeThinV);
#endif
                    VectorXd Sigma_V_comb_j = svd_V_comb_j.singularValues();

                    SV_V_comb_j=0;
                    while (Sigma_V_comb_j(SV_V_comb_j) > epsilon_rel_basis*Sigma_V_comb_j(0) && Sigma_V_comb_j(SV_V_comb_j) > machine_eps && SV_V_comb_j<(Sigma_V_comb_j.rows()-1) ){
		      SV_V_comb_j++; 
                    }
                    if (Sigma_V_comb_j(SV_V_comb_j) > epsilon_rel_basis*Sigma_V_comb_j(0) && Sigma_V_comb_j(SV_V_comb_j) > machine_eps){
                      SV_V_comb_j=Sigma_V_comb_j.rows();
                    }

                   
                    // RECOMPRESS THE SCALED (!) INTERPOLATION OPERATOR Uj
                    U_comb_j = MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank+jSV);
                    U_comb_j.block(0,0,on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank) = on->ngbr(jNgbr)->L2P_operator*on->ngbr(jNgbr)->Sigma_L2P;
                    U_comb_j.block(0,on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->nDOF,jSV) = U_j_prime*K_ji_prime;

#ifdef DISABLE_JACOBISVD
                    LapackSVD<MatrixXd> svd_U_comb_j(U_comb_j, ComputeThinU | ComputeThinV);
#else
                    JacobiSVD<MatrixXd> svd_U_comb_j(U_comb_j, ComputeThinU | ComputeThinV);
#endif
                    VectorXd Sigma_U_comb_j = svd_U_comb_j.singularValues();


                    SV_U_comb_j=0;
                    while (Sigma_U_comb_j(SV_U_comb_j) > epsilon_rel_basis*Sigma_U_comb_j(0) && Sigma_U_comb_j(SV_U_comb_j) > machine_eps && SV_U_comb_j < (Sigma_U_comb_j.rows()-1)){
		      // mexPrintf("SV_U_comb_j: %d \n",SV_U_comb_j);
		      SV_U_comb_j++; 
                    }
                    if (Sigma_U_comb_j(SV_U_comb_j) > epsilon_rel_basis*Sigma_U_comb_j(0) && Sigma_U_comb_j(SV_U_comb_j) > machine_eps){
                      SV_U_comb_j=Sigma_U_comb_j.rows();
                    }



                    // RECOMPRESS THE SCALED (!) ANTERPOLATION OPERATOR Vi
                    V_comb_i =  MatrixXd::Zero(on->ngbr(iNgbr)->nDOF,on->ngbr(iNgbr)->rank+jSV);
                    V_comb_i.block(0,0,on->ngbr(iNgbr)->nDOF,on->ngbr(iNgbr)->rank) = on->ngbr(iNgbr)->P2M_operator.transpose()*on->ngbr(iNgbr)->Sigma_P2M;
                    V_comb_i.block(0,on->ngbr(iNgbr)->rank,on->ngbr(iNgbr)->nDOF,jSV) = V_i_prime*K_ji_prime;

#ifdef DISABLE_JACOBISVD
                    LapackSVD<MatrixXd> svd_V_comb_i(V_comb_i, ComputeThinU | ComputeThinV);
#else
                    JacobiSVD<MatrixXd> svd_V_comb_i(V_comb_i, ComputeThinU | ComputeThinV);
#endif
                    VectorXd Sigma_V_comb_i = svd_V_comb_i.singularValues();

                    SV_V_comb_i=0;
                    while (Sigma_V_comb_i(SV_V_comb_i) > epsilon_rel_basis*Sigma_V_comb_i(0) && Sigma_V_comb_i(SV_V_comb_i) > machine_eps && SV_V_comb_i < (Sigma_V_comb_i.rows()-1)){
		      // while (Sigma_V_comb_i(SV_V_comb_i) > epsilon_rel*sigma0_Aextended && Sigma_V_comb_i(SV_V_comb_i) > machine_eps && SV_V_comb_i < (Sigma_V_comb_i.rows()-1)){
		      SV_V_comb_i++; 
                    }
                    if (Sigma_V_comb_i(SV_V_comb_i) > epsilon_rel_basis*Sigma_V_comb_i(0) && Sigma_V_comb_i(SV_V_comb_i) > machine_eps){
                      SV_V_comb_i=Sigma_V_comb_i.rows();
                    }



                    // MAKE SURE Ui AND Vi HAVE THE SAME RANK
                    if (SV_U_comb_i != SV_V_comb_i){
                      SV_U_comb_i=min(min(on->ngbr(iNgbr)->rank+iSV,on->ngbr(iNgbr)->rank+jSV),max(SV_U_comb_i,SV_V_comb_i));
                      SV_V_comb_i = SV_U_comb_i;
                    }

                    // MAKE SURE Uj AND Vj HAVE THE SAME RANK
                    if (SV_U_comb_j != SV_V_comb_j){
                      SV_U_comb_j=min(min(on->ngbr(jNgbr)->rank+jSV,on->ngbr(jNgbr)->rank+iSV),max(SV_U_comb_j,SV_V_comb_j));
                      SV_V_comb_j = SV_U_comb_j;
                    }




                    // RECOMPRESS Ui
                    Sigma_diag_U_comb_i=MatrixXd::Zero(SV_U_comb_i,SV_U_comb_i);
                    for (int index=0;index<SV_U_comb_i;index++){
                      Sigma_diag_U_comb_i(index,index)=Sigma_U_comb_i(index);
                    }
                    U_recomp_i = svd_U_comb_i.matrixU().block(0,0,on->ngbr(iNgbr)->nDOF,SV_U_comb_i);
                    R_recomp_U_i = Sigma_diag_U_comb_i*(svd_U_comb_i.matrixV().block(0,0,on->ngbr(iNgbr)->rank+iSV,SV_U_comb_i).transpose());

                    // RECOMPRESS Uj
                    Sigma_diag_U_comb_j=MatrixXd::Zero(SV_U_comb_j,SV_U_comb_j);
                    for (int index=0;index<SV_U_comb_j;index++){
                      Sigma_diag_U_comb_j(index,index)=Sigma_U_comb_j(index);
                    }
                    U_recomp_j = svd_U_comb_j.matrixU().block(0,0,on->ngbr(jNgbr)->nDOF,SV_U_comb_j);
                    R_recomp_U_j = Sigma_diag_U_comb_j*(svd_U_comb_j.matrixV().block(0,0,on->ngbr(jNgbr)->rank+jSV,SV_U_comb_j).transpose());
                    
                    // RECOMPRESS Vi
                    Sigma_diag_V_comb_i=MatrixXd::Zero(SV_V_comb_i,SV_V_comb_i);
                    for (int index=0;index<SV_V_comb_i;index++){
                      Sigma_diag_V_comb_i(index,index)=Sigma_V_comb_i(index);
                    }
                    V_recomp_i = svd_V_comb_i.matrixU().block(0,0,on->ngbr(iNgbr)->nDOF,SV_V_comb_i);
                    R_recomp_V_i = Sigma_diag_V_comb_i*(svd_V_comb_i.matrixV().block(0,0,on->ngbr(iNgbr)->rank+jSV,SV_V_comb_i).transpose());

                    // RECOMPRESS Vj 
                    Sigma_diag_V_comb_j=MatrixXd::Zero(SV_V_comb_j,SV_V_comb_j);
                    for (int index=0;index<SV_V_comb_j;index++){
                      Sigma_diag_V_comb_j(index,index)=Sigma_V_comb_j(index);
                    }
                    V_recomp_j = svd_V_comb_j.matrixU().block(0,0,on->ngbr(jNgbr)->nDOF,SV_V_comb_j);
                    R_recomp_V_j = Sigma_diag_V_comb_j*(svd_V_comb_j.matrixV().block(0,0,on->ngbr(jNgbr)->rank+iSV,SV_V_comb_j).transpose());


                    

		  }
