
		  // cout << "P2P_rsvd" << endl;

//20200108                  start_SVD = clock(); 


		  // int rank_RSVD = max(rank_RSVD_min,(int) min(P2P_fillin_ij.rows(),P2P_fillin_ij.cols())/10);
		  // if (rank_RSVD > min(P2P_fillin_ij.rows(),P2P_fillin_ij.cols())){
		  //    rank_RSVD = min(P2P_fillin_ij.rows(),P2P_fillin_ij.cols());
		  // }
		  int rank_RSVD = rank_RSVD_min;

                                       
		  // SVD OF P2P_ij
		  MatrixXd Y_P2P_ij = - eta_si[jNgbr].transpose()*p2pto[iNgbr];
		  gram_schmidt(Y_P2P_ij);

		  MatrixXd B_P2P_ij = - on->P2P_operator[iNgbr]*(eta_si[jNgbr]*Y_P2P_ij);
		  MatrixXd Z_P2P_ij = B_P2P_ij * P;
		  // MatrixXd Z_P2P_ij = - on->P2P_operator[iNgbr]*(eta_si[jNgbr]*(Y_P2P_ij*P));
		  gram_schmidt(Z_P2P_ij);

		  MatrixXd C_P2P_ij = Z_P2P_ij.transpose()*B_P2P_ij;
		  // MatrixXd C_P2P_ij =  (-Z_P2P_ij.transpose()*on->P2P_operator[iNgbr])*(eta_si[jNgbr]*Y_P2P_ij);
#ifdef DISABLE_JACOBISVD
		  LapackSVD<MatrixXd> svd_C_P2P_ij(C_P2P_ij, ComputeThinU | ComputeThinV);
#else
		  JacobiSVD<MatrixXd> svd_C_P2P_ij(C_P2P_ij, ComputeThinU | ComputeThinV);
#endif
		  VectorXd Sigma_i = svd_C_P2P_ij.singularValues();


		  // RedSVD<MatrixXd> rsvd_P2P_fillin_ij(P2P_fillin_ij,rank_RSVD);
		  // VectorXd Sigma_i = rsvd_P2P_fillin_ij.singularValues();
		  RedSVD<MatrixXd> rsvd_P2P_fillin_ij;


                    
		  iSV=0;
		  while (Sigma_i(iSV) > epsilon_rel_fillin*sigma0_Aextended && Sigma_i(iSV) > machine_eps && iSV< (rank_RSVD-1)){
		    iSV++; 
		  }
		  if (Sigma_i(iSV) > epsilon_rel_fillin*sigma0_Aextended && Sigma_i(iSV) > machine_eps){

		    P2P_fillin_ij= - on->P2P_operator[iNgbr]*eta_si[jNgbr];

		    // SVD OF P2P_ij
		    rsvd_P2P_fillin_ij.compute(P2P_fillin_ij,min(P2P_fillin_ij.rows(),P2P_fillin_ij.cols()));
		    Sigma_i = rsvd_P2P_fillin_ij.singularValues();

		    iSV=0;
		    while (Sigma_i(iSV) > epsilon_rel_fillin*sigma0_Aextended && Sigma_i(iSV) > machine_eps && iSV< (min(P2P_fillin_ij.rows(),P2P_fillin_ij.cols())-1)){
		      iSV++; 
		    }
		    if (Sigma_i(iSV) > epsilon_rel_fillin*sigma0_Aextended && Sigma_i(iSV) > machine_eps){
		      iSV=min(P2P_fillin_ij.rows(),P2P_fillin_ij.cols());
		    }
		  }



		  // SVD OF P2P_ji
		  MatrixXd Y_P2P_ji = - eta_si[iNgbr].transpose()*p2pto[jNgbr];
		  gram_schmidt(Y_P2P_ji);

		  MatrixXd B_P2P_ji = - on->P2P_operator[jNgbr]*(eta_si[iNgbr]*Y_P2P_ji);
		  MatrixXd Z_P2P_ji = B_P2P_ji * P;
		  // MatrixXd Z_P2P_ji =  - on->P2P_operator[jNgbr]*(eta_si[iNgbr]*(Y_P2P_ji*P));
		  gram_schmidt(Z_P2P_ji);

		  MatrixXd C_P2P_ji = Z_P2P_ji.transpose()*B_P2P_ji;
		  // MatrixXd C_P2P_ji =  (-Z_P2P_ji.transpose()*on->P2P_operator[jNgbr])*(eta_si[iNgbr]*Y_P2P_ji);

#ifdef DISABLE_JACOBISVD
                  LapackSVD<MatrixXd> svd_C_P2P_ji(C_P2P_ji, ComputeThinU | ComputeThinV);
#else
		  JacobiSVD<MatrixXd> svd_C_P2P_ji(C_P2P_ji, ComputeThinU | ComputeThinV);
#endif
		  VectorXd Sigma_j = svd_C_P2P_ji.singularValues();

		  // RedSVD<MatrixXd> rsvd_P2P_fillin_ji(P2P_fillin_ji,rank_RSVD);
		  // VectorXd Sigma_j = rsvd_P2P_fillin_ji.singularValues();
		  RedSVD<MatrixXd> rsvd_P2P_fillin_ji;



		  jSV=0;
		  while (Sigma_j(jSV) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j(jSV) > machine_eps && jSV< (rank_RSVD-1)){
		    jSV++; 
		  }
		  if (Sigma_j(jSV) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j(jSV) > machine_eps){

		    P2P_fillin_ji= - on->P2P_operator[jNgbr]*eta_si[iNgbr];

		    // SVD OF P2P_ji
		    rsvd_P2P_fillin_ji.compute(P2P_fillin_ji,min(P2P_fillin_ji.rows(),P2P_fillin_ji.cols()));
		    Sigma_j = rsvd_P2P_fillin_ji.singularValues();

		    jSV=0;
		    while (Sigma_j(jSV) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j(jSV) > machine_eps && jSV< (min(P2P_fillin_ji.rows(),P2P_fillin_ji.cols())-1)){
		      jSV++; 
		    }
		    if (Sigma_j(jSV) > epsilon_rel_fillin*sigma0_Aextended && Sigma_j(jSV) > machine_eps){
		      jSV=min(P2P_fillin_ji.rows(),P2P_fillin_ji.cols());
		    }
		  }



		  if (iSV != jSV){
		    int iSV_temp=iSV;
		    int jSV_temp=jSV;


		    iSV=min(min(on->ngbr(iNgbr)->nDOF,on->ngbr(jNgbr)->nDOF),max(iSV,jSV));
		    jSV = iSV;


		    if (iSV>rank_RSVD && iSV > iSV_temp){

		      if (iSV_temp < rank_RSVD){
			P2P_fillin_ij= - on->P2P_operator[iNgbr]*eta_si[jNgbr];
		      }

		      // RedSVD<MatrixXd> rsvd_P2P_fillin_ij_new(P2P_fillin_ij,iSV);
		      rsvd_P2P_fillin_ij.compute(P2P_fillin_ij,iSV);
		      Sigma_i = rsvd_P2P_fillin_ij.singularValues();
		    }
		    if (jSV>rank_RSVD && jSV > jSV_temp){
		      if (jSV_temp < rank_RSVD){
			P2P_fillin_ji= - on->P2P_operator[jNgbr]*eta_si[iNgbr];
		      }
		      // RedSVD<MatrixXd> rsvd_P2P_fillin_ji_new(P2P_fillin_ji,jSV);
		      rsvd_P2P_fillin_ji.compute(P2P_fillin_ji,jSV);
		      Sigma_j = rsvd_P2P_fillin_ji.singularValues();
		    }

                     

                    while (iSV>0){
                      if (Sigma_i(jSV-1) < 1000*machine_eps || Sigma_j(iSV-1)< 1000*machine_eps){
                        iSV--;
                        jSV--;
                      }
                      else{
                        break;
                      }
                    }
                   
                     
                    
		  }

//20200108		  stop_SVD = clock();    
//20200108		  TimeSVD += double(stop_SVD-start_SVD)/double(CLOCKS_PER_SEC);

                   
		  NeedUpdate = (iSV==0) ? false:true;


		  if (NeedUpdate==true && Delay_update==false){

		    // cout << "NeedUpdate" << endl;    

		    // printf("iSV,jSV: %d %d\n",iSV,jSV);


		    K_ij_prime=MatrixXd::Zero(iSV,iSV);
		    for (int index=0;index<iSV;index++){
		      K_ij_prime(index,index)=Sigma_i(index);
		    }

		    if (iSV>rank_RSVD){
		      // printf("iSV: %d\n",iSV);
		      U_i_prime = rsvd_P2P_fillin_ij.matrixU().block(0,0,on->ngbr(iNgbr)->nDOF,iSV);
		      V_j_prime = rsvd_P2P_fillin_ij.matrixV().block(0,0,on->ngbr(jNgbr)->nDOF,iSV);
		    }
		    else{
		      U_i_prime = Z_P2P_ij.block(0,0,on->ngbr(iNgbr)->nDOF,rank_RSVD)*svd_C_P2P_ij.matrixU().block(0,0,rank_RSVD,iSV);
		      V_j_prime = Y_P2P_ij.block(0,0,on->ngbr(jNgbr)->nDOF,rank_RSVD)*svd_C_P2P_ij.matrixV().block(0,0,rank_RSVD,iSV);
		    }



		    K_ji_prime=MatrixXd::Zero(jSV,jSV);
		    for (int index=0;index<jSV;index++){
		      K_ji_prime(index,index)=Sigma_j(index);
		    }

		    if (jSV>rank_RSVD){
		      U_j_prime = rsvd_P2P_fillin_ji.matrixU().block(0,0,on->ngbr(jNgbr)->nDOF,jSV);
		      V_i_prime = rsvd_P2P_fillin_ji.matrixV().block(0,0,on->ngbr(iNgbr)->nDOF,jSV);
		    }
		    else{
		      U_j_prime = Z_P2P_ji.block(0,0,on->ngbr(jNgbr)->nDOF,rank_RSVD)*svd_C_P2P_ji.matrixU().block(0,0,rank_RSVD,jSV);
		      V_i_prime = Y_P2P_ji.block(0,0,on->ngbr(iNgbr)->nDOF,rank_RSVD)*svd_C_P2P_ji.matrixV().block(0,0,rank_RSVD,jSV);
		    } 


		    //20200108		    start_ACA = clock();

		    // OLD VERSION (WITHOUT SCALING)
		    // -----------------------------------------------------------------------------------------------------------------------//
		    // // RECOMPRESS THE INTERPOLATION OPERATOR Ui
		    // U_comb_i = MatrixXd::Zero(on->ngbr(iNgbr)->nDOF,on->ngbr(iNgbr)->rank+iSV);
		    // U_comb_i.block(0,0,on->ngbr(iNgbr)->nDOF,on->ngbr(iNgbr)->rank) = on->ngbr(iNgbr)->L2P_operator;
		    // U_comb_i.block(0,on->ngbr(iNgbr)->rank,on->ngbr(iNgbr)->nDOF,iSV) = U_i_prime;
		    // // RECOMPRESS Ui
		    // MatrixXd temp_R_recomp_U_i;
		    // SV_U_comb_i=0;
		    // // SV_U_comb_i=on->ngbr(iNgbr)->rank; // FIX THE RANK !!!                
		    // // ACA_FullyPivoted(U_comb_i,U_recomp_i,temp_R_recomp_U_i,epsilon,SV_U_comb_i,SV_U_comb_i);

                      
		    // ACA_FullyPivoted(U_comb_i,U_recomp_i,temp_R_recomp_U_i,epsilon_rel_basis,SV_U_comb_i,0);
		    // R_recomp_U_i=temp_R_recomp_U_i.transpose();

                      
		    // // RECOMPRESS THE ANTERPOLATION OPERATOR Vi
		    // V_comb_i =  MatrixXd::Zero(on->ngbr(iNgbr)->nDOF,on->ngbr(iNgbr)->rank+jSV);
		    // V_comb_i.block(0,0,on->ngbr(iNgbr)->nDOF,on->ngbr(iNgbr)->rank) = on->ngbr(iNgbr)->P2M_operator.transpose();
		    // V_comb_i.block(0,on->ngbr(iNgbr)->rank,on->ngbr(iNgbr)->nDOF,jSV) = V_i_prime;
		    // // RECOMPRESS Vi
		    // MatrixXd temp_R_recomp_V_i;
		    // SV_V_comb_i=0;


		    // // SV_V_comb_i=on->ngbr(iNgbr)->rank; // FIX THE RANK !!!
		    // ACA_FullyPivoted(V_comb_i,V_recomp_i,temp_R_recomp_V_i,epsilon_rel_basis,SV_V_comb_i,SV_U_comb_i);
		    // R_recomp_V_i=temp_R_recomp_V_i.transpose();


		    // // RECOMPRESS THE INTERPOLATION OPERATOR Uj
		    // U_comb_j = MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank+jSV);
		    // U_comb_j.block(0,0,on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank) = on->ngbr(jNgbr)->L2P_operator;
		    // U_comb_j.block(0,on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->nDOF,jSV) = U_j_prime;
		    // // RECOMPRESS Uj
		    // MatrixXd temp_R_recomp_U_j;
		    // SV_U_comb_j=0;


		    // // SV_U_comb_j=on->ngbr(jNgbr)->rank; // FIX THE RANK !!!
		    // // ACA_FullyPivoted(U_comb_j,U_recomp_j,temp_R_recomp_U_j,epsilon,SV_U_comb_j,SV_U_comb_j);
		    // ACA_FullyPivoted(U_comb_j,U_recomp_j,temp_R_recomp_U_j,epsilon_rel_basis,SV_U_comb_j,0);
		    // R_recomp_U_j=temp_R_recomp_U_j.transpose();

                      
		    // // RECOMPRESS THE ANTERPOLATION OPERATOR Vj
		    // V_comb_j =  MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank+iSV);

                     
		    // V_comb_j.block(0,0,on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank) = on->ngbr(jNgbr)->P2M_operator.transpose();
		    // V_comb_j.block(0,on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->nDOF,iSV) = V_j_prime;

		    // // RECOMPRESS Vj 
		    // MatrixXd temp_R_recomp_V_j;
		    // SV_V_comb_j=0;


		    // // SV_V_comb_j=on->ngbr(jNgbr)->rank; // FIX THE RANK !!!
		    // ACA_FullyPivoted(V_comb_j,V_recomp_j,temp_R_recomp_V_j,epsilon_rel_basis,SV_V_comb_j,SV_U_comb_j);
		    // R_recomp_V_j=temp_R_recomp_V_j.transpose();



                      

		    // // MAKE SURE Ui AND Vi HAVE THE SAME RANK
		    // if (SV_U_comb_i != SV_V_comb_i){
		    //   // printf("resize P2P i %d %d \n",SV_U_comb_i,SV_V_comb_i);



		    //   // printf("SV_U_comb_i,SV_V_comb_i: %d %d \n",SV_U_comb_i,SV_V_comb_i);
		    //   int temprank=min(min(on->ngbr(iNgbr)->rank+iSV,on->ngbr(iNgbr)->rank+jSV),max(SV_U_comb_i,SV_V_comb_i));
		    //   // int temprank=min(min(on->ngbr(iNgbr)->rank+iSV,on->ngbr(iNgbr)->rank+jSV),min(SV_U_comb_i,SV_V_comb_i));

		    //   if (temprank==SV_U_comb_i){
		    //     SV_V_comb_i=temprank;

		    //     ACA_FullyPivoted(V_comb_i,V_recomp_i,temp_R_recomp_V_i,epsilon_rel_basis,SV_V_comb_i,SV_V_comb_i);
		    //     R_recomp_V_i=temp_R_recomp_V_i.transpose();
		    //   }
		    //   else if(temprank==SV_V_comb_i){
		    //     SV_U_comb_i=temprank;

		    //     ACA_FullyPivoted(U_comb_i,U_recomp_i,temp_R_recomp_U_i,epsilon_rel_basis,SV_U_comb_i,SV_U_comb_i);
		    //     R_recomp_U_i=temp_R_recomp_U_i.transpose();
		    //   }
		    //   else{
		    //     SV_U_comb_i=temprank;
		    //     SV_V_comb_i=temprank;
                          
		    //     ACA_FullyPivoted(U_comb_i,U_recomp_i,temp_R_recomp_U_i,epsilon_rel_basis,SV_U_comb_i,SV_U_comb_i);
		    //     R_recomp_U_i=temp_R_recomp_U_i.transpose();
                          
		    //     ACA_FullyPivoted(V_comb_i,V_recomp_i,temp_R_recomp_V_i,epsilon_rel_basis,SV_V_comb_i,SV_V_comb_i);
		    //     R_recomp_V_i=temp_R_recomp_V_i.transpose();
		    //   }
		    // }


		    //  // MAKE SURE Uj AND Vj HAVE THE SAME RANK
		    // if (SV_U_comb_j != SV_V_comb_j){
		    //   // printf("SV_U_comb_j,SV_V_comb_j: %d %d \n",SV_U_comb_j,SV_V_comb_j);
		    //   int temprank=min(min(on->ngbr(jNgbr)->rank+jSV,on->ngbr(jNgbr)->rank+iSV),max(SV_U_comb_j,SV_V_comb_j));
		    //   // int temprank=min(min(on->ngbr(jNgbr)->rank+jSV,on->ngbr(jNgbr)->rank+iSV),min(SV_U_comb_j,SV_V_comb_j));

		    //   if (temprank==SV_U_comb_j){
		    //     SV_V_comb_j=temprank;

		    //     ACA_FullyPivoted(V_comb_j,V_recomp_j,temp_R_recomp_V_j,epsilon_rel_basis,SV_V_comb_j,SV_V_comb_j);
		    //     R_recomp_V_j=temp_R_recomp_V_j.transpose();
		    //   }
		    //   else if(temprank==SV_V_comb_j){
		    //     SV_U_comb_j=temprank;

		    //     ACA_FullyPivoted(U_comb_j,U_recomp_j,temp_R_recomp_U_j,epsilon_rel_basis,SV_U_comb_j,SV_U_comb_j);
		    //     R_recomp_U_j=temp_R_recomp_U_j.transpose();
		    //   }
		    //   else{
		    //     SV_U_comb_j=temprank;
		    //     SV_V_comb_j=temprank;
                          
		    //     ACA_FullyPivoted(U_comb_j,U_recomp_j,temp_R_recomp_U_j,epsilon_rel_basis,SV_U_comb_j,SV_U_comb_j);
		    //     R_recomp_U_j=temp_R_recomp_U_j.transpose();
                          
		    //     ACA_FullyPivoted(V_comb_j,V_recomp_j,temp_R_recomp_V_j,epsilon_rel_basis,SV_V_comb_j,SV_V_comb_j);
		    //     R_recomp_V_j=temp_R_recomp_V_j.transpose();
		    //   }
		    // }
		    // -----------------------------------------------------------------------------------------------------------------------//

		    // NEW VERSION (WITH SCALING)

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


		    //20200108		    stop_ACA = clock();    
		    //20200108		    Time_ACA += double(stop_ACA-start_ACA)/double(CLOCKS_PER_SEC);


		  }

