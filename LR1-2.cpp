		  iSV=0;
		  // ACA_FullyPivoted(P2P_fillin_ij,U_i_prime,V_j_prime,epsilon,iSV,iSV);
		  ACA_FullyPivoted(P2P_fillin_ij,U_i_prime,V_j_prime,epsilon_rel_fillin,iSV,0);
		  K_ij_prime=MatrixXd::Identity(iSV,iSV);

		  jSV=0;
		  // ACA_FullyPivoted(P2P_fillin_ji,U_j_prime,V_i_prime,epsilon,jSV,jSV);
		  ACA_FullyPivoted(P2P_fillin_ji,U_j_prime,V_i_prime,epsilon_rel_fillin,jSV,0);
		  K_ji_prime=MatrixXd::Identity(jSV,jSV);

		  outfile << "Rank P2P_ij (" << on->id << " " << on->ngbr(iNgbr)->id << " " << on->ngbr(jNgbr)->id << ")" << " " << iSV << endl;
		  outfile << "Rank P2P_ji (" << on->id << " " << on->ngbr(iNgbr)->id << " " << on->ngbr(jNgbr)->id << ")" << " " << jSV << endl;



		  NeedUpdate = (iSV==0 && jSV==0) ? false:true;

		  if (NeedUpdate==true){



                    // RECOMPRESS THE INTERPOLATION OPERATOR Ui
                    U_comb_i = MatrixXd::Zero(on->ngbr(iNgbr)->nDOF,on->ngbr(iNgbr)->rank+iSV);
                    U_comb_i.block(0,0,on->ngbr(iNgbr)->nDOF,on->ngbr(iNgbr)->rank) = on->ngbr(iNgbr)->L2P_operator;
                    U_comb_i.block(0,on->ngbr(iNgbr)->rank,on->ngbr(iNgbr)->nDOF,iSV) = U_i_prime;
                    // RECOMPRESS Ui
                    MatrixXd temp_R_recomp_U_i;
                    SV_U_comb_i=0;
                    // SV_U_comb_i=on->ngbr(iNgbr)->rank; // FIX THE RANK !!!                
                    // ACA_FullyPivoted(U_comb_i,U_recomp_i,temp_R_recomp_U_i,epsilon,SV_U_comb_i,SV_U_comb_i);
                    ACA_FullyPivoted(U_comb_i,U_recomp_i,temp_R_recomp_U_i,epsilon_rel_basis,SV_U_comb_i,0);
                    R_recomp_U_i=temp_R_recomp_U_i.transpose();

                    // printf("iSV, jSV: %d %d \n",iSV,jSV);

                    // printf("on->ngbr(iNgbr)->rank: %d\n",on->ngbr(iNgbr)->rank);
                    // printf("on->ngbr(jNgbr)->rank: %d\n",on->ngbr(jNgbr)->rank);


                    
                    // RECOMPRESS THE ANTERPOLATION OPERATOR Vi
                    V_comb_i =  MatrixXd::Zero(on->ngbr(iNgbr)->nDOF,on->ngbr(iNgbr)->rank+jSV);
                    V_comb_i.block(0,0,on->ngbr(iNgbr)->nDOF,on->ngbr(iNgbr)->rank) = on->ngbr(iNgbr)->P2M_operator.transpose();
                    V_comb_i.block(0,on->ngbr(iNgbr)->rank,on->ngbr(iNgbr)->nDOF,jSV) = V_i_prime;
                    // RECOMPRESS Vi
                    MatrixXd temp_R_recomp_V_i;
                    SV_V_comb_i=0;
                    // SV_V_comb_i=on->ngbr(iNgbr)->rank; // FIX THE RANK !!!
                    ACA_FullyPivoted(V_comb_i,V_recomp_i,temp_R_recomp_V_i,epsilon_rel_basis,SV_V_comb_i,SV_U_comb_i);
                    R_recomp_V_i=temp_R_recomp_V_i.transpose();

                    // cout << "on->ngbr(jNgbr)->L2P_operator.rows()" << on->ngbr(jNgbr)->L2P_operator.rows() << endl;
                    // cout << "on->ngbr(jNgbr)->L2P_operator.cols()" << on->ngbr(jNgbr)->L2P_operator.cols() << endl;

                    // cout << "U_j_prime.rows()" << U_j_prime.rows() << endl;
                    // cout << "U_j_prime.cols()" << U_j_prime.cols() << endl;


                    // RECOMPRESS THE INTERPOLATION OPERATOR Uj
                    U_comb_j = MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank+jSV);
                    U_comb_j.block(0,0,on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank) = on->ngbr(jNgbr)->L2P_operator;
                    U_comb_j.block(0,on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->nDOF,jSV) = U_j_prime;
                    // RECOMPRESS Uj
                    MatrixXd temp_R_recomp_U_j;
                    SV_U_comb_j=0;


                    // SV_U_comb_j=on->ngbr(jNgbr)->rank; // FIX THE RANK !!!
                    // ACA_FullyPivoted(U_comb_j,U_recomp_j,temp_R_recomp_U_j,epsilon,SV_U_comb_j,SV_U_comb_j);
                    ACA_FullyPivoted(U_comb_j,U_recomp_j,temp_R_recomp_U_j,epsilon_rel_basis,SV_U_comb_j,0);
                    R_recomp_U_j=temp_R_recomp_U_j.transpose();

                    
                    // RECOMPRESS THE ANTERPOLATION OPERATOR Vj
                    V_comb_j =  MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank+iSV);
                    V_comb_j.block(0,0,on->ngbr(jNgbr)->nDOF,on->ngbr(jNgbr)->rank) = on->ngbr(jNgbr)->P2M_operator.transpose();
                    V_comb_j.block(0,on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->nDOF,iSV) = V_j_prime;
                    // RECOMPRESS Vj 
                    MatrixXd temp_R_recomp_V_j;
                    SV_V_comb_j=0;
                    // SV_V_comb_j=on->ngbr(jNgbr)->rank; // FIX THE RANK !!!
                    ACA_FullyPivoted(V_comb_j,V_recomp_j,temp_R_recomp_V_j,epsilon_rel_basis,SV_V_comb_j,SV_U_comb_j);
                    R_recomp_V_j=temp_R_recomp_V_j.transpose();

                    

                    // MAKE SURE Ui AND Vi HAVE THE SAME RANK
                    if (SV_U_comb_i != SV_V_comb_i){
                      // printf("SV_U_comb_i,SV_V_comb_i: %d %d \n",SV_U_comb_i,SV_V_comb_i);
                      int temprank=min(min(on->ngbr(iNgbr)->rank+iSV,on->ngbr(iNgbr)->rank+jSV),max(SV_U_comb_i,SV_V_comb_i));

                      if (temprank==SV_U_comb_i){
                        SV_V_comb_i=temprank;

                        ACA_FullyPivoted(V_comb_i,V_recomp_i,temp_R_recomp_V_i,epsilon_rel_basis,SV_V_comb_i,SV_V_comb_i);
                        R_recomp_V_i=temp_R_recomp_V_i.transpose();
                      }
                      else if(temprank==SV_V_comb_i){
                        SV_U_comb_i=temprank;

                        ACA_FullyPivoted(U_comb_i,U_recomp_i,temp_R_recomp_U_i,epsilon_rel_basis,SV_U_comb_i,SV_U_comb_i);
                        R_recomp_U_i=temp_R_recomp_U_i.transpose();
                      }
                      else{
                        SV_U_comb_i=temprank;
                        SV_V_comb_i=temprank;
                        
                        ACA_FullyPivoted(U_comb_i,U_recomp_i,temp_R_recomp_U_i,epsilon_rel_basis,SV_U_comb_i,SV_U_comb_i);
                        R_recomp_U_i=temp_R_recomp_U_i.transpose();
                        ACA_FullyPivoted(V_comb_i,V_recomp_i,temp_R_recomp_V_i,epsilon_rel_basis,SV_V_comb_i,SV_V_comb_i);
                        R_recomp_V_i=temp_R_recomp_V_i.transpose();
                      }
                    }


		    // MAKE SURE Uj AND Vj HAVE THE SAME RANK
                    if (SV_U_comb_j != SV_V_comb_j){
                      // printf("SV_U_comb_j,SV_V_comb_j: %d %d \n",SV_U_comb_j,SV_V_comb_j);
                      int temprank=min(min(on->ngbr(jNgbr)->rank+jSV,on->ngbr(jNgbr)->rank+iSV),max(SV_U_comb_j,SV_V_comb_j));

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

