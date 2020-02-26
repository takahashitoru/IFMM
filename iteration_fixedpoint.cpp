#include "iteration_fixedpoint.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;
using namespace Eigen;



void iteration_fixedpoint(H2_3D_Compute<myKernel> &compute, myKernel &Atree, vector3 * point_fmm, int &N, dvec &x, dvec &x0, dvec & b, int & max_iter, double & tol, IFMM_Matrix & ifmm, bool &precond){
	// double &L, int &l, int &n_FMM, double &eps, int &m, int & use_chebyshev, 

	VectorXd b_ = Map<VectorXd>(&b(0),N);
	VectorXd x0_ = Map<VectorXd>(&x0(0),N);
	VectorXd x_old=x0_;
	VectorXd x_new=x0_;

	dvec x_temp(N);
	dvec RHS_ifmm(N);


	// PREPARE FOR FMM-MULTIPLICATION
	double *x_fmm = new double[N]; // Source array
	double *b_fmm = new double[N]; // Field array (BBFMM calculation)

	for (int k =0; k<N;k++){
		x_fmm[k] = x0_(k);
    }
	// H2_3D_Compute<myKernel> compute(&Btree, point_fmm, point_fmm, N, N, x_fmm, m, b_fmm);
	compute.FMMCompute(&(Atree.tree),&point_fmm[0],&point_fmm[0],&x_fmm[0],&(Atree.K[0]),&(Atree.U[0]),&(Atree.VT[0]),&(Atree.Tkz[0]),&(Atree.Ktable[0]),
                       	   &(Atree.Kweights[0]),&(Atree.Cweights[0]),(double) Atree.homogen,&(Atree.cutoff),(int) Atree.n, &(Atree.dof[0]), &b_fmm[0],
                       	   (int) Atree.use_chebyshev,&(Atree.p_r2c[0]),&(Atree.p_c2r[0]));

	int i = 0;
  	VectorXd err=VectorXd::Zero(max_iter);
  	err(i)=(Map<VectorXd>(&b_fmm[0],N) - b_).norm()/b_.norm();
  	printf("Initial error:%20.5e\n",err(i));


  	// compute2.FMMCompute(&(Atree.tree),&point_fmm[0],&point_fmm[0],&qtest[0],&(Atree.K[0]),&(Atree.U[0]),&(Atree.VT[0]),&(Atree.Tkz[0]),&(Atree.Ktable[0]),
   //                     &(Atree.Kweights[0]),&(Atree.Cweights[0]),(double) Atree.homogen,&(Atree.cutoff),(int) Atree.n,&(Atree.dof[0]),&b_fmmtest[0],
   //                     (int) Atree.use_chebyshev,&(Atree.p_r2c[0]),&(Atree.p_c2r[0]));
	  // b_fmmtest_ = Map<VectorXd>(&b_fmmtest[0],N); 
	  // diff_b_fmm = b_ - b_fmmtest_;
	  // relative_error_b_fmm = diff_b_fmm.norm()/b_.norm();
	  // std::cout << "Test E:" << relative_error_b_fmm << std::endl;


	  // 	compute2.FMMCompute(&(Atree.tree),&point_fmm[0],&point_fmm[0],&qtest[0],&(Atree.K[0]),&(Atree.U[0]),&(Atree.VT[0]),&(Atree.Tkz[0]),&(Atree.Ktable[0]),
   //                     &(Atree.Kweights[0]),&(Atree.Cweights[0]),(double) Atree.homogen,&(Atree.cutoff),(int) Atree.n,&(Atree.dof[0]),&b_fmmtest[0],
   //                     (int) Atree.use_chebyshev,&(Atree.p_r2c[0]),&(Atree.p_c2r[0]));
	  // b_fmmtest_ = Map<VectorXd>(&b_fmmtest[0],N); 
	  // diff_b_fmm = b_ - b_fmmtest_;
	  // relative_error_b_fmm = diff_b_fmm.norm()/b_.norm();
	  // std::cout << "Test Y:" << relative_error_b_fmm << std::endl;

  	ofstream outfile;
    char filename [100];


	while (i < max_iter-1 && err(i)>=tol){

		printf("Iteration: %d\n",i);


	  	for (int k=0; k<N;k++)
	  		RHS_ifmm(k)=b(k)-b_fmm[k];

		// ifmm.setRHS(b_fmm_); //RESET THE RIGHT-HAND SIDE
		ifmm.setRHS(RHS_ifmm); //RESET THE RIGHT-HAND SIDE
	    ifmm.elimination_reuse(); // RE-USE THE FACTORIZATION 
	    ifmm.substitution(x_temp,outfile); // PERFORM THE SUBSTITUTION

		// x_new = x0_ + x_old - Map<VectorXd>(&x_temp(0),N);
		if (precond==true){
			x_new = x_old + Map<VectorXd>(&x_temp(0),N);
		}
		else{
			x_new = x_old + Map<VectorXd>(&RHS_ifmm(0),N);
		}


		for (int i =0; i<N;i++)
	    	x_fmm[i] = x_new(i);

	    compute.FMMCompute(&(Atree.tree),&point_fmm[0],&point_fmm[0],&x_fmm[0],&(Atree.K[0]),&(Atree.U[0]),&(Atree.VT[0]),&(Atree.Tkz[0]),&(Atree.Ktable[0]),
                       	   &(Atree.Kweights[0]),&(Atree.Cweights[0]),(double) Atree.homogen,&(Atree.cutoff),(int) Atree.n, &(Atree.dof[0]), &b_fmm[0],
                       	   (int) Atree.use_chebyshev,&(Atree.p_r2c[0]),&(Atree.p_c2r[0]));


		err(i+1)=(Map<VectorXd>(&b_fmm[0],N) - b_).norm()/b_.norm();
		printf("Error: %20.5e\n",err(i+1));
		i++;


		// relvar(i+1) = (x_new-x_old).norm()/x_old.norm();
		// printf("relvar(i+1): %5.5e\n",relvar(i+1));


		x_old=x_new;

	}


	// ITERATION COMPLETED
	for (int i =0; i<N;i++)
    	x(i)=x_new(i);

    for (int i =0; i<N;i++){
    	b_fmm[i] = 0.0;
    	x_fmm[i] = x(i);
    }
    
 //    compute.FMMCompute(&(Atree.tree),&point_fmm[0],&point_fmm[0],&x_fmm[0],&(Atree.K[0]),&(Atree.U[0]),&(Atree.VT[0]),&(Atree.Tkz[0]),&(Atree.Ktable[0]),
 //                       &(Atree.Kweights[0]),&(Atree.Cweights[0]),(double) Atree.homogen,&(Atree.cutoff),(int) Atree.n,&(Atree.dof[0]),&b_fmm[0],
 //                       (int) Atree.use_chebyshev,&(Atree.p_r2c[0]),&(Atree.p_c2r[0]));
	// VectorXd diff_fixedpoint_fmm = b_ - Map<VectorXd>(&b_fmm[0],N);
	//  printf("test test: %20.5e\n",diff_fixedpoint_fmm.norm()/b_.norm());



	//  VectorXd test2= b_ - A*Map<VectorXd>(&x(0),N);
	//  printf("test2 : %20.5e\n",test2.norm()/b_.norm());


	//  double *q  = new double[N*m];// Field array (BBFMM calculation)
	//  double *test3_ = new double[N*m];  // Source array
	//  for (int i =0; i<N;i++){
 //    	q[i] = x(i);
 //    	test3_[i]=0;
 //     }
	//  H2_3D_Compute<myKernel> compute3(&Atree, point_fmm, point_fmm, N, N, q ,m, test3_);
	//  VectorXd test3 = b_ - Map<VectorXd>(&test3_[0],N);
	//  printf("test3: %20.5e\n",test3.norm()/b_.norm());

	//  for (int k=0; k<10; k++){
	// 	printf("test3_[k],b_(k): %5.5e %5.5e\n",test3_[k],b_(k));
	// }


	// for (int i =0; i<N;i++){
 //    	q[i] = x(i);
 //    	test3_[i]=0;
 //     }
	// compute3.FMMCompute(&(Atree.tree),&point_fmm[0],&point_fmm[0],&q[0],&(Atree.K[0]),&(Atree.U[0]),&(Atree.VT[0]),&(Atree.Tkz[0]),&(Atree.Ktable[0]),
 //                       &(Atree.Kweights[0]),&(Atree.Cweights[0]),(double) Atree.homogen,&(Atree.cutoff),(int) Atree.n,&(Atree.dof[0]),&test3_[0],
 //                       (int) Atree.use_chebyshev,&(Atree.p_r2c[0]),&(Atree.p_c2r[0]));
	// test3 = b_ - Map<VectorXd>(&test3_[0],N);
	// printf("test3: %20.5e\n",test3.norm()/b_.norm());


	// compute.FMMCompute(&(Atree.tree),&point_fmm[0],&point_fmm[0],&x_fmm[0],&(Atree.K[0]),&(Atree.U[0]),&(Atree.VT[0]),&(Atree.Tkz[0]),&(Atree.Ktable[0]),
 //                       &(Atree.Kweights[0]),&(Atree.Cweights[0]),(double) Atree.homogen,&(Atree.cutoff),(int) Atree.n,&(Atree.dof[0]),&b_fmm[0],
 //                       (int) Atree.use_chebyshev,&(Atree.p_r2c[0]),&(Atree.p_c2r[0]));
	//  diff_fixedpoint_fmm = b_ - Map<VectorXd>(&b_fmm[0],N);
	//  printf("test X: %20.5e\n",diff_fixedpoint_fmm.norm()/b_.norm());





	//     myKernel Dtree(2.0,2, 2, 1e-9, 1);
	// 	Dtree.buildFMMTree(); // BUILD THE FMM-TREE

	// 	for (int i =0; i<N;i++){
	//     	b_fmm[i] = 0.0;
 //    	}
	// 	H2_3D_Compute<myKernel> computeB(&Dtree, point_fmm, point_fmm, N, N, x_fmm ,m, b_fmm);
	// 	VectorXd diff_c = b_ - Map<VectorXd>(&b_fmm[0],N);
	// 	printf("diff_c: %20.5e\n",diff_c.norm()/b_.norm());

	// 	for (int k=0;k<10;k++){
	// 		// printf("point_fmm[k].x,point_fmm[k].y,point_fmm[k].z: %f %f %f \n", point_fmm[k].x,point_fmm[k].y,point_fmm[k].z);
	// 		// printf("Atree.K[k]: %f \n", Atree.K[k]);
	// 		// printf("Atree.U[k]: %f \n", Atree.U[k]);
	// 		// printf("Atree.VT[k]: %f \n", Atree.VT[k]);
	// 		// printf("Atree.Tkz[k]: %f \n", Atree.Tkz[k]);
	// 		printf("x_fmm[k]: %f \n", x_fmm[k]);


	// 	}


	// 	for (int i =0; i<N;i++){
	//     	b_fmm[i] = 0.0;
 //    	}
	// 	compute.FMMCompute(&(Atree.tree),&point_fmm[0],&point_fmm[0],&x_fmm[0],&(Atree.K[0]),&(Atree.U[0]),&(Atree.VT[0]),&(Atree.Tkz[0]),&(Atree.Ktable[0]),
 //                       &(Atree.Kweights[0]),&(Atree.Cweights[0]),(double) Atree.homogen,&(Atree.cutoff),(int) Atree.n,&(Atree.dof[0]),&b_fmm[0],
 //                       (int) Atree.use_chebyshev,&(Atree.p_r2c[0]),&(Atree.p_c2r[0]));
	//  diff_fixedpoint_fmm = b_ - Map<VectorXd>(&b_fmm[0],N);
	//  printf("test Y: %20.5e\n",diff_fixedpoint_fmm.norm()/b_.norm());



	// 	for (int i =0; i<N;i++){
	//     	b_fmm[i] = 0.0;
 //    	}
	// 	compute.FMMCompute(&(Dtree.tree),&point_fmm[0],&point_fmm[0],&x_fmm[0],&(Dtree.K[0]),&(Dtree.U[0]),&(Dtree.VT[0]),&(Dtree.Tkz[0]),&(Dtree.Ktable[0]),
 //                       &(Dtree.Kweights[0]),&(Dtree.Cweights[0]),(double) Dtree.homogen,&(Dtree.cutoff),(int) Dtree.n,&(Dtree.dof[0]),&b_fmm[0],
 //                       (int) Dtree.use_chebyshev,&(Dtree.p_r2c[0]),&(Atree.p_c2r[0]));
	//  diff_fixedpoint_fmm = b_ - Map<VectorXd>(&b_fmm[0],N);
	//  printf("test Z: %20.5e\n",diff_fixedpoint_fmm.norm()/b_.norm());


	// 	for (int k=0;k<10;k++){
	// 		// printf("point_fmm[k].x,point_fmm[k].y,point_fmm[k].z: %f %f %f \n", point_fmm[k].x,point_fmm[k].y,point_fmm[k].z);
	// 		// printf("Atree.K[k]: %f \n", Atree.K[k]);
	// 		// printf("Atree.U[k]: %f \n", Atree.U[k]);
	// 		// printf("Atree.VT[k]: %f \n", Atree.VT[k]);
	// 		// printf("Atree.Tkz[k]: %f \n", Atree.Tkz[k]);
	// 		printf("x_fmm[k]: %f \n", x_fmm[k]);

	// 	}


	 // cout << endl;
	 // printf("Atree.L, Dtree.L: %f %f \n",Atree.L, Dtree.L);
	 // printf("Atree.level, Dtree.level: %d %d \n",Atree.level, Dtree.level);
	 // printf("Atree.n, Dtree.n: %d %d \n",Atree.n, Dtree.n);
	 // printf("Atree.epsilon, Dtree.epsilon: %f %f \n",Atree.epsilon, Dtree.epsilon);
	 // printf("Atree.homogen, Dtree.homogen: %f %f \n",Atree.homogen, Dtree.homogen);
	 // printf("Atree.symmetry, Dtree.symmetry: %d %d \n",Atree.symmetry, Dtree.symmetry);
	 // printf("Atree.alpha, Dtree.alpha: %f %f \n",Atree.alpha, Dtree.alpha);
	 // printf("Atree.use_chebyshev, Dtree.use_chebyshev: %d %d \n",Atree.use_chebyshev, Dtree.use_chebyshev);
	 // cout << endl;





	 // delete []q;
	 // delete []test3_;

	delete []x_fmm;
  	delete []b_fmm;

  	// delete []qtest;
  	// delete []b_fmmtest;



}
