#include "gmres_fmm.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;
using namespace Eigen;



void gmres_fmm(H2_3D_Compute<myKernel> &compute, myKernel &Atree, vector3 * point_fmm, int &N, dvec &x, dvec &x0, dvec & b, int & max_iter, double & tol, IFMM_Matrix & ifmm, bool &precond, ofstream& outfile){

  
  VectorXd b_ = Map<VectorXd>(&b(0),N);

  MatrixXd v=MatrixXd::Zero(N,max_iter);
  MatrixXd H=MatrixXd::Zero(max_iter+1,max_iter);
  MatrixXd z=MatrixXd::Zero(N,max_iter);


  // USE FMM FOR MATRIX-VECTOR PRODUCT
  double *x_fmm = new double[N]; // Source array
  double *b_fmm = new double[N]; // Field array (BBFMM calculation)

  for (int k =0; k<N;k++){
    x_fmm[k] = x0(k);
  }
  // H2_3D_Compute<myKernel> compute(&Atree, point_fmm, point_fmm, N, N, x_fmm, m, b_fmm);
  // CALCULATE 'AX_0'
  compute.FMMCompute(&(Atree.tree),&point_fmm[0],&point_fmm[0],&x_fmm[0],&(Atree.K[0]),&(Atree.U[0]),&(Atree.VT[0]),&(Atree.Tkz[0]),&(Atree.Ktable[0]),
                     &(Atree.Kweights[0]),&(Atree.Cweights[0]),(double) Atree.homogen,&(Atree.cutoff),(int) Atree.n, &(Atree.dof[0]), &b_fmm[0],
                     (int) Atree.use_chebyshev,&(Atree.p_r2c[0]),&(Atree.p_c2r[0]));

  VectorXd rn = b_ - Map<VectorXd>(&b_fmm[0],N); // RESIDUE
  double beta = rn.norm(); // NORMALIZE
  v.col(0)=rn/beta;

  int i = 0;
  VectorXd err=VectorXd::Zero(max_iter);
  err(i)=beta/b_.norm();
  printf("Initial error:%20.5e\n",err(i));

  if (err(i)>1){
    for (int k =0; k<N;k++)
      x0(k)=0;

    for (int k =0; k<N;k++)
      x_fmm[k] = x0(k);
    
    // CALCULATE 'AX_0'
    // compute.FMMCompute(&(Atree.tree),&point_fmm[0],&point_fmm[0],&x_fmm[0],&(Atree.K[0]),&(Atree.U[0]),&(Atree.VT[0]),&(Atree.Tkz[0]),&(Atree.Ktable[0]),
                       // &(Atree.Kweights[0]),&(Atree.Cweights[0]),(double) Atree.homogen,&(Atree.cutoff),(int) Atree.n, &(Atree.dof[0]), &b_fmm[0],
                       // (int) Atree.use_chebyshev,&(Atree.p_r2c[0]),&(Atree.p_c2r[0]));

    for (int k =0; k<N;k++)
      b_fmm[k] = 0;

    rn = b_ - Map<VectorXd>(&b_fmm[0],N); // RESIDUE
    beta = rn.norm(); // NORMALIZE
    v.col(0)=rn/beta;

    err(i)=beta/b_.norm();
    printf("Initial error:%20.5e\n",err(i));

  }

  outfile << "Iteration: " << i << endl;
  outfile << "Residue: " << err(i) << endl;


  VectorXd z0=VectorXd::Zero(N);
  VectorXd zj=VectorXd::Zero(N);
  VectorXd vn=VectorXd::Zero(N);
  VectorXd w=VectorXd::Zero(N);

  MatrixXd Ri=MatrixXd::Zero(max_iter+2,max_iter+1);
  VectorXd gi=VectorXd::Zero(max_iter+2);

  dvec vn_(N);
  dvec zj_(N);


  while (i < max_iter-1 && err(i)>=tol){

    printf("GMRES iteration: %d\n",i);

    vn = v.col(i);

    if (precond==false){    
      // NO PRECONDITIONING
      zj=vn;
    }
    else{
      // IFMM PRECONDITIONING
      zero(zj_);
      for (int k=0;k<N;k++)
        vn_(k) = vn(k);

      ifmm.setRHS(vn_); //RESET THE RIGHT-HAND SIDE
      ifmm.elimination_reuse(); // RE-USE THE FACTORIZATION 
      ifmm.substitution(zj_,outfile); // PERFORM THE SUBSTITUTION

      zj=Map<VectorXd>(&zj_(0),N);
    }

    z.col(i)=zj;

    // USE FMM FOR MATRIX-VECTOR PRODUCT
    for (int k=0; k<N;k++)
      x_fmm[k]=zj(k);

    compute.FMMCompute(&(Atree.tree),&point_fmm[0],&point_fmm[0],&x_fmm[0],&(Atree.K[0]),&(Atree.U[0]),&(Atree.VT[0]),&(Atree.Tkz[0]),&(Atree.Ktable[0]),
                       &(Atree.Kweights[0]),&(Atree.Cweights[0]),(double) Atree.homogen,&(Atree.cutoff),(int) Atree.n,&(Atree.dof[0]),&b_fmm[0],
                       (int) Atree.use_chebyshev,&(Atree.p_r2c[0]),&(Atree.p_c2r[0]));
    w = Map<VectorXd>(&b_fmm[0],N); 


    for (int t=0; t<=i ; t++){
      H(t,i) = v.col(t).adjoint()*w;
      w.noalias() += - H(t,i)*v.col(t);
    }
    H(i+1,i)=w.norm();
    v.col(i+1) = w/H(i+1,i);

    Ri.conservativeResize(i+2,i+1);
    gi.conservativeResize(i+2);

    Ri=MatrixXd::Zero(i+2,i+1);
    gi=VectorXd::Zero(i+2);

    det_approx_sol(beta,i+1,H,Ri,gi);
    err(i+1)=abs(gi(i+1))/b_.norm();
    printf("Error: %20.5e\n",err(i+1));
    outfile << "Iteration: " << i+1 << endl;
    outfile << "Residue: " << err(i+1) << endl;

    i++;



  }
  cout << "End iteration" << endl;
  cout << "Total number of iterations: " << i << endl;
  outfile << "Total_number_of_iterations: " << i << endl;

  // GET 'X'
  if (i>0){
    VectorXd ym = Ri.block(0,0,i,i).partialPivLu().solve(gi.block(0,0,i,1));
    MatrixXd x_temp = z.block(0,0,N,i)*ym;
    x_temp += Map<VectorXd>(&x0(0),N);
    for (int k=0;k<N;k++){
      x(k) = x_temp(k);
    }
  }
  else{
    x=x0; // IF NO ITERATIONS ARE PERFORMED: X=X0
  }



  delete []x_fmm;
  delete []b_fmm;
}


// AUXILARY FUNCTION
void det_approx_sol(double &beta,int j,MatrixXd &H,MatrixXd &Ri, VectorXd &gi){

  Ri = H.block(0,0,j+1,j);
  gi(0)=beta;

  double hii, hi2;
  double temp;
  double c, s;

  for (int i=0; i<j ; i++){
    hii=Ri(i,i);
    hi2=Ri(i+1,i);
    
    if (abs(hi2)>abs(hii)){
      temp = hii/hi2; 
      s=1./pow(1.+pow(abs(temp),2),0.50);
      c = -temp*s;
    }
    else{
       temp = hi2/hii;
       c=1./pow(1.+pow(abs(temp),2),0.50);
       // c = 1/sqrt(1+abs(temp)^2);
       s = -temp*c;
    }
    MatrixXd W=MatrixXd::Identity(j+1,j+1); // MAKE THIS SPARSE
    // W = speye(j+1,j+1);
    W(i,i) = c;
    W(i,i+1) = -s;
    W(i+1,i) = s;
    W(i+1,i+1) = c;
    // W(i:i+1,i:i+1) = [[c',-s'];[s,c]];

    Ri = W*Ri;
    gi = W*gi;

  }


}