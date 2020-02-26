#include "ifmm.h"
// #include "gmres_fmm.h"
// #include "iteration_fixedpoint.h"
// #include "../BBFMM3D_pieter/include/bbfmm3d.hpp"



int main(int argc, char *argv[]) {
  // Solve 'Ax=b' (approximately) with the inverse fast multipole method
  // 
  // Pieter Coulier
  // 2015
  //
  // A parallelized version
  // Toru Takahashi
  // 2020

  int l = atoi(argv[1]); // NUMBER OF LEVELS IN THE FMM-TREE (MUST BE LARGER THAN 2)
  int n = atoi(argv[2]); // NUMBER OF CHEBYSHEV POINTS (IN EACH DIRECTION) FOR THE LOW-RANK APPROXIMATIONS
  int N = atoi(argv[3]); // TOTAL NUMBER OF POINTS

  /* initialize random seed (!): */
  srand(284851);

  // DISTRIBUTE N POINTS RANDOMLY IN THE UNIT CUBE
  std::cerr << "Distribute points..." << std::endl;
  double L = 1.0; // SIZE OF THE CUBE (3D)
  Vector<vec3> point(N);
  set_points(L, point);

  // KERNEL TYPE
  benchmark_ii fmm_kernel;

  // IFMM PARAMETERS
  double epsilon_rel_fillin = 1e-5; // RELATIVE ACCURACY FOR COMPRESSING MATRIX FILL-INS
  double epsilon_rel_basis = 1e-5; // RELATIVE ACCURACY FOR UPDATING THE BASIS OF EACH NODE
  int LR_mode=3; // RECOMPRESSION TECHNIQUE FOR THE FILL-INS (LR_MODE=0:SVD, LR_MODE=1:ACA, LR_MODE=2:RSVD;  LR_MODE=3:RSVD - OPTIMIZED)
  bool Delay_update = false;
  
  // ********************************* //
  // * INVERSE FAST MULTIPOLE METHOD * //
  // ********************************* //
  
  // CHOOSE THE VECTOR 'X_EXACT'(FOR VALIDATION PURPOSES ONLY) (RANDOM VECTOR)
  dvec x_exact(N);
  for (int k=0; k<N; k++){
    x_exact(k)=frand(-1, 1);
  }

  // CALCULATE 'B=AX'(FOR VALIDATION PURPOSES ONLY)
  dvec b(N); // RIGHT-HAND SIDE
  std::cerr << "Compute RHS..." << std::endl;
  zero(b);
  for(int i=0;i<N;++i){
    for(int j=0;j<N;++j){
      vec3 r = point(j) - point(i);
      b(i)+=(fmm_kernel)(r)*x_exact(j);
    }
  }

  ofstream outfile; // OUTPUT-FILE
  char filename [200];
  //20200116  sprintf(filename,"output/3D_randompoints_benchmark_ii_a_1e-3_paper_N_%d_n_%d_l_%d_eps_rel_fillin_%1.1e_eps_rel_basis_%1.1e.txt",N,n,l,epsilon_rel_fillin,epsilon_rel_basis);
  fprintf(stderr, "Dig a subdirectory 'output' beforehand!!\n");
  sprintf(filename,"output/%s_3D_randompoints_benchmark_ii_a_1e-3_paper_N_%d_n_%d_l_%d_eps_rel_fillin_%1.1e_eps_rel_basis_%1.1e.txt",argv[0],N,n,l,epsilon_rel_fillin,epsilon_rel_basis);
  outfile.open(filename);

  outfile << "N: " << N << endl;
  outfile << "n: " << n << endl;
  outfile << "l: " << l << endl;
  outfile << "epsilon_rel_fillin: " << epsilon_rel_fillin << endl;
  outfile << "epsilon_rel_basis: " << epsilon_rel_basis << endl;
  outfile << endl;

  TICK(all);

  // NOW SOLVE FOR 'X_IFMM'

  dvec x_iFMM(N); // VECTOR OF UNKNOWNS THAT WE'RE SOLVING FOR (x_iFMM = A\b)
  zero(x_iFMM);


  std::cerr << "ifmm_intialization..." << std::endl;
  IFMM_Matrix ifmm(l, n, epsilon_rel_fillin, epsilon_rel_basis, LR_mode, Delay_update); // CREATE AN OBJECT OF 'IFMM_MATRIX'
  TICK(ifmm_initialization);
  ifmm.initialization(point,&fmm_kernel,outfile); // INITIALIZE ALL THE NECESSARY OPERATORS
  TACK(ifmm_initialization, std::cerr);

  std::cerr << "ifmm_setRHS..." << std::endl;
  TICK(ifmm_setRHS);
  ifmm.setRHS(b); // SET THE RIGHT-HAND SIDE
  TACK(ifmm_setRHS, std::cerr);  

  std::cerr << "ifmm_elimination..." << std::endl;
  TICK(ifmm_elimination);
  ifmm.elimination(outfile); // PERFORM THE ELIMINATION/FACTORIZATION
  TACK(ifmm_elimination, std::cerr);

  std::cerr << "ifmm_substitution..." << std::endl;
  TICK(ifmm_substitution);
  ifmm.substitution(x_iFMM,outfile); // PERFORM THE SUBSTITUTION
  TACK(ifmm_substitution, std::cerr);

  // CALCULATE RELATIVE ERROR (FOR VALIDATION PURPOSES ONLY)
  std::cerr << "Calculate relative error..." << std::endl;
  VectorXd x_exact_ = Map<VectorXd>(&x_exact(0),N);
  VectorXd x_iFMM_ = Map<VectorXd>(&x_iFMM(0),N);
  
  double relative_error = (x_iFMM_ - x_exact_).norm()/x_exact_.norm();
  std::cerr << "Relative error iFMM:" << relative_error << std::endl;
  outfile << "Relative_error_iFMM: " << relative_error << endl;


  // RE-USE THE FACTORIZATION TO SOLVE THE SYSTEM SEQUENTIALLY FOR 'N_RHS' OTHER RIGHT-HAND SIDES  (x_2 = A\b_2, x_3 = A\b_3, ...)

  TCREATE(ifmm_setRHS_x_N_RHS);
  TCREATE(ifmm_elimination_reuse_x_N_RHS);
  TCREATE(ifmm_substitution_x_N_RHS);

  int N_RHS=10;
  std::cerr << "N_RHS=" << N_RHS << std::endl;

  for (int k=0; k<N_RHS; k++){
    //    std::cerr << "New RHS IFMM k=" << k << " ..." << std::endl;

    for (int j=0; j<N; j++){ // NEW RIGHT-HAND SIDE (RANDOM VECTOR)
      b(j)=frand(-0.5, 0.5);
    }
    zero(x_iFMM);

    TSTART(ifmm_setRHS_x_N_RHS);
    ifmm.setRHS(b); //RESET THE RIGHT-HAND SIDE
    TSTOP(ifmm_setRHS_x_N_RHS);

    TSTART(ifmm_elimination_reuse_x_N_RHS);
    ifmm.elimination_reuse(); // RE-USE THE FACTORIZATION 
    TSTOP(ifmm_elimination_reuse_x_N_RHS);

    TSTART(ifmm_substitution_x_N_RHS);
    ifmm.substitution(x_iFMM,outfile); // PERFORM THE SUBSTITUTION
    TSTOP(ifmm_substitution_x_N_RHS);
  }

  TPRINT(ifmm_setRHS_x_N_RHS, std::cerr), TFREE(ifmm_setRHS_x_N_RHS);
  TPRINT(ifmm_elimination_reuse_x_N_RHS, std::cerr), TFREE(ifmm_elimination_reuse_x_N_RHS);
  TPRINT(ifmm_substitution_x_N_RHS, std::cerr), TFREE(ifmm_substitution_x_N_RHS);


  TACK(all, std::cerr);

  outfile.close();


  std::cerr << "Run was successful" << std::endl;
  return 0;
}
