#ifndef _IFMM_H
#define _IFMM_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <cstring>
#include <list>
#include <cassert>
#include <sys/time.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <vector>
#include <fstream>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "debugmacros.h"
#include "mathmacros.h"
#include "infomesg.h"
#include "timer.h"

#include "container.h"
#include "set_L2P_operator.h"
#include "set_L2L_operator.h"
#include "set_P2M_operator.h"
#include "set_M2M_operator.h"
#include "set_M2P_operator.h"
#include "set_P2L_operator.h"
#include "set_M2L_operator.h"
#include "set_P2P_operator.h"
#include "set_nDOF.h"
#include "set_RHS.h"
#include "set_xyz.h"
#include "set_points.h"
#include "check_rank.h"
#include "get_weights.h"

#include "compute_Tk.h"
// #include "init_M2L_entries.h"
// #include "init_M2L_block.h"
// #include "create_fmm_tree.h"

#include "transfer_M2L_to_P2P.h"
#include "transfer_M2M_to_P2M.h"
#include "transfer_L2L_to_L2P.h"
#include "transfer_RHS.h"

#include "ACA_FullyPivoted.h"
//20200108#include "/Users/ericdarve/Documents/Pieter/redsvd-0.2.0/src/RedSVD-h.hpp"
#include "RedSVD-h.hpp"
// #include "/home/pcoulier/Software/RedSVD/RedSVD-h.hpp"
// #include "../BBFMM3D_pieter/include/bbfmm3d.hpp"


// #include "gmres_fmm.h"


// class myKernel: public H2_3D_Tree {
// public:
//     myKernel(double L, int level, int n,  double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n, epsilon, use_chebyshev){};
//     virtual void setHomogen(string& kernelType,doft *dof) {
//         homogen = 0;
//         symmetry = 1;
//         kernelType = "myKernel";
//         dof->s = 1;
//         dof->f = 1;
//     }
//     virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos, double *K, doft *dof) {
//         vector3 diff;        
//         // Compute 1/r
//         diff.x = sourcepos.x - fieldpos.x;
//         diff.y = sourcepos.y - fieldpos.y;
//         diff.z = sourcepos.z - fieldpos.z;
//         double r = sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
//         double a=1e-3;

//         if (r==0){
//           *K = 1;
//         }
//         else{
//           if(r<a){
//             *K = r/a;
//           }
//           else{
//             *K = a/r;
//           }
//         }
//     }
// };
 


// using std::list;
using namespace std;
using namespace Eigen;
using namespace REDSVD;


#define debugging_enabled 1

#define DEBUG(x) do { \
    if (debugging_enabled) { std::cerr << "[" << __FILE__ << ":" << __LINE__<< "] " << \
      x << std::endl; } \
} while (0)

typedef double timeType;

/* Uniform random number generator */
#define frand(xmin,xmax)\
    (xmin + (double)(xmax-xmin) *\
        double(rand()) / (1.0+double(RAND_MAX)))

//20200116/*
//20200116 * Function: Timer
//20200116 * ----------------------------------------
//20200116 * Returns the time in seconds.
//20200116 *
//20200116 */
//20200116timeType Timer(void);

class vec3 {
public:
  double x, y, z;
  vec3() {
  }

  vec3(double x_, double y_, double z_) :
    x(x_), y(y_), z(z_) {
  }

  vec3(double c) :
    x(c), y(c), z(c) {
  }

  inline double norm() {
    return sqrt(x * x + y * y + z * z);
  }
  inline double norm2() {
    return x * x + y * y + z * z;
  }
  vec3 operator-() {
    return vec3(-x, -y, -z);
  }
};

inline vec3 operator+(vec3 r1, vec3 r2) {
  return vec3(r1.x + r2.x, r1.y + r2.y, r1.z + r2.z);
}

inline vec3 operator-(vec3 r1, vec3 r2) {
  return vec3(r1.x - r2.x, r1.y - r2.y, r1.z - r2.z);
}

inline vec3 operator*(double a, vec3 r) {
  return vec3(a * r.x, a * r.y, a * r.z);
}

inline double dot(vec3 r0, vec3 r1) {
  return r0.x * r1.x + r0.y + r1.y + r0.z + r1.z;
}

inline std::ostream &operator <<(std::ostream &os, vec3 v) {
  return os << "(" << v.x << "," << v.y << "," << v.z << ")";
}

inline double norm_max(vec3 v) {
  return std::max(std::max(fabs(v.x), fabs(v.y)), fabs(v.z));
}

inline vec3 round(vec3 v) {
  return vec3(round(v.x), round(v.y), round(v.z));
}

inline vec3 abs(vec3 v) {
  return vec3(fabs(v.x), fabs(v.y), fabs(v.z));
}

/* BLAS and LAPACK routines on a Mac using XCode.
 See
 https://developer.apple.com/library/mac/navigation/#section=Topics&topic=Mathematical%20Computation
 https://developer.apple.com/library/mac/documentation/Accelerate/Reference/BLAS_Ref/Reference/reference.html#//apple_ref/doc/uid/TP40009457
 https://developer.apple.com/library/mac/documentation/Performance/Conceptual/vecLib/Reference/reference.html#//apple_ref/doc/uid/TP40002498
 */

// #include <Accelerate/Accelerate.h>

namespace blas {

// Declaration for BLAS matrix-matrix multiply
//void
//    dgemm(char *transa, char *transb, int *m, int *n, int *k, double *alpha,
//        double *A, int *lda, double *B, int *ldb, double *beta, double *C,
//        int *ldc);

  extern "C" {
    // Declaration for BLAS matrix-vector multiply
    void dgemv_(char *trans, int *m, int *n, double *alpha, double *A, int *lda,
	       double *x, int *incx, double *beta, double *y, int *incy);
  }
//// Declaration for BLAS dot product
//double ddot(int *n, double *dx, int *incx, double *dy, int *incy);

}

class diag_plus_lr {
public:
  virtual double operator()(vec3 &, vec3 &) = 0;
};

class diag_plus_lr_test: public diag_plus_lr {
public:
  double operator()(vec3 & source, vec3 & field) {
    vec3 r =source - field;
    double a0=0;
    double a1=1;
    if (r.norm()==0){
      return 1.0 + (a0+a1*source.x)*(a0+a1*field.x);
      // return 1.0;
    }
    else{
      // double test= field.y*source.z;
      // if (test !=0){
        // cout << "test:" << field.y*source.z << endl;
        // printf("source.x,source.y,source.z: %f %f %f \n",source.x,source.y,source.z);
        // printf("field.x,field.y,field.z: %f %f %f \n",field.x,field.y,field.z);
      // }
      return (a0+a1*source.x)*(a0+a1*field.x);
      // return 0;
    }
  }
  // double a=0.0001;
  //   if (r.norm()==0)
  //     {return 1;}
  //   else{
  //     if (r.norm()<a){
  //       return r.norm()/a;
  //     }
  //   else
  //     {return a / r.norm();}
  //   }
  // }
};

class kernel {
public:
  virtual double operator()(vec3 &) = 0;
};

class laplacian: public kernel {
public:
  double operator()(vec3 & r) {
    return 1.0 / r.norm();
  }
};

class testkernel: public kernel {
public:
  double operator()(vec3 & r) {
    if (r.norm()==0)
      // {return 0;}
      {return 1;}
  else
      // {return 1.0 / (r.norm()*r.norm()*r.norm()*r.norm()*r.norm()*r.norm());}
      // {return 1.0 / (r.norm()*r.norm());}
      // {return 1.0 / (r.norm()*r.norm()*r.norm()*r.norm());}
      {return 1.0 / (r.norm());}
  }
};

class benchmark_ii: public kernel {
public:
  double operator()(vec3 & r) {
  double a=1e-3;

  if (r.norm()==0)
      {return 1.0;}
  else{
    if (r.norm()<a){
      return r.norm()/a;
      }
    else
      {return a / r.norm();}
    }
  }
};

class benchmark_ii_adjusted: public kernel {
public:
  double operator()(vec3 & r) {
    double a=1e-3;
    if (r.norm()<a)
      {return 1;}
    else
      {return a / r.norm();}
    }
};


class benchmark_i: public kernel {
public:
  double operator()(vec3 & r) {
    double a=1e-3;
    if (r.norm()==0)
      {return 1;}
  else{
    if (r.norm()<a)
      {return (r.norm()*(log10(r.norm()) - 1))/(a*(log10(a) - 1));}
    else
      {return log(r.norm())/log(a);}
    }
  }
};

class benchmark_iii: public kernel {
public:
  double operator()(vec3 & r) {
    double a=1e-3;
    if (r.norm()==0)
      {return 1;}
      // {return 0;}
  else{
    if (r.norm()<a)
      {return (r.norm()*r.norm()*r.norm()*(3*log10(r.norm()) - 1))/(a*a*a*(3*log10(a) - 1));}
    else
      {return (r.norm()*r.norm()*log(r.norm()))/(a*a*log(a));}
    }
  }
};

class testkernel_identity: public kernel {
public:
  double operator()(vec3 & r) {
    if (r.norm()==0)
      {return 1;}
  else
      {return 0;}
  }
};

class laplacian_force: public kernel {
public:
  double operator()(vec3 & r) {
    double rn = r.norm();
    return r.x / (rn * rn * rn);
  }
};

class one_over_rfour: public kernel {
public:
  double operator()(vec3 & r) {
    double rn2 = r.norm2();
    return 1.0 / (rn2 * rn2);
  }
};

class poly0: public kernel {
public:
  double operator()(vec3 & r) {
    (void) (sizeof(r));
    /* remove annoying compiler warnings about unused variable. */
    return 1.0;
  }
};

class poly1: public kernel {
public:
  double operator()(vec3 & r) {
    return r.x;
  }
};

class poly2: public kernel {
public:
  double operator()(vec3 & r) {
    return r.x * r.x;
  }
};

class poly_diag: public kernel {
public:
  double operator()(vec3 & r) {
    double d=1e-3;
    if (r.norm()==0)
      {return d+r.x * r.y;}
      // {return d+pow(r.x, 3) * pow(r.y, 2) * r.z;}
      // {return d+(1+r.x +r.y * r.z);}
      // {return d+(1+r.x +r.y + r.z);}
    else{
      return r.x * r.y;
      // return pow(r.x, 3) * pow(r.y, 2) * r.z;
      // return (1+r.x +r.y + r.z);
    }
}
};


class poly3: public kernel {
public:
  double operator()(vec3 & r) {
    return pow(r.x, 3) * pow(r.y, 2) * r.z;
  }
};

class poly4: public kernel {
public:
  double operator()(vec3 & r) {
    return pow(r.x, 4) * pow(r.y, 4) * pow(r.z, 3);
  }
};

struct oct_node {
  int id;
  oct_node *parent;
  oct_node *child[8];
  Vector<oct_node *> ngbr; // List of neighbor nodes
  Vector<oct_node *> Ilist; // List of interaction nodes
  Vector<double *> K; // Pointer to kernel data
  list<int> pidx; // Global index of point
  Vector<vec3> pnt; // Vector of point coordinates
  dvec sigma;
  dvec x_leaf_check;
  vec3 ctr; // center
  double S; // size
  dvec M, L; // M and L data
  
  int nDOF;
  int rank;
  ivec rank_children_cum;


  MatrixXd P2M_operator;  // P2M (M2M at non-leaf levels)
  vector<MatrixXd> M2L_operator; // MatrixXd M2L_operator;  // M2L
  MatrixXd L2P_operator;  // L2P (L2L at non-leaf levels)
  vector<MatrixXd> P2P_operator; // MatrixXd P2P_operator;  // P2P (based on M2L operators of its children at non-leaf levels)
  vector<MatrixXd> P2L_operator; // MatrixXd P2L_operator;  // P2L
  vector<MatrixXd> M2P_operator; // MatrixXd M2P_operator;  // M2P
  vector<MatrixXd> M2M_operator; // M2M
  vector<MatrixXd> L2L_operator; // L2L


  MatrixXd Sigma_L2P;
  MatrixXd Sigma_P2M;

  VectorXd RHS;
  VectorXd RHS_leaf;

  PartialPivLU<MatrixXd> P2P_operator_self_LU;
  PartialPivLU<MatrixXd> alpha_self_LU;

  VectorXd x;
  VectorXd y;
  VectorXd z;


  int IsEliminated;

  int dim_cum_index_x;
  int dim_cum_index_y;
  int dim_cum_index_z;

#if defined(PARA_CHECK_RANK)
  //20200108  MatrixXcd R_P2M_operator;
  //20200108  MatrixXcd R_L2P_operator;
  MatrixXd R_P2M_operator;
  MatrixXd R_L2P_operator;
#endif

  Vector<oct_node *> concurrent; // for IFMM_PARALLELIZE

  oct_node() {
    parent = NULL;
    for (int i = 0; i < 8; ++i)
      child[i] = NULL;
  }

  ~oct_node() {
    for (int i = 0; i < 8; ++i)
      delete child[i];
  }
};

struct IFMM_Matrix {
private:
  int l; // Number of levels in the tree
  int n; // Order of Chebyshev polynomials
  double epsilon_rel_fillin;
  double epsilon_rel_basis;
  int LR_mode;
  bool Delay_update;
  bool Cubic_lattice;
  oct_node root; // Root node of tree
  dvec Sr; // Sr operator, stored as a vector
  dmat Tkmat; // Tk pre-computed data
  bool reuse;
  PartialPivLU<MatrixXd> DenseMatrix_top_level_LU;

  // void M2M();
  // void M2L();
  // void L2L(dvec & phi, Stokes_mobility * mfun);

public:
  ivec lvl_idx; // Starting index for nodes at each level
  Vector<oct_node *> lvl_p; // Vector of node pointers

  IFMM_Matrix(int l_, int n_, double epsilon_rel_fillin_, double epsilon_rel_basis_, int LR_mode_, bool Delay_update_): 
    l(l_), n(n_), epsilon_rel_fillin(epsilon_rel_fillin_), epsilon_rel_basis(epsilon_rel_basis_), LR_mode(LR_mode_), Delay_update(Delay_update_), Sr(n_ * n_), Tkmat(n_, n_) {
  }

  void init_tree_Ilist(vec3 & ctr, double S, Vector<vec3> & pnt);
  // void init_tree_Ilist_spheres(vec3 & ctr, double S, MatrixXd & sphereCenters, int & numPointsPerSphere, double & diam_sphere, Vector<vec3> & pnt);
  void create_lvl_indexing();
  // void setup_leaf_data(Vector<vec3> & point, dvec & b);
  void setup_leaf_data(Vector<vec3> & point);
  void setRHS(dvec & b);
  void setup_x_leaf_check(Vector<vec3> & point, dvec & x_leaf_check);
  void init_M2L_op(kernel * kfun, ivec & Kidx, double S, dvec & M2LOp);
  void init_M2L_p(ivec & Kidx, dvec & M2LOp);
  void compute(kernel * kfun, dvec & phi);

  
  // void initialization(Vector<vec3> & point, kernel * kfun);
  void initialization(Vector<vec3> & point, kernel * kfun,ofstream& outfile);
  // void initialize_ifmm_operators(Stokes_mobility * mfun,ofstream& outfile,double &DX);
  // void initialize_ifmm_operators(kernel * kfun,ofstream& outfile);
  void initialize_ifmm_operators(kernel * kfun,ofstream& outfile);
  // void elimination(double & sigma0_Aextended, ofstream& outfile);
  void elimination(ofstream& outfile);
  void eliminate_level(double & sigma0_Aextended, int curr_level, bool IsLeaf, ofstream& outfile);
  // void eliminate_level(double & sigma0_Aextended, int curr_level, bool IsLeaf);
  void solve_top_level();
  void downward_pass(dvec & x_iFMM,ofstream& outfile);
  // void downward_pass(dvec & x_iFMM);
  void substitution(dvec & x_iFMM,ofstream& outfile);
  // void substitution(dvec & x_iFMM);


  // void elimination_reuse(ofstream& outfile);
  void elimination_reuse();
  // void eliminate_level_reuse(int curr_level, bool IsLeaf, ofstream& outfile);
  void eliminate_level_reuse(int curr_level, bool IsLeaf);


  void estimate_lsv(double & sigma0_Aextended, int & leafLevel,ofstream& outfile);
  void A_extended_multiplication(MatrixXd &X, MatrixXd &AX, ivec &map_index_to_id, ivec &dimensions_cum, int & leafLevel, ofstream& outfile);
  void A_extended_transpose_multiplication(MatrixXd &X, MatrixXd &AX, ivec &map_index_to_id, ivec &dimensions_cum, int & leafLevel);

  // Following members are for parallelization

  void eliminate_level_node(double & sigma0_Aextended, int curr_level, ofstream& outfile, const double machine_eps, const int rank_RSVD_min, MatrixXd &P, const int k);
#if defined(IFMM_PARALLELIZE)
  void seek2_concurrent_nodes(const int curr_level, const int m);
  void statistics_concurrent_nodes();
#endif

  void downward_pass_node(dvec &x_iFMM, ofstream &outfile, const int curr_level, const int k);
  void eliminate_level_node_reuse(const int k);

};

// void set_potential(double L, dvec & b, Vector<vec3> & point);

// void compute_Tk(dvec & nodes, dmat & T);

void compute_Sr(Vector<vec3> & point, dmat & Sr_x, dmat & Sr_y, dmat & Sr_z,
    dmat & Tkmat);

void init_M2L_block(double b_size, dvec & nodes, vec3 & c_f, vec3 & c_s,dmat & Kmat, kernel * kfun);
// void init_M2L_block(double b_size, dvec & nodes, vec3 & c_f, vec3 & c_s,dmat & Kmat, diag_plus_lr * dlr_fun);


// void P2P(vec3 *fieldpos, vec3 *sourcepos, double *q, int Nf, int Ns, double *fieldval, kernel * kfun);
// void P2P(vec3 *fieldpos, vec3 *sourcepos, double *q, int Nf, int Ns, double *fieldval, diag_plus_lr * dlr_fun);

// void P2M(dvec & w, dvec & sigma, dmat & Sr_x, dmat & Sr_y, dmat & Sr_z);

// void L2P(dvec & phi, dvec & f, dmat & Sr_x, dmat & Sr_y, dmat & Sr_z);

void assemble_full(vec3 *location,int N, kernel * kfun,MatrixXd& A_full);
// void assemble_full(vec3 *location,int N, diag_plus_lr * dlr_fun,MatrixXd& A_full);

void ifmm(int l, int n, double epsilon, double epsilon_pre_SVD, int LR_mode, Vector<vec3> & point, dvec & q, kernel * kfun, dvec & phi, bool & Delay_update);
// void ifmm(int l, int n, double epsilon, int LR_mode, Vector<vec3> & point, dvec & q, diag_plus_lr * dlr_fun, dvec & phi, dvec & xdirect_check);


// void set_nDOF(oct_node * node, int &n, bool &IsLeaf);

// void set_xyz(oct_node * node);

// void set_P2M_operator(MatrixXd & P2M_operator,dmat & Sr_x, dmat & Sr_y, dmat & Sr_z);

// void set_M2M_operator(oct_node * node, MatrixXd & P2M_operator, int &n, dvec & Sr);

// void set_L2P_operator(MatrixXd & L2P_operator,dmat & Sr_x, dmat & Sr_y, dmat & Sr_z);

// void set_L2L_operator(oct_node * node, MatrixXd & L2P_operator, int &n,dvec & Sr);

// void set_M2L_operator(oct_node * node, int &n);

// void set_P2L_operator(oct_node * node, int &n, bool &IsLeaf);

// void set_M2P_operator(oct_node * node, int &n, bool &IsLeaf);

// void set_RHS(oct_node * node, bool &IsLeaf);

// void set_P2P_operator(oct_node * node,kernel * kfun, bool &IsLeaf);

template<class T>
void copy_list_to_vector(list<T> & l, Vector<T> & vec) {
  if (l.size() > 0) {
    vec.resize(l.size());
    typename list<T>::iterator it = l.begin();
    int i = 0;
    for (; it != l.end(); ++it, ++i)
      vec(i) = *it;
  }
}

double compute_error(double *phi, double *phidir, int Nf, int * chksum);

#endif
