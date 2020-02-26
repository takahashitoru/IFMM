#ifndef _FMM_H
#define _FMM_H

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


#include "container.h"

// using std::list;
using namespace std;
using namespace Eigen;

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

/*
 * Function: Timer
 * ----------------------------------------
 * Returns the time in seconds.
 *
 */
timeType Timer(void);

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
      {return 0;}
  else
      {return 1.0 / r.norm();}
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
  vec3 ctr; // center
  double S; // size
  dvec M, L; // M and L data

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

void set_potential(double L, dvec & b, Vector<vec3> & point);

void compute_Tk(dvec & nodes, dmat & T);
void compute_Sr(Vector<vec3> & point, dmat & Sr_x, dmat & Sr_y, dmat & Sr_z,
    dmat & Tkmat);

void init_M2L_block(double b_size, dvec & nodes, vec3 & c_f, vec3 & c_s,
    dmat & Kmat, kernel * kfun);

void P2P(vec3 *fieldpos, vec3 *sourcepos, double *q, int Nf, int Ns,
    double *fieldval, kernel * kfun);

void P2M(dvec & w, dvec & sigma, dmat & Sr_x, dmat & Sr_y, dmat & Sr_z);

void L2P(dvec & phi, dvec & f, dmat & Sr_x, dmat & Sr_y, dmat & Sr_z);

void assemble_full(vec3 *location,int N,kernel * kfun,MatrixXd& A_full);

void fmm(int l, int n, Vector<vec3> & point, dvec & q, kernel * kfun,
    dvec & phi);

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
