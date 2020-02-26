/*
  This software "LapackSVD.hpp" provides a wrapper of LAPACK's
  singular value decomposition (SVD) routines, i.e. {d,z}gesvd, in the
  fashion of Eigen's SVD routine, i.e. JacobiSVD.

  LapackSVD is a class to perform the SVD of a given matrix A with the
  help of Eigen and LAPACK+BLAS libraries. A real or complex matrix is
  decomposed as follows:

  1) m-by-n real matrix A is decomposed to U*S*V^t, where ^t stands for transpose,
  2) m-by-n complex matrix A is decomposed to U*S*V^*, where ^* stands for adjoint (conjugate transpose).

  Here is an example in the case of real:
  - - - - - - - - - - - - - - - - - - 

  #include "lapacksvd.hpp"

  int m = 3, n = 4;
  MatrixXd A = MatrixXd::Random(m, n); // give a real matrix as, for example, a random matrix

  LapackSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV); // second argument can be omitted; then 'ComputeFullU | ComputeFULLV' is assumed

  MatrixXd U = svd.matrixU(); // get the left singular matrix U
  MatrixXd V = svd.matrixV(); // get the right singular matrix V (not V^t)
  VectorXd s = svd.singularValues(); // a real vector of length MIN(m,n)

  MatrixXd P = U * s.asDiagonal() * V.transpose(); // reconstruct the matrix
  MatrixXd D = A - P; // difference from the reconstructed matrix P from the original matrix A
  cout << D.norm() / A.norm() << endl; // confirm that this relative norm is sufficiently small

  - - - - - - - - - - - - - - - - - - 
  
  To validate LapackSVD, you may compare it with Eigen's JacobiSVD. To
  this end, you may replace "LapackSVD" with "JacobiSVD".

  Note 1. Singular values are uniquely determined for every matrix. On
  the other hand, neither U nor V is uniquely determined.  Therefore,
  U and V of LapackSVD can be different from those of JacobiSVD.

  Note 2. Eigen includes a wrapper of Intel(TM) Math Kernel Library
  Lapack. See Eigin/src/SVD/JacobiSVD_MKL.h.

  Author: Toru Takahashi at Nagoya University, Japan

*/

#ifndef LAPACKSVD_HPP
#define LAPACKSVD_HPP

#include <complex>
#include <typeinfo>

#include <Eigen/Core>
#include <Eigen/Dense>


using namespace std;
using namespace Eigen;

#ifdef LAPACKSVD_DISABLE_UNDERSCORING

#define DGESVD dgesvd
#define ZGESVD zgesvd
#define DGESDD dgesdd
#define ZGESDD zgesdd

#else // default

#define DGESVD dgesvd_
#define ZGESVD zgesvd_
#define DGESDD dgesdd_
#define ZGESDD zgesdd_

#endif


#if !defined(EIGEN_USE_MKL)

extern "C" {
  void DGESVD(const char *JOBU, const char *JOBVT, const int *M, const int *N, double *A, const int *LDA, double *S, double *U, const int *LDU, double *VT, const int *LDVT, double *WORK, const int *LWORK, int *INFO);
  void ZGESVD(const char *JOBU, const char *JOBVT, const int *M, const int *N, complex<double> *A, const int *LDA, double *S, complex<double> *U, const int *LDU, complex<double> *VT, const int *LDVT, complex<double> *WORK, const int *LWORK, double *RWORK, int *INFO);
  void DGESDD(const char *JOBZ, const int *M, const int *N, double *A, const int *LDA, double *S, double *U, const int *LDU, double *VT, const int *LDVT, double *WORK, const int *LWORK, int *IWORK, int *INFO);
  void ZGESDD(const char *JOBZ, const int *M, const int *N, complex<double> *A, const int *LDA, double *S, complex<double> *U, const int *LDU, complex<double> *VT, const int *LDVT, complex<double> *WORK, const int *LWORK, double *RWORK, int *IWORK, int *INFO);
}

#endif

#include <Eigen/src/Core/util/Constants.h> // defines Compute{Full,Thin}{U,V}
#define FULL_U  ( ComputeFullU ) // Output U as m-by-m matrix
#define THIN_U  ( ComputeThinU ) // Output U as m-by-MIN(m,n) matrix
#define FULL_V  ( ComputeFullV ) // Output V as n-by-n matrix
#define THIN_V  ( ComputeThinV ) // Output V as n-by-MIN(m,n) matrix

#ifndef MAX
#define MAX(a, b) ( (a) > (b) ? (a) : (b) )
#endif
#ifndef MIN
#define MIN(a, b) ( (a) < (b) ? (a) : (b) )
#endif

template<typename T> class LapackSVD {

  typedef typename T::Scalar Scalar; // Scalar is double or complex<double>

private:

  int m, n; // number of rows and columns
  int min, max; // minimum and maximum of m and n
  int info; // return flag of {d,z}gesvd
  double *s; // singular values
  void *u, *vt; // U and V^t
  int Ucol, Vcol; // #columns of U and V
  
#ifdef LAPACKSVD_USE_XGESDD

  char jobz[1];
  int *iwork;

  void mydgesdd(void *A_, void **u, void **vt, double **s)
  {
    MatrixXd A = *((MatrixXd *)A_);
    
    *u = (void *)malloc(m * Ucol * sizeof(double)); // U is m-by-Ucol in column major
    *vt = (void *)malloc(Vcol * n * sizeof(double)); // V^t is Vcol-by-n in column major
    *s = new double [min]; // singular values (always real)
    
    int lwork = min * (6 + 4 * min) + max; // for JOBZ='A' or 'S'; not for 'N' or 'O'
    double *work = new double [lwork];
    int *iwork = new int [8 * min]; // dgesdd
    
    double *a = new double [m * n]; // A is m-by-n in column major and will be destroyed by dgesvd
    for (int i = 0; i < m; i ++) {
      for (int j = 0; j < n; j ++) {
	a[i + m * j] = A(i, j);
      }
    }
    
    DGESDD(jobz, &m, &n, a, &m, *s, (double *)*u, &m, (double *)*vt, &Vcol, work, &lwork, iwork, &info);
    
    delete [] a;
    delete [] work;
    delete [] iwork;
  }

  void myzgesdd(void *A_, void **u, void **vt, double **s)
  {
    MatrixXcd A = *((MatrixXcd *)A_);

    *u = (void *)malloc(m * Ucol * sizeof(complex<double>)); // U is m-by-Ucol in column major
    *vt = (void *)malloc(Vcol * n * sizeof(complex<double>)); // V^* is Vcol-by-n in column major
    *s = new double [min]; // singular values (always real)

    int lwork = 2 * min * min + 2 * min + max; // this is for JOBZ='O'; sufficient even for 'N', 'S', or 'A'
    complex<double> *work = new complex<double> [lwork];

    double *rwork = new double [min * MAX(5 * min + 7, 2 * max + 2 * min + 1)]; // sufficient even for JOBZ='N'

    int *iwork = new int [8 * min]; // zgesdd

    complex<double> *a = new complex<double> [m * n]; // A is m-by-n in column major and will be destroyed by zgesvd
    for (int i = 0; i < m; i ++) {
      for (int j = 0; j < n; j ++) {
	a[i + m * j] = A(i, j);
      }
    }

    ZGESDD(jobz, &m, &n, a, &m, *s, (complex<double> *)*u, &m, (complex<double> *)*vt, &Vcol, work, &lwork, rwork, iwork, &info);
    
    complex<double> *tmp = (complex<double> *)*vt;
    for (int i = 0; i < n * Vcol; i ++) {
      tmp[i] = conj(tmp[i]); // conjugate; then vt truely represents V^t (not V^*)
    }

    delete [] a;
    delete [] work;
    delete [] rwork;
    delete [] iwork;

  }

#else

  char jobu[1], jobvt[1];

  void mydgesvd(void *A_, void **u, void **vt, double **s)
  {
    MatrixXd A = *((MatrixXd *)A_);
    
    *u = (void *)malloc(m * Ucol * sizeof(double)); // U is m-by-Ucol in column major
    *vt = (void *)malloc(Vcol * n * sizeof(double)); // V^t is Vcol-by-n in column major
    *s = new double [min]; // singular values (always real)
    
    int lwork = 3 * min + max > 5 * min ? 3 * min + max : 5 * min;
    double *work = new double [lwork];
    
    double *a = new double [m * n]; // A is m-by-n in column major and will be destroyed by dgesvd
    for (int i = 0; i < m; i ++) {
      for (int j = 0; j < n; j ++) {
	a[i + m * j] = A(i, j);
      }
    }
    
    DGESVD(jobu, jobvt, &m, &n, a, &m, *s, (double *)*u, &m, (double *)*vt, &Vcol, work, &lwork, &info);
    
    delete [] a;
    delete [] work;
    
  }


  void myzgesvd(void *A_, void **u, void **vt, double **s)
  {
    MatrixXcd A = *((MatrixXcd *)A_);

    *u = (void *)malloc(m * Ucol * sizeof(complex<double>)); // U is m-by-Ucol in column major
    *vt = (void *)malloc(Vcol * n * sizeof(complex<double>)); // V^* is Vcol-by-n in column major
    *s = new double [min]; // singular values (always real)

    int lwork = 2 * min + max;
    complex<double> *work = new complex<double> [lwork];

    double *rwork = new double [5 * min];

    complex<double> *a = new complex<double> [m * n]; // A is m-by-n in column major and will be destroyed by zgesvd
    for (int i = 0; i < m; i ++) {
      for (int j = 0; j < n; j ++) {
	a[i + m * j] = A(i, j);
      }
    }

    ZGESVD(jobu, jobvt, &m, &n, a, &m, *s, (complex<double> *)*u, &m, (complex<double> *)*vt, &Vcol, work, &lwork, rwork, &info);
    
    complex<double> *tmp = (complex<double> *)*vt;
    for (int i = 0; i < n * Vcol; i ++) {
      tmp[i] = conj(tmp[i]); // conjugate; then vt truely represents V^t (not V^*)
    }

    delete [] a;
    delete [] work;
    delete [] rwork;

  }

#endif // ! LAPACKSVD_USE_XGESDD

  void compute(T &A, unsigned int options) {
    
    m = A.rows();
    n = A.cols();
    
    min = m < n ? m : n;
    max = m > n ? m : n;
    
#ifdef LAPACKSVD_USE_XGESDD

    if ((options & ComputeFullU) && (options & ComputeFullV)) {

      jobz[0] = 'a';
      Ucol = m;
      Vcol = n;

    } else if ((options & ComputeThinU) && (options & ComputeThinV)) {

      jobz[0] = 's';
      Ucol = min;
      Vcol = min;

    } else {

      cerr << "Invalid option. Exit." << endl;
      exit(EXIT_FAILURE);

    }

    if (typeid(Scalar) == typeid(double)) {

      mydgesdd((void *)&A, &u, &vt, &s);

    } else if (typeid(Scalar) == typeid(complex<double>)) {

      myzgesdd((void *)&A, &u, &vt, &s);

    } else {

      cerr << "Invalid type. Exit." << endl;
      exit(EXIT_FAILURE);

    }

#else // use {d,z}gesvd

    if ((options & ComputeFullU) != 0) {
      jobu[0] = 'a';
      Ucol = m;
    } else { // THIN_U
      jobu[0] = 's';
      Ucol = min;
    }

    if ((options & ComputeFullV) != 0) {
      jobvt[0] = 'a';
      Vcol = n;
    } else { // THIN_V
      jobvt[0] = 's';
      Vcol = min;
    }

    if (typeid(Scalar) == typeid(double)) {

      mydgesvd((void *)&A, &u, &vt, &s);

    } else if (typeid(Scalar) == typeid(complex<double>)) {

      myzgesvd((void *)&A, &u, &vt, &s);

    } else {

      cerr << "Invalid type. Exit." << endl;
      exit(EXIT_FAILURE);

    }

#endif
    
  }
  
public:
  
  LapackSVD(T &A) {
    compute(A, FULL_U | FULL_V);
  }
  
  LapackSVD(T &A, unsigned int options) {
    compute(A, options);
  }
  
  VectorXd singularValues() {
    VectorXd S = VectorXd::Zero(min);
    for (int i = 0; i < min; i ++) {
      S(i) = s[i];
    }
    return S;
  }
  
  T matrixU() {
    T U = T::Zero(m, Ucol);
    for (int i = 0; i < m; i ++) {
      for (int j = 0; j < Ucol; j ++) {
    	U(i, j) = ((Scalar *)u)[i + m * j];
      }
    }
    return U;
  }
  
  T matrixV() {
    T V = T::Zero(n, Vcol);
    for (int i = 0; i < n; i ++) {
      for (int j = 0; j < Vcol; j ++) {
	V(i, j) = ((Scalar *)vt)[j + Vcol * i]; // transpose; conjugation was already considered in the case of complex
      }
    }
    return V;
  }

  int getMax() {
    return max; // otherwise, let max be a public member
  }

  int getMin() {
    return min;
  }

  int getInfo() {
    return info;
  }

  ~LapackSVD() {
    delete [] s;
    //    delete [] u; // deleting ‘void*’ is undefined
    //    delete [] vt; // deleting ‘void*’ is undefined
    free(u);
    free(vt);
  }

};


#endif /* LAPACKSVD_HPP */
