#ifndef _container_h
#define _container_h

template<class T> class Vector {
 public:
  int m; // Size
  T * p; // Pointer to data

 Vector() :
  // m(-1), p(NULL), memory_allocated(false) {
  m(0), p(NULL), memory_allocated(false) {

  }
  ;

 Vector(int m_) :
  m(m_) {
#ifndef NDEBUG
    //20200108    assert(m > 0);
    assert(m >= 0); // N42, p11
#endif
    p = new T[m];
    memory_allocated = true;
#ifndef NDEBUG
    assert(p != NULL);
#endif
  }

 Vector(int m_, T * p_) :
  m(m_), p(p_) {
#ifndef NDEBUG
    //20200108    assert(m > 0);
    assert(m >= 0);
    assert(p != NULL);
    memory_allocated = false;
#endif    
  }

  ~Vector() {
    if (memory_allocated) {
      delete[] p;
      p = NULL;
    }
  }

  T& operator()(int i) {
#ifndef NDEBUG
    assert(i >= 0);
    assert(i < m);
#endif
    return p[i];
  }

  void resize(int m_) {
    m = m_;
#ifndef NDEBUG
    //20200108    assert(m > 0);
    assert(m >= 0);
#endif    
    assert((p == NULL) || (p != NULL && memory_allocated));
    delete[] p;
    p = new T[m];
    memory_allocated = true;
#ifndef NDEBUG
    assert(p != NULL);
#endif    
  }

 private:
  bool memory_allocated;
};

template<class T>
inline void zero(Vector<T> & v) {
  for (int i = 0; i < v.m; ++i) {
    v.p[i] = T(0);
  }
}

/* We use a column major format like Matlab and LAPACK */
template<class T> class Matrix {
 public:
  int m; // Number of rows
  int n; // Number of columns
  int ld; // Leading dimension
  T * p; // Pointer to data


 Matrix() :
  m(-1), n(-1), ld(-1), p(NULL), memory_allocated(false) {
  }
  ;

 Matrix(int m_, int n_) :
  m(m_), n(n_), ld(m) {
#ifndef NDEBUG
    assert(m > 0);
    assert(n > 0);
    assert(ld > 0);
    // printf("m,n: %d %d \n", m,n);
#endif
    p = new T[m * n];
    memory_allocated = true;
#ifndef NDEBUG
    assert(p != NULL);
#endif
  }

 Matrix(int m_, int n_, T * p_) :
  m(m_), n(n_), ld(m_), p(p_) {
#ifndef NDEBUG
    assert(m > 0);
    assert(n > 0);
    assert(ld > 0);
    assert(p != NULL);
    memory_allocated = false;
#endif    
  }

  ~Matrix() {
    if (p != NULL && memory_allocated)
      delete[] p;
  }

  T& operator()(int i, int j) {
#ifndef NDEBUG
    assert(i >= 0);
    // printf("i,m: %d %d \n", i,m);
    assert(i < m);
    assert(j >= 0);
    assert(j < n);
#endif
    return p[j * ld + i];
  }

  void resize(int m_, int n_) {
    m = m_;
    n = n_;
    ld = m;
#ifndef NDEBUG
    assert(m > 0);
    assert(n > 0);
#endif    
    assert((p == NULL) || (p != NULL && memory_allocated));
    delete[] p;
    p = new T[m * n];
    memory_allocated = true;
#ifndef NDEBUG
    assert(p != NULL);
#endif    
  }

 private:
  bool memory_allocated;
};

typedef Vector<int> ivec;
typedef Vector<double> dvec;
typedef Matrix<int> imat;
typedef Matrix<double> dmat;

#endif
