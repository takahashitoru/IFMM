#include "fmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

//20200116#define evaltime(timeval_time) (double)timeval_time.tv_sec	\
//20200116    + (double)timeval_time.tv_usec*1e-6
//20200116
//20200116/*
//20200116 * Function: Timer
//20200116 * ----------------------------------------
//20200116 * Returns the time in seconds.
//20200116 *
//20200116 */
//20200116timeType Timer(void) {
//20200116  struct timeval timeval_time;
//20200116  gettimeofday(&timeval_time, NULL);
//20200116  return evaltime(timeval_time);
//20200116}

/*
// Declaration for BLAS matrix-vector multiply
void blas::dgemv(char *trans, int *m, int *n, double *alpha, double *A,
    int *lda, double *x, int *incx, double *beta, double *y, int *incy) {

  char letterN[] = "n";
  CBLAS_TRANSPOSE TransA;

  if (strcmp(trans, letterN) == 0) {
    TransA = CblasNoTrans;
  } else {
    TransA = CblasTrans;
  }

  cblas_dgemv(CblasColMajor, TransA, *m, *n, *alpha, A, *lda, x, *incx, *beta,
      y, *incy);
}
*/

void set_potential(double L, dvec & b, Vector<vec3> & point) {
  for (int i = 0; i < point.m; i++) {
    b(i) = frand(-1, 1);

    point(i).x = frand(-0.5, 0.5) * L;
    point(i).y = frand(-0.5, 0.5) * L;
    point(i).z = frand(-0.5, 0.5) * L;
  }
}

void compute_Tk(dvec & nodes, dmat & T) {
  const int r = nodes.m;
  for (int i = 0; i < r; ++i) {
    const double x = nodes(i);
    T(0, i) = 1;
    T(1, i) = x;
    for (int k = 2; k < r; k++)
      T(k, i) = 2.0 * x * T(k - 1, i) - T(k - 2, i);
  }
}

double clenshaw(int n, double x, double *Tk) {
  double d0, d1, d2;
  d0 = d1 = 0;
  for (int j = n - 1; j > 0; j--) {
    d2 = d0;
    d0 = 2.0 * x * d0 - d1 + Tk[j];
    d1 = d2;
  }
  return x * d0 - d1 + 0.5 * Tk[0];
}

void compute_Sr(Vector<vec3> & point, dmat & Sr_x, dmat & Sr_y, dmat & Sr_z,
    dmat & Tkmat) {
  const int n = Sr_x.m; // Number of rows
  const int r = Sr_x.n; // Number of columns
  double pfac = 2.0 / r;

  // Loop over Chebyshev node m
  for (int m = 0; m < r; m++) {
    // Compute S_n for each direction independently using Clenshaw
    double * Tk = &Tkmat(0, m);
    for (int i = 0; i < n; i++) {
      Sr_x(i, m) = pfac * clenshaw(r, point(i).x, Tk);
      Sr_y(i, m) = pfac * clenshaw(r, point(i).y, Tk);
      Sr_z(i, m) = pfac * clenshaw(r, point(i).z, Tk);
    }
  }
}

int BSD_checksum(int* data, int n) {
  int checksum = 0;
  for (int i = 0; i < n; ++i) {
    checksum = (checksum >> 1) + ((checksum & 1) << 15);
    checksum += data[i];
    checksum &= 0xffff; /* Keep it within bounds. */
  }

  return checksum;
}

double compute_error(double *phi, double *phidir, int N, int * checksum) {

  double sum1, sum2;

  dvec diff(N);
  for (int i = 0; i < N; i++) {
    diff(i) = phi[i] - phidir[i];
    //    DEBUG(phi[i] << " " << phidir[i]);
  }

  *checksum = BSD_checksum((int*) (diff.p), sizeof(double) / sizeof(int) * N);

  sum1 = 0;
  sum2 = 0;
  for (int i = 0; i < N; i++) {
    sum1 += diff(i) * diff(i);
    sum2 += phidir[i] * phidir[i];
  }

  return sqrt(sum1 / sum2);
}

struct fmm_data {
private:
  int l; // Number of levels in the tree
  int n; // Order of Chebyshev polynomials
  oct_node root; // Root node of tree
  dvec Sr; // Sr operator, stored as a vector
  dmat Tkmat; // Tk pre-computed data

  void M2M();
  void M2L();
  void L2L(dvec & phi, kernel * kfun);

public:
  ivec lvl_idx; // Starting index for nodes at each level
  Vector<oct_node *> lvl_p; // Vector of node pointers

  fmm_data(int l_, int n_) : 
    l(l_), n(n_), Sr(n_ * n_), Tkmat(n_, n_) {
  }

  void init_tree_Ilist(vec3 & ctr, double S, Vector<vec3> & pnt);
  void create_lvl_indexing();
  void setup_leaf_data(Vector<vec3> & point, dvec & q);
  void init_M2L_op(kernel * kfun, ivec & Kidx, double S, dvec & M2LOp);
  void init_M2L_p(ivec & Kidx, dvec & M2LOp);
  void compute(kernel * kfun, dvec & phi);
};

void nnode_at_level(bool save, int l, oct_node * node, ivec & nnode_level,
    Vector<oct_node *> & nptr) {
  for (int i = 0; i < 8; ++i)
    if (node->child[i] != NULL)
      nnode_at_level(save, l + 1, node->child[i], nnode_level, nptr);

  if (save) {
    node->id = nnode_level(l);
    // std::cout << "node->id: %d" << node->id << std::endl;
    nptr(node->id) = node;
  }
  ++nnode_level(l);

  // printf("Level, nnode_level: %d %d\n",l,nnode_level(l));

}

void cum_sum(ivec & p) {
  int tmp0, tmp1;
  tmp0 = tmp1 = 0;
  for (int i = 0; i < p.m - 1; ++i) {
    tmp1 += p(i);
    p(i) = tmp0;
    tmp0 = tmp1;
  }
  p(p.m - 1) = tmp1; // ???
}

inline int M2L_index_map(vec3 & i) {
  return int(49 * (i.x + 3) + 7 * (i.y + 3) + i.z + 3);
}

void init_X2X_op(dmat & Tkmat, dvec & Sr, int n) { // ???
  dvec nodes(n);
  const int Nc = 2 * n;
  Vector<vec3> fieldt(Nc); // Chebyshev-transformed coordinates
  double pi = M_PI;

  // Compute the n Chebyshev nodes of T_n(x)
  for (int m = 0; m < n; m++)
    nodes(m) = cos(pi * ((double) m + 0.5) / (double) n);

  // Evaluate the Chebyshev polynomials of degree 0 to n-1 at the nodes
  compute_Tk(nodes, Tkmat);

  // Map all Chebyshev nodes from the children cells to the parent
  double z = -1.;
  for (int l3 = 0; l3 < n; l3++) {
    fieldt(l3).x = 0.;
    fieldt(l3).y = 0.;
    fieldt(l3).z = 0.5 * (nodes(l3) + z);
  }

  // Compute Sr, the mapping function for the field points
  dmat Sr_x(Nc, n);
  dmat Sr_y(Nc, n);
  dmat Sr_z(Nc, n);
  compute_Sr(fieldt, Sr_x, Sr_y, Sr_z, Tkmat);

  // Extract out the Chebyshev weights
  for (int l1 = 0; l1 < n; l1++)
    for (int l2 = 0; l2 < n; l2++)
      Sr(l1 * n + l2) = Sr_z(l2, l1);
}

void init_M2L_entries(int n, int l, double S, ivec & Kidx, kernel * kfun,
    dvec & M2LOp) {
  int n3 = n * n * n; // n3 = n^3
  int n6 = n3 * n3;

  dvec nodes(n);

  // Compute the n Chebyshev nodes of T_n(x)
  double pi = M_PI;
  for (int m = 0; m < n; m++)
    nodes(m) = cos(pi * ((double) m + 0.5) / (double) n);

  // Compute the kernel values for interactions with all 316 cells at all levels
  S *= (1.0 / (1 << l));
  for (int lvl = l; lvl >= 2; --lvl) {
    for (int k1 = -3; k1 < 4; k1++) {
      vec3 fcenter(0, 0, 0);
      vec3 scenter;
      scenter.x = k1 * S;
      for (int k2 = -3; k2 < 4; k2++) {
        scenter.y = k2 * S;
        for (int k3 = -3; k3 < 4; k3++) {
          scenter.z = k3 * S;
          if (abs(k1) > 1 || abs(k2) > 1 || abs(k3) > 1) {
            // Compute the kernel at each of the Chebyshev nodes
            vec3 kidx(k1, k2, k3);
            int ncell = M2L_index_map(kidx);
            ncell = Kidx(ncell) * n6 + (l - lvl) * n6 * 316;
            dmat Kmat0(n3, n3, M2LOp.p + ncell);
            init_M2L_block(S, nodes, fcenter, scenter, Kmat0, kfun);
          }
        }
      }
    }
    S *= 2.0;
  }
}

void fmm_data::init_M2L_op(kernel * kfun, ivec & Kidx, double S, dvec & M2LOp) {

  // Compute the Chebyshev weights and sets up the lookup table
  init_X2X_op(Tkmat, Sr, n);

  // Initialize lookup table
  int ncell = 0;
  int ninteract = 0;
  for (int l1 = -3; l1 < 4; l1++) {
    for (int l2 = -3; l2 < 4; l2++) {
      for (int l3 = -3; l3 < 4; l3++, ncell++) {
        if (abs(l1) > 1 || abs(l2) > 1 || abs(l3) > 1) {
          Kidx(ncell) = ninteract;
          ninteract++;
        }
      }
    }
  }

  // for (int k=0; k<343;k++)
  // {
  //   printf("Kidx: %d \n",Kidx(k));
  // }


  // Precompute the M2L interaction matrices
  int n3 = n * n * n;
  int n6 = n3 * n3;
  if (l > 1)
    M2LOp.resize(n6 * 316 * (l - 1));
  init_M2L_entries(n, l, S, Kidx, kfun, M2LOp);
}

// Debug function to test that the tree was built correctly
int count_points(oct_node * node) {
  int n = 0;
  for (int i = 0; i < 8; ++i)
    if (node->child[i] != NULL)
      n += count_points(node->child[i]);

  n += node->pidx.size();

  return n;
}

void create_FMM_tree(int levels, int idx, vec3 & point, oct_node * node) {
  vec3 ctr = node->ctr;
  double S = node->S;

  for (int l_ = levels; l_ >= 1; --l_) {

    // Index of child cell containing point
    int icell = ((point.z < ctr.z) ? 0 : 1) + 2 * ((point.y < ctr.y) ? 0 : 1)
        + 4 * ((point.x < ctr.x) ? 0 : 1);

    // printf("Level, iCell: %d %d \n",l_,icell);

    S = 0.5 * S;

    ctr.z += 0.5 * S * ((icell & 1) ? 1 : -1);
    ctr.y += 0.5 * S * ((icell & 2) ? 1 : -1);
    ctr.x += 0.5 * S * ((icell & 4) ? 1 : -1);

    // Check whether child node exists or not
    if (node->child[icell] == NULL) {
      oct_node * child = (node->child[icell] = new oct_node);
      child->parent = node;
      child->ctr = ctr;
      child->S = S;
    }
    node = node->child[icell];
  }

  // This is a leaf node; add the current point
  
  node->pidx.push_back(idx);
}

bool is_well_separated(double S, vec3 & p) {
  vec3 r = round(abs((1.0 / S) * p));
  return r.x > 1 || r.y > 1 || r.z > 1;
}

void create_Ilist(oct_node * node) {
  double S = node->S;
  vec3 & ctr = node->ctr;

  if (node->parent != NULL) {
    list<oct_node*> Ilist_, ngbr_;
    // Loop over parent's neighbors
    Vector<oct_node*> & pn = node->parent->ngbr;
    for (int j = 0; j < pn.m; ++j) {
      // Loop over children clusters
      for (int i = 0; i < 8; ++i) {
        oct_node * child = pn(j)->child[i];
        if (child != NULL) {
          // Test to determine whether this node is a neighbor
          // or is well-separated
          vec3 r = child->ctr - ctr;
          if (is_well_separated(S, r))
            Ilist_.push_back(child);
          else
            ngbr_.push_back(child);
        }
      }
    }

    // Copy back into vectors
    copy_list_to_vector<oct_node*> (Ilist_, node->Ilist);
    copy_list_to_vector<oct_node*> (ngbr_, node->ngbr);
  }

  // for (int k=0;k<1;k++)
  // {
  //   printf("Ilist of node %d  \n",node->Ilist(0)->id);
  // }


  for (int i = 0; i < 8; ++i)
    if (node->child[i] != NULL)
      create_Ilist(node->child[i]);
}

void fmm_data::init_tree_Ilist(vec3 & ctr, double S, Vector<vec3> & pnt) {
  root.ctr = ctr;
  root.S = S;

  printf("root.ctr: %f %f %f\n",root.ctr.x,root.ctr.y,root.ctr.z);


  printf("point0: %f %f %f\n",pnt(0).x,pnt(0).y,pnt(0).z);



  for (int i = 0; i < pnt.m; ++i)
    create_FMM_tree(l, i, pnt(i), &root);

  root.ngbr.resize(1);
  root.ngbr(0) = &root;
  create_Ilist(&root);
}

void init_M2L_block(double b_size, dvec & nodes, vec3 & c_f, vec3 & c_s,
    dmat & Kmat, kernel * kfun) {
  const int n = nodes.m;
  b_size /= 2.;

  // Location of source points
  vec3 source;
  int i, j;
  j = 0;
  for (int l1 = 0; l1 < n; l1++) {
    source.x = c_s.x + b_size * nodes(l1);
    for (int l2 = 0; l2 < n; l2++) {
      source.y = c_s.y + b_size * nodes(l2);
      for (int l3 = 0; l3 < n; l3++, j++) {
        source.z = c_s.z + b_size * nodes(l3);

        // Location of field points
        vec3 field;
        i = 0;
        for (int k1 = 0; k1 < n; k1++) {
          field.x = c_f.x + b_size * nodes(k1);
          for (int k2 = 0; k2 < n; k2++) {
            field.y = c_f.y + b_size * nodes(k2);
            for (int k3 = 0; k3 < n; k3++, i++) {
              field.z = c_f.z + b_size * nodes(k3);
              vec3 r = source - field;
              Kmat(i, j) = (*kfun)(r);
            }
          }
        }
      }
    }
  }
}

void P2P(vec3 *field, vec3 *source, double *q, int Nf, int Ns,
    double *fieldval, kernel * kfun) {
  int i, j;

  for (i = 0; i < Nf; i++) {
    // Compute the interaction between each field point and source
    for (j = 0; j < Ns; j++) {
      if (source[j].x != field[i].x || source[j].y != field[i].y || source[j].z
          != field[i].z) {
        vec3 r = source[j] - field[i];
        fieldval[i] += q[j] * (*kfun)(r);
      }
    }
  }
}

void P2M(dvec & w, dvec & sigma, dmat & Sr_x, dmat & Sr_y, dmat & Sr_z) {
  const int nx = Sr_x.m;
  const int r = Sr_x.n;
  int k = 0;
  for (int l1 = 0; l1 < r; l1++) {
    for (int l2 = 0; l2 < r; l2++) {
      for (int l3 = 0; l3 < r; l3++, k++) {
        for (int j = 0; j < nx; j++) {
          w(k) += sigma(j) * Sr_x(j, l1) * Sr_y(j, l2) * Sr_z(j, l3);
        }
      }
    }
  }
}

void L2P(dvec & phi, dvec & f, dmat & Sr_x, dmat & Sr_y, dmat & Sr_z) {
  const int nx = Sr_x.m;
  const int r = Sr_x.n;
  for (int i = 0; i < nx; i++) {
    int k = 0;
    for (int l1 = 0; l1 < r; l1++) {
      for (int l2 = 0; l2 < r; l2++) {
        double tmp = Sr_x(i, l1) * Sr_y(i, l2);
        for (int l3 = 0; l3 < r; l3++, k++)
          phi(i) += f(k) * tmp * Sr_z(i, l3);
      }
    }
  }
}


void assemble_full(vec3 *location,int N,kernel * kfun,MatrixXd& A_full){
  A_full = MatrixXd::Zero(N,N);
  for(unsigned long i=0;i<N;++i){
    for(unsigned long j=0;j<N;++j){
        vec3 r = location[j] - location[i];
        A_full(i,j) =   (*kfun)(r);
        // printf("r - A(r): %f %f %f %f \n",r.x,r.y,r.z,(*kfun)(r));
    }
  }

}

void fmm_data::M2L() {
  int n3 = n * n * n;
  char trans[] = "n";
  double alpha = 1.;
  int incr = 1;
  for (int lvl = 2; lvl < lvl_idx.m - 1; ++lvl) {
    for (int k = lvl_idx(lvl); k < lvl_idx(lvl + 1); ++k) {

      // printf("M2l - lvl, k: %d %d\n",lvl,k);


      oct_node * on = lvl_p(k);
      zero(on->L);



      for (int i = 0; i < on->Ilist.m; ++i) {
      // printf("Ilist - lvl, k: %d \n",i);

        double * M = on->Ilist(i)->M.p;
        blas::dgemv_(trans, &n3, &n3, &alpha, on->K(i), &n3, M, &incr, &alpha,
            on->L.p, &incr);
      }
    }
  }
}

void init_Sr_cell(oct_node* on, dmat & Tkmat, dmat & Sr_x, dmat & Sr_y,
    dmat & Sr_z) {
  vec3 fcenter = on->ctr; // Center of field cell

  int Nf = on->pnt.m; // Number of points
  Vector<vec3> fieldt(Nf); // Chebyshev-transformed coordinates

  // Map all of the field points to the box ([-1 1])^3
  const double ihalfL = 2.0 / on->S;
  for (int i = 0; i < Nf; ++i)
    fieldt(i) = ihalfL * (on->pnt(i) - fcenter);

  // Compute Sf, the mapping function for the field points
  compute_Sr(fieldt, Sr_x, Sr_y, Sr_z, Tkmat);
}

void Sr_indexing(oct_node * on, int & xcount, int * xindex, int & ycount,
    int * yindex, int & zcount, int * zindex) {

  xcount = ycount = zcount = 0;

  // Determine which children cells contain field points
  for (int i = 0; i < 8; i++) {
    if (on->child[i] != NULL) {
      zindex[zcount] = i;
      zcount++;
      if (ycount == 0 || yindex[ycount - 1] != i / 2) {
        yindex[ycount] = i / 2;
        ycount++;
      }
      if (xcount == 0 || xindex[xcount - 1] != i / 4) {
        xindex[xcount] = i / 4;
        xcount++;
      }
    }
  }
}

void Sr_idx_get(int idx, int & j1, int & j2, int & incr) {
  j1 = idx;
  j2 = (int) (j1 / 2);
  (j1 % 2 == 0) ? incr = 1 : incr = -1;
}

void Sr_xgemv(char * trans, int n, int incr, double * Sr, double * pin,
    double * pout) {
  double alpha = 1.;
  double beta = 1.;
  const int n2 = n * n;

  for (int l2 = 0; l2 < n2; l2 += n) {
    for (int l3 = l2; l3 < l2 + n; l3++) {
      blas::dgemv_(trans, &n, &n, &alpha, Sr, &n, pin + l3, &incr, &beta, pout
          + l3, &incr);
    }
  }
}

void Sr_ygemv(char * trans, int n, int incr, double * Sr, double * pin,
    double * pout) {
  double alpha = 1.;
  double beta = 1.;
  const int n2 = n * n;
  const int n3 = n2 * n;

  for (int l1 = 0; l1 < n3; l1 += n2) {
    for (int l3 = l1; l3 < l1 + n; l3++)
      blas::dgemv_(trans, &n, &n, &alpha, Sr, &n, pin + l3, &incr, &beta, pout
          + l3, &incr);
  }
}

void Sr_zgemv(char * trans, int n, int incr, double * Sr, double * pin,
    double * pout) {
  double alpha = 1.;
  double beta = 1.;
  const int n3 = n * n * n;

  for (int l1 = 0; l1 < n3; l1 += n)
    blas::dgemv_(trans, &n, &n, &alpha, Sr, &n, pin + l1, &incr, &beta, pout
        + l1, &incr);
}

void fmm_data::M2M() {

  int zindex[8], zcount, yindex[4], ycount, xindex[2], xcount;

  int n2 = n * n; // n2 = n^2
  int n3 = n * n * n; // n3 = n^3

  // Initialize leaf M expansions
  int lvl = lvl_idx.m - 2;

  for (int k = lvl_idx(lvl); k < lvl_idx(lvl + 1); ++k) {

    // printf("k %d\n",k);

    oct_node * on = lvl_p(k);
    int N = on->pnt.m; // Number of points
    dmat Sr_x(N, n), Sr_y(N, n), Sr_z(N, n);
    init_Sr_cell(on, Tkmat, Sr_x, Sr_y, Sr_z);

    // Calling P2M
    zero(on->M);
    P2M(on->M, on->sigma, Sr_x, Sr_y, Sr_z);
  }

  dvec Sy(2 * n3), Sz(4 * n3);

  // Loop over all levels
  for (int lvl = lvl_idx.m - 3; lvl >= 2; --lvl) {
    for (int k = lvl_idx(lvl); k < lvl_idx(lvl + 1); ++k) {
      // printf("lvl, k: %d %d\n",lvl,k);

      oct_node * on = lvl_p(k);

      Sr_indexing(on, xcount, xindex, ycount, yindex, zcount, zindex);

      char trans[] = "t";
      int incr, j_in, j_out;

      zero(on->M);
      zero(Sy);
      zero(Sz);

      // Transform along the z-component
      for (int i = 0; i < zcount; i++) {
        Sr_idx_get(zindex[i], j_in, j_out, incr);
        Sr_zgemv(trans, n, incr, Sr.p, on->child[j_in]->M.p, Sz.p + j_out * n3);
      }

      // Transform along the y-component
      for (int i = 0; i < ycount; i++) {
        Sr_idx_get(yindex[i], j_in, j_out, incr);
        Sr_ygemv(trans, n, incr * n, Sr.p, Sz.p + j_in * n3, Sy.p + j_out * n3);
      }

      // Transform along the x-component
      for (int i = 0; i < xcount; i++) {
        Sr_idx_get(xindex[i], j_in, j_out, incr);
        Sr_xgemv(trans, n, incr * n2, Sr.p, Sy.p + j_in * n3, on->M.p);
      }
    }
  }
}

void fmm_data::L2L(dvec & phi, kernel * kfun) {

  int xindex[2], xcount, yindex[4], ycount, zindex[8], zcount;

  int n2 = n * n; // n2 = n^2
  int n3 = n * n * n; // n3 = n^3

  // Loop over all levels
  dvec Fx(2 * n3), Fy(4 * n3);

  for (int lvl = 2; lvl < lvl_idx.m - 2; ++lvl) {
    for (int k = lvl_idx(lvl); k < lvl_idx(lvl + 1); ++k) {

      oct_node * on = lvl_p(k);

      Sr_indexing(on, xcount, xindex, ycount, yindex, zcount, zindex);

      char trans[] = "n";

      int incr, j_in, j_out;

      zero(Fx);
      zero(Fy);

      // Interpolate the parent field along the x-component
      for (int i = 0; i < xcount; i++) {
        Sr_idx_get(xindex[i], j_out, j_in, incr);
        Sr_xgemv(trans, n, incr * n2, Sr.p, on->L.p, Fx.p + j_out * n3);
      }

      // Interpolate the parent field along the y-component
      for (int i = 0; i < ycount; i++) {
        Sr_idx_get(yindex[i], j_out, j_in, incr);
        Sr_ygemv(trans, n, incr * n, Sr.p, Fx.p + j_in * n3, Fy.p + j_out * n3);
      }

      /* Interpolate the parent field along the z-component and add
       to child field */
      for (int i = 0; i < zcount; i++) {
        Sr_idx_get(zindex[i], j_out, j_in, incr);
        Sr_zgemv(trans, n, incr, Sr.p, Fy.p + j_in * n3, on->child[j_out]->L.p);
      }
    }
  }

  // Process leaf L expansions
  int lvl = lvl_idx.m - 2;

  // printf("lvl_index.m %d \n", lvl_idx.m);
  // printf("lvl_index.m %f \n", (double)lvl_idx.m);

  for (int k = lvl_idx(lvl); k < lvl_idx(lvl + 1); ++k) {
    oct_node * on = lvl_p(k);
    int Nf = on->pnt.m; // Number of points
    dmat Sr_x(Nf, n), Sr_y(Nf, n), Sr_z(Nf, n);
    init_Sr_cell(on, Tkmat, Sr_x, Sr_y, Sr_z);

    dvec phi_vec(Nf);
    zero(phi_vec);

    // Compute the values at the field points
    L2P(phi_vec, on->L, Sr_x, Sr_y, Sr_z);

    // Due to near field interactions
    for (int m = 0; m < on->ngbr.m; m++) {
      oct_node *B = on->ngbr(m);
      int Ns = B->pnt.m;
      P2P(on->pnt.p, B->pnt.p, B->sigma.p, Nf, Ns, phi_vec.p, kfun);
    }

    list<int>::iterator it = on->pidx.begin();
    for (int i = 0; i < Nf; ++it, ++i)
      phi(*it) = phi_vec(i);
  }
}

void fmm_data::compute(kernel * kfun, dvec & phi) {

  // Upward pass
  M2M();

  // Compute all the cell interactions
  M2L();

  // Downward pass
  L2L(phi, kfun);
}

void fmm_data::create_lvl_indexing() { // ???
  // Number of nodes at each level
  lvl_idx.resize(l + 2);
  zero(lvl_idx);
  int curr_l = 0;
  nnode_at_level(false, curr_l, &root, lvl_idx, lvl_p);

  for (int k=0; k<l+2;k++)
  {
    printf("k, lvl_idx: %d %d\n",k,lvl_idx(k) );
  }

  // Cumulative sum
  cum_sum(lvl_idx);

  for (int k=0; k<l+2;k++)
  {
    printf("k, lvl_idx: %d %d\n",k,lvl_idx(k) );
  }


  int n_nodes = lvl_idx(lvl_idx.m - 1); // Total number of nodes

  printf("n_nodes: %d\n",n_nodes);

  lvl_p.resize(n_nodes);

  // Node IDs
  ivec nnode_(l + 1);
  for (int i = 0; i < lvl_idx.m - 1; ++i)
    nnode_(i) = lvl_idx(i);
  curr_l = 0;
  nnode_at_level(true, curr_l, &root, nnode_, lvl_p);


for (int k=0; k<l+1;k++)
  {
    printf("k, nnode_: %d %d\n",k,nnode_(k));
  }


}



void allocate_M_L(int n3, Vector<oct_node*> & lvl_p) {
  // Allocate memory for M and L data
  for (int i = 0; i < lvl_p.m; ++i) {
    oct_node * on = lvl_p(i);
    on->M.resize(n3);
    on->L.resize(n3);
  }
}

void fmm_data::setup_leaf_data(Vector<vec3> & point, dvec & q) {
  // Initialize the vector with point coordinates for leaf nodes

  // printf("lvl_idx.m %d\n", lvl_idx.m);

  for (int i = lvl_idx(lvl_idx.m - 2); i < lvl_idx(lvl_idx.m - 1); ++i) {
    oct_node * on = lvl_p(i);


    list<int> & l = on->pidx;

    // for (std::list<int>::iterator it=l.begin(); it != l.end(); ++it)
    // {
    //   std::cout << ' ' << *it;
    // }

    // std::cout << '\n'; 

    // printf("l.size() %d\n",l.size());
    

    on->pnt.resize(l.size());
    on->sigma.resize(l.size());
    list<int>::iterator it = l.begin();

    // for (std::list<int>::iterator it=l.begin(); it != l.end(); ++it)
    // {
    //   std::cout << ' ' << *it;
    // }


    int j = 0;
    for (; it != l.end(); ++it, ++j) {
      on->pnt(j) = point(*it);
      on->sigma(j) = q(*it);
    }
  }
}

void fmm_data::init_M2L_p(ivec & Kidx, dvec & M2LOp) {
  int n3 = n * n * n;
  int n6 = n3 * n3;
  // Calculate all pointers for the M2L operator
  // Loop over all levels and cells
  int l = lvl_idx.m - 2;
  for (int lvl = l; lvl >= 2; --lvl) {
    for (int i = lvl_idx(lvl); i < lvl_idx(lvl + 1); ++i) {
      oct_node * on = lvl_p(i);
      if (on->Ilist.m > 0) {
        vec3 & ctr = on->ctr;
        double S = on->S;
        // Loop over all nodes in the interaction
        on->K.resize(on->Ilist.m);
        for (int j = 0; j < on->Ilist.m; ++j) {
          // Center of cluster
          vec3 ctrj = round((1.0 / S) * (on->Ilist(j)->ctr - ctr));
          // Index in table
          int idx = M2L_index_map(ctrj);
          int idx_table = n6 * (316 * (l - lvl) + Kidx(idx));
          on->K(j) = M2LOp.p + idx_table;
        }
      }
    }
  }
}

void fmm(int l, int n, Vector<vec3> & point, dvec & q, kernel * kfun,
    dvec & phi) {

  assert(l>0);
  assert(n>0);
  assert(point.m > 0);
  assert(point.m == q.m);

  const int N = point.m;

  timeType t0, t1, t2; // Time variables

  int n3 = n * n * n; // n3 = n^3

  // Begin pre-computation timing
  t0 = Timer();

  // Calculate the center and size of the box containing all the points
  vec3 xmin(1e32, 1e32, 1e32);
  vec3 xmax(-1e32, -1e32, -1e32);
  for (int i = 0; i < N; ++i) {
    xmin.x = min(xmin.x, point(i).x);
    xmin.y = min(xmin.y, point(i).y);
    xmin.z = min(xmin.z, point(i).z);

    xmax.x = max(xmax.x, point(i).x);
    xmax.y = max(xmax.y, point(i).y);
    xmax.z = max(xmax.z, point(i).z);
  }

  printf("xmin %f %f %f\n",xmin.x,xmin.y,xmin.z);
  printf("xmax %f %f %f\n",xmax.x,xmax.y,xmax.z);


    
  vec3 ctr = 0.5 * (xmin + xmax);
  double S = norm_max(xmax - xmin);
  
  fmm_data fmmd(l, n);

  // Distribute the sources and field points and set up the interaction and neighbor lists
  fmmd.init_tree_Ilist(ctr, S, point);

  // Initialize the vectors of nodes at each level
  fmmd.create_lvl_indexing();

  // Copy data into leaf nodes
  fmmd.setup_leaf_data(point, q);

  // Allocate memory for M and L
  allocate_M_L(n3, fmmd.lvl_p);

  // Initialize the M2L operator
  ivec Kidx(343);
  dvec M2LOp;
  fmmd.init_M2L_op(kfun, Kidx, S, M2LOp);

  // Initialize M2L pointers
  fmmd.init_M2L_p(Kidx, M2LOp);

  // End pre-computation timing
  t1 = Timer();

  // Compute the field using FMM
  fmmd.compute(kfun, phi);

  // End FMM timing
  t2 = Timer();

  printf("Pre-compution time: %.3fs,\n"
	 "fmm time: %.3fs.\n", t1-t0, t2-t1);
}
