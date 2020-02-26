#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

#ifdef DISABLE_JACOBISVD
#include "lapacksvd.hpp"
#endif

using namespace std;

//20200117#define evaltime(timeval_time) (double)timeval_time.tv_sec + (double)timeval_time.tv_usec*1e-6

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

// void set_potential(double L, dvec & b, Vector<vec3> & point);

// void compute_Tk(dvec & nodes, dmat & T);

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

// struct fmm_data {
// private:
//   int l; // Number of levels in the tree
//   int n; // Order of Chebyshev polynomials
//   double epsilon;
//   double epsilon_rel;
//   int LR_mode;
//   oct_node root; // Root node of tree
//   dvec Sr; // Sr operator, stored as a vector
//   dmat Tkmat; // Tk pre-computed data



//   void M2M();
//   void M2L();
//   void L2L(dvec & phi, kernel * kfun);
//   // void L2L(dvec & phi, diag_plus_lr * dlr_fun);


// public:
//   ivec lvl_idx; // Starting index for nodes at each level
//   Vector<oct_node *> lvl_p; // Vector of node pointers

//   fmm_data(int l_, int n_, double epsilon_, double epsilon_rel_, int LR_mode_) : 
//     l(l_), n(n_), epsilon(epsilon_), epsilon_rel(epsilon_rel_), LR_mode(LR_mode_), Sr(n_ * n_), Tkmat(n_, n_) {
//   }

//   void init_tree_Ilist(vec3 & ctr, double S, Vector<vec3> & pnt);
//   void create_lvl_indexing();
//   void setup_leaf_data(Vector<vec3> & point, dvec & b);
//   void setup_x_leaf_check(Vector<vec3> & point, dvec & x_leaf_check);
//   void init_M2L_op(kernel * kfun, ivec & Kidx, double S, dvec & M2LOp);
//   // void init_M2L_op(diag_plus_lr * dlr_fun, ivec & Kidx, double S, dvec & M2LOp);
//   void init_M2L_p(ivec & Kidx, dvec & M2LOp);
//   void compute(kernel * kfun, dvec & phi);
//   // void compute(diag_plus_lr * dlr_fun, dvec & phi);

  
//   void initialize_ifmm_operators(kernel * kfun,ofstream& outfile);
//   // void initialize_ifmm_operators(diag_plus_lr * dlr_fun,ofstream& outfile);
//   void elimination(ofstream& outfile);
//   void eliminate_level(int curr_level, bool IsLeaf, ofstream& outfile);
//   void solve_top_level();
//   void downward_pass(dvec & x_iFMM,ofstream& outfile);
//   void substitution(dvec & x_iFMM,ofstream& outfile);
// };

void nnode_at_level(bool save, int l, oct_node * node, ivec & nnode_level,
    Vector<oct_node *> & nptr) {
  for (int i = 0; i < 8; ++i)
    if (node->child[i] != NULL)
      nnode_at_level(save, l + 1, node->child[i], nnode_level, nptr);

  if (save) {
    node->id = nnode_level(l);
    nptr(node->id) = node;
  }
  ++nnode_level(l);

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

void init_M2L_entries(int n, int l, double S, ivec & Kidx, kernel * kfun, dvec & M2LOp){
// void init_M2L_entries(int n, int l, double S, ivec & Kidx, diag_plus_lr * dlr_fun, dvec & M2LOp){
  int n3 = n * n * n; // n3 = n^3
  int n6 = n3 * n3;

  dvec nodes(n);

  // Compute the n Chebyshev nodes of T_n(x)
  double pi = M_PI;
  for (int m = 0; m < n; m++){
    nodes(m) = cos(pi * ((double) m + 0.5) / (double) n);
    }

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
            // init_M2L_block(S, nodes, fcenter, scenter, Kmat0, dlr_fun);
          }
        }
      }
    }
    S *= 2.0;
  }
};

void IFMM_Matrix::init_M2L_op(kernel * kfun, ivec & Kidx, double S, dvec & M2LOp) {
// void fmm_data::init_M2L_op(diag_plus_lr * dlr_fun, ivec & Kidx, double S, dvec & M2LOp) {

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

  // Precompute the M2LOp interaction matrices
  int n3 = n * n * n;
  int n6 = n3 * n3;
  if (l > 1)
    M2LOp.resize(n6 * 316 * (l - 1));
  init_M2L_entries(n, l, S, Kidx, kfun, M2LOp);
  // init_M2L_entries(n, l, S, Kidx, dlr_fun, M2LOp);

}

// Debug function to test that the tree was built correctly
// int count_points(oct_node * node) {
//   int n = 0;
//   for (int i = 0; i < 8; ++i)
//     if (node->child[i] != NULL)
//       n += count_points(node->child[i]);

//   n += node->pidx.size();

//   return n;
// }


void create_FMM_tree(int levels, int idx, vec3 & point, oct_node * node){
  vec3 ctr = node->ctr;
  double S = node->S;

  for (int l_ = levels; l_ >= 1; --l_) {

    // Index of child cell containing point
    int icell = ((point.z < ctr.z) ? 0 : 1) + 2 * ((point.y < ctr.y) ? 0 : 1)
        + 4 * ((point.x < ctr.x) ? 0 : 1);

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
};


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

  for (int i = 0; i < 8; ++i)
    if (node->child[i] != NULL)
      create_Ilist(node->child[i]);
}

void IFMM_Matrix::init_tree_Ilist(vec3 & ctr, double S, Vector<vec3> & pnt) {
  root.ctr = ctr;
  root.S = S;

  for (int i = 0; i < pnt.m; ++i)
    create_FMM_tree(l, i, pnt(i), &root);

  root.ngbr.resize(1);
  root.ngbr(0) = &root;
  create_Ilist(&root);
}

void init_M2L_block(double b_size, dvec & nodes, vec3 & c_f, vec3 & c_s, dmat & Kmat, kernel * kfun) {
// void init_M2L_block(double b_size, dvec & nodes, vec3 & c_f, vec3 & c_s, dmat & Kmat, diag_plus_lr * dlr_fun) {
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
              // Kmat(i, j) = (*dlr_fun)(source,field);
              // if (i==j){
              //   if (c_s.x == c_f.x){
              //     Kmat(i,j) = 1.0+ source.x*field.x;
              // //     // cout << "fcb" << endl;
              //   // }
              // }
              // else{
              //   Kmat(i,j) =c_s.x*c_f.x;
              // }
            }
          }
        }
      }
    }
  }
};

// void P2P(vec3 *field, vec3 *source, double *q, int Nf, int Ns, double *fieldval, kernel * kfun) {
// // void P2P(vec3 *field, vec3 *source, double *q, int Nf, int Ns, double *fieldval, diag_plus_lr * dlr_fun) {
//   int i, j;

//   for (i = 0; i < Nf; i++) {
//     // Compute the interaction between each field point and source
//     for (j = 0; j < Ns; j++) {
//       if (source[j].x != field[i].x || source[j].y != field[i].y || source[j].z
//           != field[i].z) {
//         vec3 r = source[j] - field[i];
//         fieldval[i] += q[j] * (*kfun)(r);
//         // fieldval[i] += q[j] * (*dlr_fun)(field[i],source[j]); 
//       }
//     }
//   }
// }

// void P2M(dvec & w, dvec & sigma, dmat & Sr_x, dmat & Sr_y, dmat & Sr_z) {
//   const int nx = Sr_x.m;
//   const int r = Sr_x.n;

//   int k = 0;
//   for (int l1 = 0; l1 < r; l1++) {
//     for (int l2 = 0; l2 < r; l2++) {
//       for (int l3 = 0; l3 < r; l3++, k++) {
//         for (int j = 0; j < nx; j++) {         
//           w(k) += sigma(j) * Sr_x(j, l1) * Sr_y(j, l2) * Sr_z(j, l3);
//         }
//       }
//     }
//   }
// }

// void L2P(dvec & phi, dvec & f, dmat & Sr_x, dmat & Sr_y, dmat & Sr_z) {
//   const int nx = Sr_x.m;
//   const int r = Sr_x.n;

//   for (int i = 0; i < nx; i++) {
//     int k = 0;
//     for (int l1 = 0; l1 < r; l1++) {
//       for (int l2 = 0; l2 < r; l2++) {
//         double tmp = Sr_x(i, l1) * Sr_y(i, l2);
//         for (int l3 = 0; l3 < r; l3++, k++)
//         {
//           phi(i) += f(k) * tmp * Sr_z(i, l3);
//         }
//       }
//     }
//   }
// }


void assemble_full(vec3 *location,int N,kernel * kfun,MatrixXd& A_full){
// void assemble_full(vec3 *location,int N,diag_plus_lr * dlr_fun,MatrixXd& A_full){
  A_full = MatrixXd::Zero(N,N);
  for(int i=0;i<N;++i){
    for(int j=0;j<N;++j){
        vec3 r = location[j] - location[i];
        A_full(i,j) =   (*kfun)(r);
        // A_full(i,j) =   (*dlr_fun)(location[j],location[i]);
    }
  }
}

// void fmm_data::M2L() {
//   int n3 = n * n * n;
//   char trans[] = "n";
//   double alpha = 1.;
//   int incr = 1;
//   for (int lvl = 2; lvl < lvl_idx.m - 1; ++lvl) {
//     for (int k = lvl_idx(lvl); k < lvl_idx(lvl + 1); ++k) {

//       oct_node * on = lvl_p(k);
//       zero(on->L);

//       for (int i = 0; i < on->Ilist.m; ++i) {
//         double * M = on->Ilist(i)->M.p;
//         blas::dgemv_(trans, &n3, &n3, &alpha, on->K(i), &n3, M, &incr, &alpha,
//             on->L.p, &incr);
//       }
//     }
//   }
// }

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
    for (int l3 = l1; l3 < l1 + n; l3++){
      blas::dgemv_(trans, &n, &n, &alpha, Sr, &n, pin + l3, &incr, &beta, pout
          + l3, &incr);
    }
  }
}

void Sr_zgemv(char * trans, int n, int incr, double * Sr, double * pin,
    double * pout) {
  double alpha = 1.;
  double beta = 1.;
  const int n3 = n * n * n;

  for (int l1 = 0; l1 < n3; l1 += n){
    blas::dgemv_(trans, &n, &n, &alpha, Sr, &n, pin + l1, &incr, &beta, pout
        + l1, &incr);
    }
}

// void fmm_data::M2M() {

//   int zindex[8], zcount, yindex[4], ycount, xindex[2], xcount;

//   int n2 = n * n; // n2 = n^2
//   int n3 = n * n * n; // n3 = n^3

//   // Initialize leaf M expansions
//   int lvl = lvl_idx.m - 2;

//   for (int k = lvl_idx(lvl); k < lvl_idx(lvl + 1); ++k) {

//     oct_node * on = lvl_p(k);
//     int N = on->pnt.m; // Number of points
//     dmat Sr_x(N, n), Sr_y(N, n), Sr_z(N, n);
//     init_Sr_cell(on, Tkmat, Sr_x, Sr_y, Sr_z);

//     // Calling P2M
//     zero(on->M);
//     P2M(on->M, on->sigma, Sr_x, Sr_y, Sr_z);
//   }

//   dvec Sy(2 * n3), Sz(4 * n3);

//   // Loop over all levels
//   for (int lvl = lvl_idx.m - 3; lvl >= 2; --lvl) {
//     for (int k = lvl_idx(lvl); k < lvl_idx(lvl + 1); ++k) {
      
//       oct_node * on = lvl_p(k);

//       Sr_indexing(on, xcount, xindex, ycount, yindex, zcount, zindex);

//       char trans[] = "t";
//       int incr, j_in, j_out;

//       zero(on->M);
//       zero(Sy);
//       zero(Sz);

//       // Transform along the z-component
//       for (int i = 0; i < zcount; i++) {
//         Sr_idx_get(zindex[i], j_in, j_out, incr);
//         Sr_zgemv(trans, n, incr, Sr.p, on->child[j_in]->M.p, Sz.p + j_out * n3);
//       }

//       // Transform along the y-component
//       for (int i = 0; i < ycount; i++) {
//         Sr_idx_get(yindex[i], j_in, j_out, incr);
//         Sr_ygemv(trans, n, incr * n, Sr.p, Sz.p + j_in * n3, Sy.p + j_out * n3);
        
//       }

//       // Transform along the x-component
//       for (int i = 0; i < xcount; i++) {
//         Sr_idx_get(xindex[i], j_in, j_out, incr);
//         Sr_xgemv(trans, n, incr * n2, Sr.p, Sy.p + j_in * n3, on->M.p);
//       }

//     }
//   }

// }

// void fmm_data::L2L(dvec & phi, kernel * kfun) {
// // void fmm_data::L2L(dvec & phi, diag_plus_lr * dlr_fun) {


//   int xindex[2], xcount, yindex[4], ycount, zindex[8], zcount;

//   int n2 = n * n; // n2 = n^2
//   int n3 = n * n * n; // n3 = n^3

//   // Loop over all levels
//   dvec Fx(2 * n3), Fy(4 * n3);

//   for (int lvl = 2; lvl < lvl_idx.m - 2; ++lvl) {
//     for (int k = lvl_idx(lvl); k < lvl_idx(lvl + 1); ++k) {
//       oct_node * on = lvl_p(k);

//       Sr_indexing(on, xcount, xindex, ycount, yindex, zcount, zindex);

//       char trans[] = "n";

//       int incr, j_in, j_out;

//       zero(Fx);
//       zero(Fy);


//       // Interpolate the parent field along the x-component
//       for (int i = 0; i < xcount; i++) {
//         Sr_idx_get(xindex[i], j_out, j_in, incr);
//         Sr_xgemv(trans, n, incr * n2, Sr.p, on->L.p, Fx.p + j_out * n3);
//       }

//       // Interpolate the parent field along the y-component
//       for (int i = 0; i < ycount; i++) {
//         Sr_idx_get(yindex[i], j_out, j_in, incr);
//         Sr_ygemv(trans, n, incr * n, Sr.p, Fx.p + j_in * n3, Fy.p + j_out * n3);
//       }

//        Interpolate the parent field along the z-component and add
//        to child field 
//       for (int i = 0; i < zcount; i++) {
//         Sr_idx_get(zindex[i], j_out, j_in, incr);
//         Sr_zgemv(trans, n, incr, Sr.p, Fy.p + j_in * n3, on->child[j_out]->L.p);
//       }
//     }
//   }

//   // Process leaf L expansions
//   int lvl = lvl_idx.m - 2;

//   for (int k = lvl_idx(lvl); k < lvl_idx(lvl + 1); ++k) {
//     oct_node * on = lvl_p(k);
//     int Nf = on->pnt.m; // Number of points
//     dmat Sr_x(Nf, n), Sr_y(Nf, n), Sr_z(Nf, n);
//     init_Sr_cell(on, Tkmat, Sr_x, Sr_y, Sr_z);

//     dvec phi_vec(Nf);
//     zero(phi_vec);

//     // Compute the values at the field points
//     L2P(phi_vec, on->L, Sr_x, Sr_y, Sr_z);

//     // Due to near field interactions
//     for (int m = 0; m < on->ngbr.m; m++) {
//       oct_node *B = on->ngbr(m);
//       int Ns = B->pnt.m;
//       P2P(on->pnt.p, B->pnt.p, B->sigma.p, Nf, Ns, phi_vec.p, kfun);
//       // P2P(on->pnt.p, B->pnt.p, B->sigma.p, Nf, Ns, phi_vec.p, dlr_fun);
//     }
    
//     list<int>::iterator it = on->pidx.begin();
//     for (int i = 0; i < Nf; ++it, ++i)
//       phi(*it) = phi_vec(i);
//   }
// }

// void fmm_data::compute(kernel * kfun, dvec & phi) {
// // void fmm_data::compute(diag_plus_lr * dlr_fun, dvec & phi) {


//   // Upward pass
//   M2M();

//   // Compute all the cell interactions
//   M2L();

//   // Downward pass
//   L2L(phi, kfun);
//   // L2L(phi, dlr_fun);
// }

// void IFMM_Matrix::initialize_ifmm_operators(kernel * kfun,ofstream& outfile){
void IFMM_Matrix::initialize_ifmm_operators(kernel * kfun,ofstream& outfile){
// void fmm_data::initialize_ifmm_operators(diag_plus_lr * dlr_fun,ofstream& outfile){


  // LEAF LEVEL
  int iLevel=l;
  bool IsLeaf=true;

  int n3=n*n*n;

#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
    oct_node * on = lvl_p(k);
    on->IsEliminated=0;


    // set number of degrees of freedom
    set_nDOF(on,n,IsLeaf);

    dmat Sr_x(on->nDOF, n), Sr_y(on->nDOF, n), Sr_z(on->nDOF, n);
    init_Sr_cell(on, Tkmat, Sr_x, Sr_y, Sr_z);

    // set_P2M_operator
    on->P2M_operator = MatrixXd::Zero(n3,on->nDOF);
    set_P2M_operator(on->P2M_operator,Sr_x, Sr_y, Sr_z);

    // set_L2P_operator
    on->L2P_operator = MatrixXd::Zero(on->nDOF,n3);
    set_L2P_operator(on->L2P_operator,Sr_x, Sr_y, Sr_z);

    // set M2L_operator
    set_M2L_operator(on,n);

    // set P2P_operator
    set_P2P_operator(on,kfun,IsLeaf);
    // set_P2P_operator(on,dlr_fun,IsLeaf);

      
    // set P2L_operator
    set_P2L_operator(on,n,IsLeaf);
    
    // set M2P_operator
    set_M2P_operator(on,n,IsLeaf);

    // set RHS
    set_RHS(on,IsLeaf);

  }


  // NON-LEAF LEVELS (UP TO LEVEL 2)
  for (int iLevel=l-1;iLevel>=2;iLevel--){
    IsLeaf=false;
    // int n3=n*n*n;

#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
      oct_node * on = lvl_p(k);
      on->IsEliminated=0;

      on->rank = n3;

      // on->P2M_operator = MatrixXd::Zero(n3,n3);
      set_M2M_operator(on,n,Sr); // P2M at non-leaf levels corresponds to M2M

      // on->L2P_operator = MatrixXd::Zero(n3,n3);
      set_L2L_operator(on,n,Sr); // L2P at non-leaf levels corresponds to L2L

    }
  }

  // CHECK IF (nDOF > RANK) AT THE LEAF LEVEL
  iLevel=l;

#if defined(_OPENMP) && defined(PARA_CHECK_RANK) // #if defined(CHECK_RANK_SELF)

#pragma omp parallel for
  for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++ k) {
    oct_node *on = lvl_p(k);
    check_rank_self_setup(on, epsilon_rel_basis);
  }

#pragma omp parallel for
  for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++ k) {
    oct_node *on = lvl_p(k);
    check_rank_self(on, iLevel);
  }

#else // ! _OPENMP || ! PARA_CHECK_RANK

  for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {

    oct_node * on = lvl_p(k);
    outfile << "Node: " << on->id << ", Rank (inital): " << on->rank << endl;
    check_rank(on,iLevel,epsilon_rel_basis);
    outfile << "Node: " << on->id << ", Rank (after SVD): " << on->rank << endl;

  }

#endif // ! _OPENMP || ! PARA_CHECK_RANK

  for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {   
    oct_node * on = lvl_p(k);   
    get_weights(on,iLevel,epsilon_rel_basis);   
  }

  // NON-LEAF LEVELS (UP TO LEVEL 2)
  for (int iLevel=l-1;iLevel>=2;iLevel--){
    IsLeaf=false;

#if defined(_OPENMP)
#pragma omp parallel for
#endif
   for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
      oct_node * on = lvl_p(k);
      on->IsEliminated=0;

      // int n3=n*n*n;

      // set number of degrees of freedom
      set_nDOF(on,n,IsLeaf);

      // set M2L_operator
      set_M2L_operator(on,n);

      // no initial P2P_operator at non-leaf levels (will be constructed based on M2L of the children)

      // set RHS
      set_RHS(on,IsLeaf);
    }
  }

}

void IFMM_Matrix::estimate_lsv(double & sigma0_Aextended, int & leafLevel, ofstream& outfile){

  // NUMBER OF NODES PER LEVEL  
  ivec nNodes_per_level(leafLevel-1);
  for (int iLevel=2; iLevel<=leafLevel; iLevel++){
    nNodes_per_level(iLevel-2)=lvl_idx(iLevel+1)-lvl_idx(iLevel);
  }


  // TOTAL NUMBER OF VARIABLES
  int nVariables = 0;
  for (int iLevel=2; iLevel<leafLevel; iLevel++){
    nVariables += 2*nNodes_per_level(iLevel-2);
  }
  nVariables += 3*nNodes_per_level(leafLevel-2);

  // ivec nVariables_cum(leafLevel+1);
  // nVariables_cum(0)=0;
  // nVariables_cum(1)=2*nNodes_per_level(leafLevel-2);
  // for (int iLevel=leafLevel-1; iLevel >=2 ; iLevel--){
  //   nVariables_cum(leafLevel-iLevel+1) = nVariables_cum(leafLevel-iLevel) + nNodes_per_level(iLevel-2)+nNodes_per_level(iLevel+1-2);
  // }
  // nVariables_cum(leafLevel)=nVariables_cum(leafLevel-1)+nNodes_per_level(0);



  // for (int iLevel=0;iLevel<lvl_idx.m;iLevel++){
  //   printf("lvl_idx(iLevel): %d \n",lvl_idx(iLevel));
  // }
  // for (int iLevel=0;iLevel<nNodes_per_level.m;iLevel++){
  //   printf("nNodes_per_level(iLevel): %d \n",nNodes_per_level(iLevel));
  // }
  // printf("nVariables: %d\n",nVariables);

  // CUMULATIVE VECTOR
  ivec dimensions_cum(nVariables+1);
  dimensions_cum(0)=0;

  ivec map_index_to_id(nVariables);


  // LEAF LEVEL
  int iLevel=leafLevel;
  int nNodes_before=lvl_idx(iLevel);
  // int nNodes_before=0;
  // for (int jLevel=iLevel-1;jLevel>=2;jLevel--){
    // nNodes_before+=nNodes_per_level(jLevel);
  // }
  // printf("nNodes_before: %d\n",nNodes_before);


  for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
    // CURRENT NODE
    oct_node * on = lvl_p(k);
    // printf("k,2*(k-nNodes_before) + 2: %d %d\n",k,2*(k-nNodes_before) + 2);
    dimensions_cum(2*(k-nNodes_before) + 1) =  dimensions_cum(2*(k-nNodes_before)) + on->nDOF;
    dimensions_cum(2*(k-nNodes_before) + 2) =  dimensions_cum(2*(k-nNodes_before)+1) + on->rank;
    map_index_to_id(2*(k-nNodes_before)) = on->id;
    map_index_to_id(2*(k-nNodes_before)+1) = on->id;


    on->dim_cum_index_x = 2*(k-nNodes_before);
    on->dim_cum_index_z = 2*(k-nNodes_before)+1;
  }
  if (iLevel==2){
    for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
      // CURRENT NODE
      oct_node * on = lvl_p(k);
      dimensions_cum(2*nNodes_per_level(iLevel-2) + k-nNodes_before + 1)= dimensions_cum(2*nNodes_per_level(iLevel-2) + k-nNodes_before) + on->rank;
      map_index_to_id(2*nNodes_per_level(iLevel-2) + k-nNodes_before) = on->id;

      on->dim_cum_index_y=2*nNodes_per_level(iLevel-2) + k-nNodes_before;

    }
  }

  



  int index = 2*(lvl_idx(iLevel + 1) -1 -nNodes_before) + 2;
  // NON-LEAF LEVELS
  for (int iLevel=leafLevel-1;iLevel>=2;iLevel--){
    for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
      // CURRENT NODE
      oct_node * on = lvl_p(k);
      for (int iChild=0;iChild<8;iChild++){
        if (on->child[iChild] != NULL){
          dimensions_cum(index+1) = dimensions_cum(index) + on->child[iChild]->rank;
          map_index_to_id(index) = on->child[iChild]->id;
          on->child[iChild]->dim_cum_index_y=index;
          index++;
        }
      }
      dimensions_cum(index+1) = dimensions_cum(index) + on->rank;
      map_index_to_id(index) = on->id;
      on->dim_cum_index_z=index;
      index++;

    }

    if (iLevel==2){
      for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
        // CURRENT NODE
        oct_node * on = lvl_p(k);
        dimensions_cum(index + 1)= dimensions_cum(index) + on->rank;
        map_index_to_id(index) = on->id;
        on->dim_cum_index_y=index;
        index++;
      }
    }
  }

  // for (int iLevel=leafLevel;iLevel>=2;iLevel--){
  //   printf("Level: %d\n",iLevel);
  //   for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
  //     // CURRENT NODE
  //     oct_node * on = lvl_p(k);

  //     printf("on->dim_cum_index_x: %d\n",on->dim_cum_index_x);
  //     printf("on->dim_cum_index_y: %d\n",on->dim_cum_index_y);
  //     printf("on->dim_cum_index_z: %d\n",on->dim_cum_index_z);
  //   }
  // }


  // for (int i=0; i<nVariables;i++){
  //   printf("i,map_index_to_id(i): %d %d \n",i,map_index_to_id(i));
  // }

  // for (int i=0; i<nVariables+1;i++){
  //   printf("dimensions_cum(i): %d \n",dimensions_cum(i));
  // }
  // for (int i=0; i<nVariables;i++){
  //   printf("map_index_to_id: %d \n",map_index_to_id(i));
  // }
  // for (int k = lvl_idx(curr_level); k < lvl_idx(curr_level + 1); ++k) {
  //   oct_node * on = lvl_p(k);
  //   rank_cum(k-nNodes_before+1)=rank_cum(k-nNodes_before)+on->rank;
  // }

  // // // // // // // // // // //
  // ESTIMATE SINGULAR VALUES //
  // // // // // // // // // // //
  int rank_RSVD = 10;
  MatrixXd X = MatrixXd::Random(dimensions_cum(nVariables),rank_RSVD);
  MatrixXd AX = MatrixXd::Zero(dimensions_cum(nVariables),rank_RSVD);
  // A_extended_multiplication(X,AX,map_index_to_id,dimensions_cum);
  // cout << "AX.norm()" << endl << AX.norm() << endl;

  
  // Gaussian Random Matrix for A_extended^T
  MatrixXd O(dimensions_cum(nVariables), rank_RSVD);
  sample_gaussian(O);
  
  // Compute Sample Matrix of A^T
  // MatrixXd Y = A.transpose() * O;
  MatrixXd Y = MatrixXd::Zero(dimensions_cum(nVariables),rank_RSVD);
  A_extended_transpose_multiplication(O,Y,map_index_to_id,dimensions_cum,leafLevel);

  // Orthonormalize Y
  gram_schmidt(Y);
  
  // Range(B) = Range(A^T)
  // DenseMatrix B = A * Y;
  MatrixXd B = MatrixXd::Zero(dimensions_cum(nVariables),rank_RSVD);
  A_extended_multiplication(Y,B,map_index_to_id,dimensions_cum,leafLevel,outfile);
  
  // Gaussian Random Matrix
  // DenseMatrix P(B.cols(), r);
  MatrixXd P(rank_RSVD, rank_RSVD);
  sample_gaussian(P);
  
  // Compute Sample Matrix of B
  MatrixXd Z = B * P;
  
  // Orthonormalize Z
  gram_schmidt(Z);
  
  // Range(C) = Range(B)
  MatrixXd C = Z.transpose() * B; 
  
  // JacobiSVD<MatrixXd> svd_C(C, ComputeThinU | ComputeThinV);
#ifdef DISABLE_JACOBISVD 
  LapackSVD<MatrixXd> svd_C(C);
#else
  JacobiSVD<MatrixXd> svd_C(C);
#endif
  VectorXd Sigma_A_extended = svd_C.singularValues();
  // cout << "Sigma_A_extended:" << endl << Sigma_A_extended << endl; 
  sigma0_Aextended = Sigma_A_extended(0);

 

  



  // CHECK THE IMPLEMENTATION
  // int rank_RSVD = 1;
  // MatrixXd X  = MatrixXd::Zero(dimensions_cum(nVariables),rank_RSVD);
  // MatrixXd AX = MatrixXd::Zero(dimensions_cum(nVariables),rank_RSVD);
  // MatrixXd B  = MatrixXd::Zero(dimensions_cum(nVariables),rank_RSVD);


  // iLevel=l;
  // bool IsLeaf=true;
  // for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
  //   // CURRENT NODE
  //   oct_node * on = lvl_p(k);
  //   set_RHS(on,IsLeaf);

  //   int index_m2i_iSelf_x=0;
  //   for (int i=0; i<nVariables; i++){
  //     if (map_index_to_id(i)==on->id){
  //       index_m2i_iSelf_x = i;
  //       break;
  //     }
  //   }
  //   int index_m2i_iSelf_z = index_m2i_iSelf_x +1;
    
  //   X.block(dimensions_cum(index_m2i_iSelf_x),0,on->nDOF,rank_RSVD) = on->x;
  //   B.block(dimensions_cum(index_m2i_iSelf_x),0,on->nDOF,rank_RSVD) = on->RHS_leaf; //  (has been changed!)
  //   X.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,rank_RSVD) = on->z;
  
  //   if (iLevel==2){
  //     int index_m2i_iSelf_y=0;
  //     for (int i=0; i<nVariables; i++){
  //       if (map_index_to_id(i)==on->id){
  //         index_m2i_iSelf_y = i;
  //       }
  //     }
  //     X.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,rank_RSVD) = on->y;
  //   }
  // }

  // for (int iLevel=l-1;iLevel>=2;iLevel--){
  //   for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
  //     // CURRENT NODE
  //     oct_node * on = lvl_p(k);

  //     int index_m2i_iSelf_z=0;
  //     for (int i=0; i<nVariables; i++){
  //       if (map_index_to_id(i)==on->id){
  //         index_m2i_iSelf_z = i;
  //         break;
  //       }
  //     }

  //     for (int iChild=0;iChild<8;iChild++){
  //       if (on->child[iChild] != NULL){

  //         int index_m2i_iChild_y=0;
  //         for (int i=0; i<nVariables; i++){
  //           if (map_index_to_id(i)==on->child[iChild]->id){
  //             index_m2i_iChild_y = i;
  //           }
  //         }
  //         X.block(dimensions_cum(index_m2i_iChild_y),0,on->child[iChild]->rank,rank_RSVD) = on->child[iChild]->y;
  //       }
  //     }
  //     X.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,rank_RSVD) = on->z;

  //     if (iLevel==2){
  //       int index_m2i_iSelf_y=0;
  //       for (int i=0; i<nVariables; i++){
  //         if (map_index_to_id(i)==on->id){
  //           index_m2i_iSelf_y = i;
  //         }
  //       }
  //       X.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,rank_RSVD) = on->y;
  //     }
  //   }
  // }



  // A_extended_multiplication(X,AX,map_index_to_id,dimensions_cum);

  // int test=0;
  // for (int i=0; i<dimensions_cum(nVariables);i++){
  //   if (abs(AX(i)-B(i))<1e-10){
  //     test=1;
  //   }
  //   else{
  //     test=0;
  //   }
  //   printf("X(i),AX(i),B(i),AX(i)/B(i),AX(i)-B(i),test: %20.5e %20.5e %20.5e %20.5e %20.5e %d\n",X(i),AX(i),B(i),AX(i)/B(i),AX(i)-B(i),test);
  // }

  // VectorXd diff = B-AX;
  // cout << "diff.norm()/B.norm(): " << diff.norm()/B.norm() << endl;



}

void IFMM_Matrix::A_extended_multiplication(MatrixXd &X, MatrixXd &AX, ivec &map_index_to_id, ivec &dimensions_cum, int & leafLevel, ofstream& outfile){

  // cout << "Start A_extended_multiplication..." << endl;

  // clock_t start_index, stop_index; double Time_index=0;
  // clock_t start_oper, stop_oper; double Time_oper=0;

  // int nVariables = map_index_to_id.m;
  int nCol = X.cols();

  // LEAF LEVEL
  int iLevel=leafLevel; 
  // printf("iLevel: %d\n",iLevel);
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
    // CURRENT NODE
    oct_node * on = lvl_p(k);

    // start_index = clock();
    // int index_m2i_iSelf_x=0;
    // // for (int i=0; i<nVariables; i++){
    // for (int i=nVariables_cum(leafLevel-iLevel); i<nVariables_cum(leafLevel-iLevel+1); i++){
    //   if (map_index_to_id(i)==on->id){
    //     index_m2i_iSelf_x = i;
    //     break;
    //   }
    // }
    // int index_m2i_iSelf_z = index_m2i_iSelf_x+1;
    // int index_m2i_iSelf_y=0;
    // // for (int i=0; i<nVariables; i++){
    // // int found=0;
    // for (int i=nVariables_cum(leafLevel-iLevel+1); i<nVariables_cum(leafLevel-iLevel+2); i++){
    //   if (map_index_to_id(i)==on->id){
    //     index_m2i_iSelf_y = i;
    //     // found++;
    //     // if (found==3){
    //       // break;
    //     // }
    //   }
    // }
    // stop_index = clock();
    // Time_index += double(stop_index-start_index)/double(CLOCKS_PER_SEC);

    // EQUATION XI
    for (int iNgbr = 0; iNgbr < on->ngbr.m; ++iNgbr) {
      // start_index = clock();
      int index_iNgbr_to_iSelf=0;
      for (int i=0; i < on->ngbr(iNgbr)->ngbr.m;i++){
        if (on->ngbr(iNgbr)->ngbr(i)->id == on->id)
        {
          index_iNgbr_to_iSelf=i; // index to current node (in the neighbour list of i)
          break;
         }
      }
      // int index_m2i_iNgbr_to_iSelf=0;
      // // for (int i=0; i<nVariables; i++){
      // for (int i=nVariables_cum(leafLevel-iLevel); i<nVariables_cum(leafLevel-iLevel+1); i++){
      //   if (map_index_to_id(i)==on->ngbr(iNgbr)->id){
      //     index_m2i_iNgbr_to_iSelf = i;
      //     break;
      //   }
      // }
      // stop_index = clock();
      // Time_index += double(stop_index-start_index)/double(CLOCKS_PER_SEC);

      // start_oper = clock();
      // AX.block(dimensions_cum(index_m2i_iSelf_x),0,on->nDOF,nCol).noalias() += on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf]*X.block(dimensions_cum(index_m2i_iNgbr_to_iSelf),0,on->ngbr(iNgbr)->nDOF,nCol);
      // printf("index_m2i_iSelf_x,on->dim_cum_index_x: %d %d\n",index_m2i_iSelf_x,on->dim_cum_index_x);
      // printf("index_m2i_iNgbr_to_iSelf,on->ngbr(iNgbr)->dim_cum_index_x: %d %d\n",index_m2i_iNgbr_to_iSelf,on->ngbr(iNgbr)->dim_cum_index_x);

      AX.block(dimensions_cum(on->dim_cum_index_x),0,on->nDOF,nCol).noalias() += on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf]*X.block(dimensions_cum(on->ngbr(iNgbr)->dim_cum_index_x),0,on->ngbr(iNgbr)->nDOF,nCol);
      // stop_oper = clock();
      // Time_oper += double(stop_oper-start_oper)/double(CLOCKS_PER_SEC);

    }
    // start_oper = clock();
    // AX.block(dimensions_cum(index_m2i_iSelf_x),0,on->nDOF,nCol).noalias() += on->L2P_operator*X.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol);
    // printf("index_m2i_iSelf_z,on->dim_cum_index_z: %d %d\n",index_m2i_iSelf_z,on->dim_cum_index_z);
    AX.block(dimensions_cum(on->dim_cum_index_x),0,on->nDOF,nCol).noalias() += on->L2P_operator*X.block(dimensions_cum(on->dim_cum_index_z),0,on->rank,nCol);
    // stop_oper = clock();
    // Time_oper += double(stop_oper-start_oper)/double(CLOCKS_PER_SEC);

    // EQUATION ZI
    // start_oper = clock();
    // AX.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol).noalias() = -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,nCol) + on->P2M_operator*X.block(dimensions_cum(index_m2i_iSelf_x),0,on->nDOF,nCol);
    AX.block(dimensions_cum(on->dim_cum_index_z),0,on->rank,nCol).noalias() = -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(on->dim_cum_index_y),0,on->rank,nCol) + on->P2M_operator*X.block(dimensions_cum(on->dim_cum_index_x),0,on->nDOF,nCol);
    // stop_oper = clock();
    // Time_oper += double(stop_oper-start_oper)/double(CLOCKS_PER_SEC);
  
    if (iLevel==2){
      // EQUATION YI
      for (int lInteract=0; lInteract < on->Ilist.m; lInteract++){
        // start_index = clock();
        int index_lInteract_to_iSelf=0;
        for (int i=0; i < on->Ilist(lInteract)->Ilist.m;i++){
          if (on->Ilist(lInteract)->Ilist(i)->id == on->id)    
          {
            index_lInteract_to_iSelf=i; // index to current node (in the interaction list)
            break;
          }
        }
        // int index_m2i_lInteract_to_iSelf=0;
        // for (int i=0; i<nVariables; i++){
        //   if (map_index_to_id(i)==on->Ilist(lInteract)->id){
        //     index_m2i_lInteract_to_iSelf = i;
        //   }
        // }
        // stop_index = clock();
        // Time_index += double(stop_index-start_index)/double(CLOCKS_PER_SEC);

        // start_oper = clock();
        // AX.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,nCol).noalias() += on->Ilist(lInteract)->M2L_operator[index_lInteract_to_iSelf]*X.block(dimensions_cum(index_m2i_lInteract_to_iSelf),0,on->Ilist(lInteract)->rank,nCol);
        AX.block(dimensions_cum(on->dim_cum_index_y),0,on->rank,nCol).noalias() += on->Ilist(lInteract)->M2L_operator[index_lInteract_to_iSelf]*X.block(dimensions_cum(on->Ilist(lInteract)->dim_cum_index_y),0,on->Ilist(lInteract)->rank,nCol);
        // stop_oper = clock();
        // Time_oper += double(stop_oper-start_oper)/double(CLOCKS_PER_SEC);

      }

      // start_oper = clock();
      // AX.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,nCol).noalias() += -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol);
      AX.block(dimensions_cum(on->dim_cum_index_y),0,on->rank,nCol).noalias() += -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(on->dim_cum_index_z),0,on->rank,nCol);
      // stop_oper = clock();
      // Time_oper += double(stop_oper-start_oper)/double(CLOCKS_PER_SEC);
    }
  }


  
  // NON-LEAF LEVELS
  for (int iLevel=leafLevel-1;iLevel>=2;iLevel--){ 

    // printf("iLevel: %d\n",iLevel);
#if defined(_OPENMP)
#pragma omp parallel for
#endif
 for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {

      // CURRENT NODE
      oct_node * on = lvl_p(k);

      // printf("k: %d\n",k);
      // start_index = clock();
      // int index_m2i_iSelf_z=0;
      // // for (int i=0; i<nVariables; i++){
      // for (int i=nVariables_cum(leafLevel-iLevel); i<nVariables_cum(leafLevel-iLevel+1); i++){
      //   if (map_index_to_id(i)==on->id){
      //     index_m2i_iSelf_z = i;
      //     break;
      //   }
      // }
      // int index_m2i_iSelf_y=0;
      // for (int i=0; i<nVariables; i++){
      //   if (map_index_to_id(i)==on->id){
      //     index_m2i_iSelf_y = i;
      //   }
      // }
      // stop_index = clock();
      // Time_index += double(stop_index-start_index)/double(CLOCKS_PER_SEC);

      // EQUATION YI OF CHILDREN
      // bool firstChild=false;
      // int index_m2i_iChild_y_firstChild = 0;
      for (int iChild=0;iChild<8;iChild++){
        if (on->child[iChild] != NULL){

          // printf("on->id,on->child[iChild]->id: %d %d\n",on->id,on->child[iChild]->id);
          // start_index = clock();
          // int index_m2i_iChild_z=0;
          // // for (int i=0; i<nVariables; i++){
          // //   if (map_index_to_id(i)==on->child[iChild]->id){
          // //     if (iLevel+1==l){
          // //       index_m2i_iChild_z = i+1;
          // //       break;
          // //     }
          // //     else{
          // //       index_m2i_iChild_z = i;
          // //     }
          // //   }
          // // }
          // if (iLevel+1==l){
          //   for (int i=nVariables_cum(leafLevel-iLevel-1); i<nVariables_cum(leafLevel-iLevel); i++){
          //     if (map_index_to_id(i)==on->child[iChild]->id){
          //       index_m2i_iChild_z = i+1;
          //       break;
          //     }
          //   }
          // }
          // else{
          //   for (int i=nVariables_cum(leafLevel-iLevel); i<nVariables_cum(leafLevel-iLevel+1); i++){
          //      if (map_index_to_id(i)==on->child[iChild]->id){
          //         index_m2i_iChild_z = i;
          //         break;
          //      }
          //   }
          // }

          // int index_m2i_iChild_y=0;
          // // for (int i=0; i<nVariables; i++){
          // for (int i=nVariables_cum(leafLevel-iLevel); i<nVariables_cum(leafLevel-iLevel+1); i++){
          //   if (map_index_to_id(i)==on->child[iChild]->id){
          //     index_m2i_iChild_y = i;
          //     break;
          //   }
          // }
          // stop_index = clock();
          // Time_index += double(stop_index-start_index)/double(CLOCKS_PER_SEC);

          // start_index = clock();
          // if (firstChild==false){
          //   firstChild=true;
          //   index_m2i_iChild_y_firstChild = index_m2i_iChild_y;
          // }
          // stop_index = clock();
          // Time_index += double(stop_index-start_index)/double(CLOCKS_PER_SEC);

          for (int lInteract=0; lInteract < on->child[iChild]->Ilist.m; lInteract++){
            // start_index = clock();
            int index_lInteract_to_iChild=0;
            for (int i=0; i < on->child[iChild]->Ilist(lInteract)->Ilist.m;i++){
              if (on->child[iChild]->Ilist(lInteract)->Ilist(i)->id == on->child[iChild]->id)    
              {
                index_lInteract_to_iChild=i; // index to current node (in the interaction list)
                break;
              }
            }

            // int index_m2i_lInteract_to_iChild=0;
            // // for (int i=0; i<nVariables; i++){
            // for (int i=nVariables_cum(leafLevel-iLevel); i<nVariables_cum(leafLevel-iLevel+1); i++){
            //   if (map_index_to_id(i)==on->child[iChild]->Ilist(lInteract)->id){
            //     index_m2i_lInteract_to_iChild = i;
            //     break;
            //   }
            // }
            // stop_index = clock();
            // Time_index += double(stop_index-start_index)/double(CLOCKS_PER_SEC);

            // start_oper = clock();
            // AX.block(dimensions_cum(index_m2i_iChild_y),0,on->child[iChild]->rank,nCol).noalias() += on->child[iChild]->Ilist(lInteract)->M2L_operator[index_lInteract_to_iChild]*X.block(dimensions_cum(index_m2i_lInteract_to_iChild),0,on->child[iChild]->Ilist(lInteract)->rank,nCol);
            AX.block(dimensions_cum(on->child[iChild]->dim_cum_index_y),0,on->child[iChild]->rank,nCol).noalias() += on->child[iChild]->Ilist(lInteract)->M2L_operator[index_lInteract_to_iChild]*X.block(dimensions_cum(on->child[iChild]->Ilist(lInteract)->dim_cum_index_y),0,on->child[iChild]->Ilist(lInteract)->rank,nCol);
            // stop_oper = clock();
            // Time_oper += double(stop_oper-start_oper)/double(CLOCKS_PER_SEC);
          }

          // start_oper = clock();
          // printf("index_m2i_iChild_z,on->child[iChild]->dim_cum_index_z: %d %d\n",index_m2i_iChild_z,on->child[iChild]->dim_cum_index_z);
          // AX.block(dimensions_cum(index_m2i_iChild_y),0,on->child[iChild]->rank,nCol).noalias() +=  -MatrixXd::Identity(on->child[iChild]->rank,on->child[iChild]->rank)*X.block(dimensions_cum(index_m2i_iChild_z),0,on->child[iChild]->rank,nCol);
          AX.block(dimensions_cum(on->child[iChild]->dim_cum_index_y),0,on->child[iChild]->rank,nCol).noalias() +=  -MatrixXd::Identity(on->child[iChild]->rank,on->child[iChild]->rank)*X.block(dimensions_cum(on->child[iChild]->dim_cum_index_z),0,on->child[iChild]->rank,nCol);
          // stop_oper = clock();
          // Time_oper += double(stop_oper-start_oper)/double(CLOCKS_PER_SEC);

          // start_index = clock();
          // int index_m2i_parent_to_iSelf=0;
          // for (int i=0; i<nVariables; i++){
          //   if (map_index_to_id(i)==on->parent->id){
          //     index_m2i_parent_to_iSelf = i;
          //     break;
          //   }
          // }
          // stop_index = clock();
          // Time_index += double(stop_index-start_index)/double(CLOCKS_PER_SEC);
          // cout << "a5" << endl;
          // cout << "on->child[iChild]->rank" << on->child[iChild]->rank << endl;
          // cout << "on->L2L_operator[iChild].rows()" << on->L2L_operator[iChild].rows() << endl;
          // cout << "on->L2L_operator[iChild].cols()" << on->L2L_operator[iChild].cols() << endl;
          // cout << "on->rank" << on->rank << endl;
          
          // start_oper = clock();
          // AX.block(dimensions_cum(index_m2i_iChild_y),0,on->child[iChild]->rank,nCol).noalias() += on->L2L_operator[iChild]*X.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol);
          AX.block(dimensions_cum(on->child[iChild]->dim_cum_index_y),0,on->child[iChild]->rank,nCol).noalias() += on->L2L_operator[iChild]*X.block(dimensions_cum(on->dim_cum_index_z),0,on->rank,nCol);
          // stop_oper = clock();
          // Time_oper += double(stop_oper-start_oper)/double(CLOCKS_PER_SEC);
          // MatrixXd temp = on->L2P_operator*X.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol);
          // AX.block(dimensions_cum(index_m2i_iChild_y),0,on->child[iChild]->rank,nCol).noalias() += temp.block(dimensions_cum(index_m2i_iChild_y)-dimensions_cum(index_m2i_iChild_y_firstChild),0,on->child[iChild]->rank,nCol);
        }
      }

      

      // EQUATION ZI
      // start_oper = clock();
      // AX.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol).noalias() = -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,nCol);
      AX.block(dimensions_cum(on->dim_cum_index_z),0,on->rank,nCol).noalias() = -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(on->dim_cum_index_y),0,on->rank,nCol);
      // stop_oper = clock();
      // Time_oper += double(stop_oper-start_oper)/double(CLOCKS_PER_SEC);

      // int dimension_local=0;
      // for (int iChild=0; iChild<8; iChild++){
      //   if (on->child[iChild] != NULL){
      //     start_index = clock();
      //     int index_m2i_iChild_to_iSelf_y=0;
      //     // for (int i=0; i<nVariables; i++){
      //     for (int i=nVariables_cum(leafLevel-iLevel-1); i<nVariables_cum(leafLevel-iLevel); i++){
      //       if (map_index_to_id(i)==on->child[iChild]->id){
      //         index_m2i_iChild_to_iSelf_y = i;
      //         break;
      //         }
      //     }
      //     stop_index = clock();
      //     Time_index += double(stop_index-start_index)/double(CLOCKS_PER_SEC);

      //     dimension_local+=on->child[iChild]->rank;
      //   }
      // }
      // MatrixXd temp = MatrixXd::Zero(dimension_local,nCol);
      for (int iChild=0; iChild<8; iChild++){
        if (on->child[iChild] != NULL){
          // start_index = clock();
          // int index_m2i_iChild_to_iSelf_y=0;
          // // for (int i=0; i<nVariables; i++){
          // for (int i=nVariables_cum(leafLevel-iLevel); i<nVariables_cum(leafLevel-iLevel+1); i++){
          //   if (map_index_to_id(i)==on->child[iChild]->id){
          //     index_m2i_iChild_to_iSelf_y = i;
          //     break;
          //     }
          // }
          // stop_index = clock();
          // Time_index += double(stop_index-start_index)/double(CLOCKS_PER_SEC);

          // start_oper = clock();
          // AX.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol).noalias() += on->M2M_operator[iChild]*X.block(dimensions_cum(index_m2i_iChild_to_iSelf_y),0,on->child[iChild]->rank,nCol);
          AX.block(dimensions_cum(on->dim_cum_index_z),0,on->rank,nCol).noalias() += on->M2M_operator[iChild]*X.block(dimensions_cum(on->child[iChild]->dim_cum_index_y),0,on->child[iChild]->rank,nCol);
          // stop_oper = clock();
          // Time_oper += double(stop_oper-start_oper)/double(CLOCKS_PER_SEC);
          // temp.block(dimensions_cum(index_m2i_iChild_to_iSelf_y)-dimensions_cum(index_m2i_iChild_y_firstChild),0,on->child[iChild]->rank,nCol) = X.block(dimensions_cum(index_m2i_iChild_to_iSelf_y),0,on->child[iChild]->rank,nCol);
        }
      }
      // AX.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol).noalias() += on->P2M_operator*temp;

      

      if (iLevel==2){
        // EQUATION YI
        for (int lInteract=0; lInteract < on->Ilist.m; lInteract++){
          // start_index = clock();
          int index_lInteract_to_iSelf=0;
          for (int i=0; i < on->Ilist(lInteract)->Ilist.m;i++){
            if (on->Ilist(lInteract)->Ilist(i)->id == on->id)    
            {
              index_lInteract_to_iSelf=i; // index to current node (in the interaction list)
              break;
            }
          }

          // int index_m2i_lInteract_to_iSelf=0;
          // // for (int i=0; i<nVariables; i++){
          // for (int i=nVariables_cum(leafLevel-1); i<nVariables_cum(leafLevel); i++){
          //   if (map_index_to_id(i)==on->Ilist(lInteract)->id){
          //     index_m2i_lInteract_to_iSelf = i;
          //     break;
          //   }
          // }
          // stop_index = clock();
          // Time_index += double(stop_index-start_index)/double(CLOCKS_PER_SEC);

          // start_oper = clock();
          // AX.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,nCol).noalias() += on->Ilist(lInteract)->M2L_operator[index_lInteract_to_iSelf]*X.block(dimensions_cum(index_m2i_lInteract_to_iSelf),0,on->Ilist(lInteract)->rank,nCol);
          AX.block(dimensions_cum(on->dim_cum_index_y),0,on->rank,nCol).noalias() += on->Ilist(lInteract)->M2L_operator[index_lInteract_to_iSelf]*X.block(dimensions_cum(on->Ilist(lInteract)->dim_cum_index_y),0,on->Ilist(lInteract)->rank,nCol);
          // stop_oper = clock();
          // Time_oper += double(stop_oper-start_oper)/double(CLOCKS_PER_SEC);
        }

        // start_oper = clock();
        // AX.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,nCol).noalias() += -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol);
        AX.block(dimensions_cum(on->dim_cum_index_y),0,on->rank,nCol).noalias() += -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(on->dim_cum_index_z),0,on->rank,nCol);
        // stop_oper = clock();
        // Time_oper += double(stop_oper-start_oper)/double(CLOCKS_PER_SEC);
      }

    }

  }

  // outfile << "Time_index: " << Time_index << endl;
  // outfile << "Time_oper: " << Time_oper << endl;

}


void IFMM_Matrix::A_extended_transpose_multiplication(MatrixXd &X, MatrixXd &AX, ivec &map_index_to_id, ivec &dimensions_cum, int & leafLevel){

  // cout << "Start A_extended_multiplication..." << endl;

  // int nVariables = map_index_to_id.m;
  int nCol = X.cols();

  // LEAF LEVEL
  int iLevel=leafLevel; 
  // printf("iLevel: %d\n",iLevel);
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
    // CURRENT NODE
    oct_node * on = lvl_p(k);


    // int index_m2i_iSelf_x=0;
    // // for (int i=0; i<nVariables; i++){
    // for (int i=nVariables_cum(leafLevel-iLevel); i<nVariables_cum(leafLevel-iLevel+1); i++){
    //   if (map_index_to_id(i)==on->id){
    //     index_m2i_iSelf_x = i;
    //     break;
    //   }
    // }
    // int index_m2i_iSelf_z = index_m2i_iSelf_x+1;
    // int index_m2i_iSelf_y=0;
    // // for (int i=0; i<nVariables; i++){
    // for (int i=nVariables_cum(leafLevel-iLevel+1); i<nVariables_cum(leafLevel-iLevel+2); i++){
    //   if (map_index_to_id(i)==on->id){
    //     index_m2i_iSelf_y = i;
    //   }
    // }

    // EQUATION XI
    for (int iNgbr = 0; iNgbr < on->ngbr.m; ++iNgbr) {
      int index_iNgbr_to_iSelf=0;
      for (int i=0; i < on->ngbr(iNgbr)->ngbr.m;i++){
        if (on->ngbr(iNgbr)->ngbr(i)->id == on->id)
        {
          index_iNgbr_to_iSelf=i; // index to current node (in the neighbour list of i)
          break;
        }
      }
      // int index_m2i_iNgbr_to_iSelf=0;
      // // for (int i=0; i<nVariables; i++){
      // for (int i=nVariables_cum(leafLevel-iLevel); i<nVariables_cum(leafLevel-iLevel+1); i++){
      //   if (map_index_to_id(i)==on->ngbr(iNgbr)->id){
      //     index_m2i_iNgbr_to_iSelf = i;
      //     break;
      //   }
      // }
      // AX.block(dimensions_cum(index_m2i_iSelf_x),0,on->nDOF,nCol).noalias() += on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf]*X.block(dimensions_cum(index_m2i_iNgbr_to_iSelf),0,on->ngbr(iNgbr)->nDOF,nCol);
      // AX.block(dimensions_cum(index_m2i_iSelf_x),0,on->nDOF,nCol).noalias() += on->P2P_operator[iNgbr].transpose()*X.block(dimensions_cum(index_m2i_iNgbr_to_iSelf),0,on->ngbr(iNgbr)->nDOF,nCol);
      AX.block(dimensions_cum(on->dim_cum_index_x),0,on->nDOF,nCol).noalias() += on->P2P_operator[iNgbr].transpose()*X.block(dimensions_cum(on->ngbr(iNgbr)->dim_cum_index_x),0,on->ngbr(iNgbr)->nDOF,nCol);
    }
    // AX.block(dimensions_cum(index_m2i_iSelf_x),0,on->nDOF,nCol).noalias() += on->L2P_operator*X.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol);
    // AX.block(dimensions_cum(index_m2i_iSelf_x),0,on->nDOF,nCol).noalias() += on->P2M_operator.transpose()*X.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol);
    AX.block(dimensions_cum(on->dim_cum_index_x),0,on->nDOF,nCol).noalias() += on->P2M_operator.transpose()*X.block(dimensions_cum(on->dim_cum_index_z),0,on->rank,nCol);

    // EQUATION ZI
    // AX.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol).noalias() = -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,nCol) + on->P2M_operator*X.block(dimensions_cum(index_m2i_iSelf_x),0,on->nDOF,nCol);
    // AX.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol).noalias() = -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,nCol) + on->L2P_operator.transpose()*X.block(dimensions_cum(index_m2i_iSelf_x),0,on->nDOF,nCol);
    AX.block(dimensions_cum(on->dim_cum_index_z),0,on->rank,nCol).noalias() = -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(on->dim_cum_index_y),0,on->rank,nCol) + on->L2P_operator.transpose()*X.block(dimensions_cum(on->dim_cum_index_x),0,on->nDOF,nCol);


  
    if (iLevel==2){
      // EQUATION YI
      for (int lInteract=0; lInteract < on->Ilist.m; lInteract++){
        int index_lInteract_to_iSelf=0;
        for (int i=0; i < on->Ilist(lInteract)->Ilist.m;i++){
          if (on->Ilist(lInteract)->Ilist(i)->id == on->id)    
          {
            index_lInteract_to_iSelf=i; // index to current node (in the interaction list)
            break;
          }
        }
        // int index_m2i_lInteract_to_iSelf=0;
        // for (int i=0; i<nVariables; i++){
        //   if (map_index_to_id(i)==on->Ilist(lInteract)->id){
        //     index_m2i_lInteract_to_iSelf = i;
        //   }
        // }
        // AX.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,nCol).noalias() += on->Ilist(lInteract)->M2L_operator[index_lInteract_to_iSelf]*X.block(dimensions_cum(index_m2i_lInteract_to_iSelf),0,on->Ilist(lInteract)->rank,nCol);
        // AX.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,nCol).noalias() += on->M2L_operator[lInteract].transpose()*X.block(dimensions_cum(index_m2i_lInteract_to_iSelf),0,on->Ilist(lInteract)->rank,nCol);
        AX.block(dimensions_cum(on->dim_cum_index_y),0,on->rank,nCol).noalias() += on->M2L_operator[lInteract].transpose()*X.block(dimensions_cum(on->Ilist(lInteract)->dim_cum_index_y),0,on->Ilist(lInteract)->rank,nCol);

      }
      // AX.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,nCol).noalias() += -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol);
      AX.block(dimensions_cum(on->dim_cum_index_y),0,on->rank,nCol).noalias() += -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(on->dim_cum_index_z),0,on->rank,nCol);
    }
  }


  // NON-LEAF LEVELS
  for (int iLevel=leafLevel-1;iLevel>=2;iLevel--){ 
    // printf("iLevel: %d\n",iLevel);
#if defined(_OPENMP)
#pragma omp parallel for
#endif
   for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
      // CURRENT NODE
      oct_node * on = lvl_p(k);

      // printf("k: %d\n",k);

      // int index_m2i_iSelf_z=0;
      // // for (int i=0; i<nVariables; i++){
      // for (int i=nVariables_cum(leafLevel-iLevel); i<nVariables_cum(leafLevel-iLevel+1); i++){
      //   if (map_index_to_id(i)==on->id){
      //     index_m2i_iSelf_z = i;
      //     break;
      //   }
      // }
      // int index_m2i_iSelf_y=0;
      // for (int i=0; i<nVariables; i++){
      //   if (map_index_to_id(i)==on->id){
      //     index_m2i_iSelf_y = i;
      //   }
      // }

      // EQUATION YI OF CHILDREN
      // bool firstChild=false;
      // int index_m2i_iChild_y_firstChild = 0;
      for (int iChild=0;iChild<8;iChild++){
        if (on->child[iChild] != NULL){

          // printf("on->id,on->child[iChild]->id: %d %d\n",on->id,on->child[iChild]->id);
          // int index_m2i_iChild_z=0;
          // // for (int i=0; i<nVariables; i++){
          // //   if (map_index_to_id(i)==on->child[iChild]->id){
          // //     if (iLevel+1==l){
          // //       index_m2i_iChild_z = i+1;
          // //       break;
          // //     }
          // //     else{
          // //       index_m2i_iChild_z = i;
          // //     }
          // //   }
          // // }
          // if (iLevel+1==l){
          //   for (int i=nVariables_cum(leafLevel-iLevel-1); i<nVariables_cum(leafLevel-iLevel); i++){
          //     if (map_index_to_id(i)==on->child[iChild]->id){
          //       index_m2i_iChild_z = i+1;
          //       break;
          //     }
          //   }
          // }
          // else{
          //   for (int i=nVariables_cum(leafLevel-iLevel); i<nVariables_cum(leafLevel-iLevel+1); i++){
          //      if (map_index_to_id(i)==on->child[iChild]->id){
          //         index_m2i_iChild_z = i;
          //         break;
          //      }
          //   }
          // }
          // int index_m2i_iChild_y=0;
          // // for (int i=0; i<nVariables; i++){
          // for (int i=nVariables_cum(leafLevel-iLevel); i<nVariables_cum(leafLevel-iLevel+1); i++){
          //   if (map_index_to_id(i)==on->child[iChild]->id){
          //     index_m2i_iChild_y = i;
          //     break;
          //   }
          // }

          // if (firstChild==false){
            // firstChild=true;
            // index_m2i_iChild_y_firstChild = index_m2i_iChild_y;
          // }

          for (int lInteract=0; lInteract < on->child[iChild]->Ilist.m; lInteract++){
            int index_lInteract_to_iChild=0;
            for (int i=0; i < on->child[iChild]->Ilist(lInteract)->Ilist.m;i++){
              if (on->child[iChild]->Ilist(lInteract)->Ilist(i)->id == on->child[iChild]->id)    
              {
                index_lInteract_to_iChild=i; // index to current node (in the interaction list)
                break;
              }
            }
            // int index_m2i_lInteract_to_iChild=0;
            // // for (int i=0; i<nVariables; i++){
            // for (int i=nVariables_cum(leafLevel-iLevel); i<nVariables_cum(leafLevel-iLevel+1); i++){
            //   if (map_index_to_id(i)==on->child[iChild]->Ilist(lInteract)->id){
            //     index_m2i_lInteract_to_iChild = i;
            //     break;
            //   }
            // }
            // AX.block(dimensions_cum(index_m2i_iChild_y),0,on->child[iChild]->rank,nCol).noalias() += on->child[iChild]->Ilist(lInteract)->M2L_operator[index_lInteract_to_iChild]*X.block(dimensions_cum(index_m2i_lInteract_to_iChild),0,on->child[iChild]->Ilist(lInteract)->rank,nCol);
            // AX.block(dimensions_cum(index_m2i_iChild_y),0,on->child[iChild]->rank,nCol).noalias() += on->child[iChild]->M2L_operator[lInteract].transpose()*X.block(dimensions_cum(index_m2i_lInteract_to_iChild),0,on->child[iChild]->Ilist(lInteract)->rank,nCol);
             AX.block(dimensions_cum(on->child[iChild]->dim_cum_index_y),0,on->child[iChild]->rank,nCol).noalias() += on->child[iChild]->M2L_operator[lInteract].transpose()*X.block(dimensions_cum(on->child[iChild]->Ilist(lInteract)->dim_cum_index_y),0,on->child[iChild]->Ilist(lInteract)->rank,nCol);
          }
          // AX.block(dimensions_cum(index_m2i_iChild_y),0,on->child[iChild]->rank,nCol).noalias() +=  -MatrixXd::Identity(on->child[iChild]->rank,on->child[iChild]->rank)*X.block(dimensions_cum(index_m2i_iChild_z),0,on->child[iChild]->rank,nCol);
          AX.block(dimensions_cum(on->child[iChild]->dim_cum_index_y),0,on->child[iChild]->rank,nCol).noalias() +=  -MatrixXd::Identity(on->child[iChild]->rank,on->child[iChild]->rank)*X.block(dimensions_cum(on->child[iChild]->dim_cum_index_z),0,on->child[iChild]->rank,nCol);


          // int index_m2i_parent_to_iSelf=0;
          // for (int i=0; i<nVariables; i++){
          //   if (map_index_to_id(i)==on->parent->id){
          //     index_m2i_parent_to_iSelf = i;
          //     break;
          //   }
          // }
          // cout << "on->child[iChild]->rank" << on->child[iChild]->rank << endl;
          // cout << "on->L2L_operator[iChild].rows()" << on->L2L_operator[iChild].rows() << endl;
          // cout << "on->L2L_operator[iChild].cols()" << on->L2L_operator[iChild].cols() << endl;
          // cout << "on->rank" << on->rank << endl;

          // AX.block(dimensions_cum(index_m2i_iChild_y),0,on->child[iChild]->rank,nCol).noalias() += on->L2L_operator[iChild]*X.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol);
          // AX.block(dimensions_cum(index_m2i_iChild_y),0,on->child[iChild]->rank,nCol).noalias() += on->M2M_operator[iChild].transpose()*X.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol);
          AX.block(dimensions_cum(on->child[iChild]->dim_cum_index_y),0,on->child[iChild]->rank,nCol).noalias() += on->M2M_operator[iChild].transpose()*X.block(dimensions_cum(on->dim_cum_index_z),0,on->rank,nCol);
          // MatrixXd temp = on->L2P_operator*X.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol);
          // AX.block(dimensions_cum(index_m2i_iChild_y),0,on->child[iChild]->rank,nCol).noalias() += temp.block(dimensions_cum(index_m2i_iChild_y)-dimensions_cum(index_m2i_iChild_y_firstChild),0,on->child[iChild]->rank,nCol);
        }
      }

      // EQUATION ZI
      // AX.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol).noalias() = -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,nCol);
      AX.block(dimensions_cum(on->dim_cum_index_z),0,on->rank,nCol).noalias() = -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(on->dim_cum_index_y),0,on->rank,nCol);

      
      // int dimension_local=0;
      // for (int iChild=0; iChild<8; iChild++){
      //   if (on->child[iChild] != NULL){
      //     int index_m2i_iChild_to_iSelf_y=0;
      //     for (int i=0; i<nVariables; i++){
      //       if (map_index_to_id(i)==on->child[iChild]->id){
      //         index_m2i_iChild_to_iSelf_y = i;
      //         }
      //     }
      //     dimension_local+=on->child[iChild]->rank;
      //   }
      // }
      // MatrixXd temp = MatrixXd::Zero(dimension_local,nCol);
      for (int iChild=0; iChild<8; iChild++){
        if (on->child[iChild] != NULL){
          // int index_m2i_iChild_to_iSelf_y=0;
          // for (int i=0; i<nVariables; i++){
          // for (int i=nVariables_cum(leafLevel-iLevel); i<nVariables_cum(leafLevel-iLevel+1); i++){
          //   if (map_index_to_id(i)==on->child[iChild]->id){
          //     index_m2i_iChild_to_iSelf_y = i;
          //     break;
          //     }
          // }
          // AX.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol).noalias() += on->M2M_operator[iChild]*X.block(dimensions_cum(index_m2i_iChild_to_iSelf_y),0,on->child[iChild]->rank,nCol);
          // AX.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol).noalias() += on->L2L_operator[iChild].transpose()*X.block(dimensions_cum(index_m2i_iChild_to_iSelf_y),0,on->child[iChild]->rank,nCol);
          AX.block(dimensions_cum(on->dim_cum_index_z),0,on->rank,nCol).noalias() += on->L2L_operator[iChild].transpose()*X.block(dimensions_cum(on->child[iChild]->dim_cum_index_y),0,on->child[iChild]->rank,nCol);
          // temp.block(dimensions_cum(index_m2i_iChild_to_iSelf_y)-dimensions_cum(index_m2i_iChild_y_firstChild),0,on->child[iChild]->rank,nCol) = X.block(dimensions_cum(index_m2i_iChild_to_iSelf_y),0,on->child[iChild]->rank,nCol);
        }
      }
      // AX.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol).noalias() += on->P2M_operator*temp;


      if (iLevel==2){
        // EQUATION YI
        for (int lInteract=0; lInteract < on->Ilist.m; lInteract++){
          int index_lInteract_to_iSelf=0;
          for (int i=0; i < on->Ilist(lInteract)->Ilist.m;i++){
            if (on->Ilist(lInteract)->Ilist(i)->id == on->id)    
            {
              index_lInteract_to_iSelf=i; // index to current node (in the interaction list)
              break;
            }
          }
          // int index_m2i_lInteract_to_iSelf=0;
          // // for (int i=0; i<nVariables; i++){
          // for (int i=nVariables_cum(leafLevel-1); i<nVariables_cum(leafLevel); i++){
          //   if (map_index_to_id(i)==on->Ilist(lInteract)->id){
          //     index_m2i_lInteract_to_iSelf = i;
          //     break;
          //   }
          // }
          // AX.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,nCol).noalias() += on->Ilist(lInteract)->M2L_operator[index_lInteract_to_iSelf]*X.block(dimensions_cum(index_m2i_lInteract_to_iSelf),0,on->Ilist(lInteract)->rank,nCol);
          // AX.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,nCol).noalias() += on->M2L_operator[lInteract].transpose()*X.block(dimensions_cum(index_m2i_lInteract_to_iSelf),0,on->Ilist(lInteract)->rank,nCol);
          AX.block(dimensions_cum(on->dim_cum_index_y),0,on->rank,nCol).noalias() += on->M2L_operator[lInteract].transpose()*X.block(dimensions_cum(on->Ilist(lInteract)->dim_cum_index_y),0,on->Ilist(lInteract)->rank,nCol);

        }
        // AX.block(dimensions_cum(index_m2i_iSelf_y),0,on->rank,nCol).noalias() += -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(index_m2i_iSelf_z),0,on->rank,nCol);
        AX.block(dimensions_cum(on->dim_cum_index_y),0,on->rank,nCol).noalias() += -MatrixXd::Identity(on->rank,on->rank)*X.block(dimensions_cum(on->dim_cum_index_z),0,on->rank,nCol);

      }
    }
  }
}

void IFMM_Matrix::elimination(ofstream& outfile){
  
  //20200116  clock_t start_transfer, stop_transfer; double Time_transfer=0;
  //20200116  clock_t start_estimate_lsv, stop_estimate_lsv; double Time_estimate_lsv=0;

  double sigma0_Aextended;

  //20200116  start_estimate_lsv = clock();
  TICK(timer_estimate_lsv);
  estimate_lsv(sigma0_Aextended,l,outfile);
  fprintf(stderr, "sigma0_Aextended: %f\n",sigma0_Aextended);
  //20200116  stop_estimate_lsv = clock();
  //20200116  Time_estimate_lsv += double(stop_estimate_lsv-start_estimate_lsv)/double(CLOCKS_PER_SEC);
  //20200116  outfile << "estimate_lsv: " << Time_estimate_lsv << endl;
  TACK(timer_estimate_lsv, outfile);

  reuse=false;

  for (int iLevel=l;iLevel>=2;iLevel--){
    printf("Level %d\n",iLevel );
    outfile << "Level: " << iLevel << endl;


    bool IsLeaf = (iLevel==l) ? true:false;

    //20200116    clock_t start_elim_level, stop_elim_level;
    //20200116    double time_elim_level=0;
    TICK(timer_elim_level);

    // Loop over all levels
    //20200116    start_elim_level =clock();
    eliminate_level(sigma0_Aextended,iLevel,IsLeaf,outfile);
    //20200116    stop_elim_level =clock();

    //20200116    time_elim_level += double(stop_elim_level-start_elim_level)/double(CLOCKS_PER_SEC);
    //20200116    printf("time_elim_level: %f\n",time_elim_level);
    //20200116    outfile << "time_elim_level: " << time_elim_level << endl;
    TACK(timer_elim_level, outfile);



    //20200116    start_transfer=clock();
    TICK(timer_transfer);
    // PREPARE PARENT LEVEL
    if (iLevel>2){
      bool isl_p=false;
      int lvl_local=iLevel-1;
      for (int k = lvl_idx(iLevel-1); k < lvl_idx(iLevel ); ++k) {
        oct_node * on = lvl_p(k);
        
        // SET nDOF (BASED ON RANK OF CHILDREN)
        set_nDOF(on,n,isl_p);
      }
#if defined(_OPENMP) && defined(PARA_CHECK_RANK) // then, perform check_rank() separately.
#pragma omp parallel for
#endif
      for (int k = lvl_idx(iLevel-1); k < lvl_idx(iLevel ); ++k) {
        oct_node * on = lvl_p(k);

        // set P2L_operator
        set_P2L_operator(on,n,isl_p);
    
        // set M2P_operator
        set_M2P_operator(on,n,isl_p);

        // transfer RHS
        transfer_RHS(on);

        // set M2L of level as P2P of higher level
        transfer_M2L_to_P2P(on);
        transfer_M2M_to_P2M(on);
        transfer_L2L_to_L2P(on);

        // check rank
#if !defined(_OPENMP) || !defined(PARA_CHECK_RANK)
        check_rank(on,lvl_local,epsilon_rel_basis);
#endif
      }

#if defined(_OPENMP) && defined(PARA_CHECK_RANK)
#pragma omp parallel for
      for (int k = lvl_idx(iLevel-1); k < lvl_idx(iLevel ); ++k) {
        oct_node * on = lvl_p(k);
        check_rank_self_setup(on, epsilon_rel_basis);
      }
#pragma omp parallel for
      for (int k = lvl_idx(iLevel-1); k < lvl_idx(iLevel ); ++k) {
        oct_node * on = lvl_p(k);
        check_rank_self(on, lvl_local);
      }
#endif


      for (int k = lvl_idx(iLevel-1); k < lvl_idx(iLevel ); ++k) {   
        oct_node * on = lvl_p(k);   
        get_weights(on,lvl_local,epsilon_rel_basis);    
      }
      
      for (int k = lvl_idx(iLevel-1); k < lvl_idx(iLevel ); ++k) {
        oct_node * on = lvl_p(k);
        on->L2L_operator.clear();
        on->M2M_operator.clear();
        for (int iChild=0; iChild<8; iChild++){
          if (on->child[iChild] != NULL){
            on->child[iChild]->M2L_operator.clear();
          }
        }
      }
      //20200116      start_estimate_lsv = clock();
      estimate_lsv(sigma0_Aextended,lvl_local,outfile);
      //20200116      stop_estimate_lsv = clock();
      //20200116      Time_estimate_lsv += double(stop_estimate_lsv-start_estimate_lsv)/double(CLOCKS_PER_SEC);
      //20200116      outfile << "estimate_lsv: " << Time_estimate_lsv << endl;

    }
    //20200116    stop_transfer=clock();
    //20200116    Time_transfer+= double(stop_transfer-start_transfer)/double(CLOCKS_PER_SEC);
    TACK(timer_transfer, outfile);
  }

  //20200116  printf("Time_transfer: %f\n",Time_transfer);
  //20200116  outfile << "transfer time: " << Time_transfer << endl;
  //20200116  outfile << "estimate_lsv: " << Time_estimate_lsv << endl;
}


void didIdoJ(oct_node *on, int iSelf, int iNgbr, int jNgbr, int index_iNgbr_to_jNgbr, int index_jNgbr_to_iNgbr, MatrixXd &rho_self, double &epsilon_rel_fillin, double &epsilon_rel_basis, double &sigma0_Aextended, double machine_eps, int rank_RSVD_min, MatrixXd &P, int LR_mode, bool Delay_update, ofstream &outfile, int curr_level, bool WellSeparated, MatrixXd &alpha_self, MatrixXd &alpha_self_inv, VectorXd &kappa_s, VectorXd &iota_s, vector<MatrixXd> &beta_is, vector<MatrixXd> &xi_si, vector<MatrixXd> &psi_si, vector<MatrixXd> &mu_is, vector<MatrixXd> &zeta_si, vector<MatrixXd> &eta_si, vector<MatrixXd> &p2lto, vector<MatrixXd> &p2pto)
{
  // cout << "ij: 1 0"<< endl;

  // // MatrixXd omicron_sj = P2P_operator_self_LU.solve(on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iSelf]);

  if (iNgbr==iSelf){

    // MatrixXd nu_sj = - on->P2M_operator*(on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iSelf]));
    // MatrixXd mu_js = - on->P2P_operator[jNgbr]*(on->P2P_operator[iSelf].partialPivLu().solve(on->L2P_operator));
    // MatrixXd chi_jj = - on->P2P_operator[jNgbr]*(on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iSelf]));

    // MatrixXd nu_sj = - on->P2M_operator*(P2P_operator_self_LU.solve(on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iSelf]));
    // // MatrixXd nu_sj = - on->P2M_operator*omicron_sj;
    // MatrixXd mu_js = - on->P2P_operator[jNgbr]*(P2P_operator_self_LU.solve(on->L2P_operator));
    // //  MatrixXd mu_js = - on->P2P_operator[jNgbr]*rho_self;
    // MatrixXd chi_jj = - on->P2P_operator[jNgbr]*(P2P_operator_self_LU.solve(on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iSelf]));
    // // MatrixXd chi_jj = - on->P2P_operator[jNgbr]*omicron_sj;


    if (on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_jNgbr].cols() != on->ngbr(iNgbr)->rank){
      on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_jNgbr] = MatrixXd::Zero(on->ngbr(jNgbr)->nDOF,on->ngbr(iNgbr)->rank);
      on->ngbr(jNgbr)->P2L_operator[index_jNgbr_to_iNgbr] = MatrixXd::Zero(on->ngbr(iNgbr)->rank,on->ngbr(jNgbr)->nDOF);
    }

    //20200108	    start_operations=clock();
    // M2P_ji
    // on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_jNgbr] += mu_js*alpha_self.partialPivLu().solve(MatrixXd::Identity(on->rank,on->rank));
    // on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_jNgbr] += mu_js*alpha_self_LU.solve(MatrixXd::Identity(on->rank,on->rank));
    // on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_jNgbr] += mu_js*alpha_self_inv;
    on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_jNgbr].noalias() += mu_is[jNgbr]*alpha_self_inv;


    // P2L_ij
    // on->ngbr(jNgbr)->P2L_operator[index_jNgbr_to_iNgbr] += alpha_self.partialPivLu().solve(nu_sj);
    // on->ngbr(jNgbr)->P2L_operator[index_jNgbr_to_iNgbr] += alpha_self_LU.solve(nu_sj);
    on->ngbr(jNgbr)->P2L_operator[index_jNgbr_to_iNgbr] += zeta_si[jNgbr];

    //20200108	    stop_operations=clock();
    //20200108	    Time_operations+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);

    //20200108	    number_operations_f++;
    //20200108	    Time_operations_f+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);

  }
  else{

    // MatrixXd nu_sj = - on->P2M_operator*(on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iSelf]));
    // MatrixXd mu_js = - on->P2P_operator[jNgbr]*(on->P2P_operator[iSelf].partialPivLu().solve(on->L2P_operator));
              

    // MatrixXd nu_sj = - on->P2M_operator*(P2P_operator_self_LU.solve(on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iSelf]));
    // // MatrixXd nu_sj = - on->P2M_operator*omicron_sj;
    // MatrixXd mu_js = - on->P2P_operator[jNgbr]*(P2P_operator_self_LU.solve(on->L2P_operator));
    // // MatrixXd mu_js = - on->P2P_operator[jNgbr]*rho_self;


    // M2P_ji
    // MatrixXd M2P_fillin_ji = - on->P2P_operator[jNgbr]*on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_iSelf]) - mu_js*alpha_self.partialPivLu().solve(gamma_si);
    // MatrixXd M2P_fillin_ji = - on->P2P_operator[jNgbr]*P2P_operator_self_LU.solve(on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_iSelf]) - mu_js*alpha_self_LU.solve(gamma_si);
    // MatrixXd M2P_fillin_ji = - on->P2P_operator[jNgbr]*tau_si - mu_js*alpha_self_LU.solve(gamma_si);
    // MatrixXd M2P_fillin_ji = - on->P2P_operator[jNgbr]*tau_si[iNgbr] - mu_is[jNgbr]*xi_si[iNgbr];
    // TEST MatrixXd M2P_fillin_ji = - on->P2P_operator[jNgbr]*psi_si[iNgbr];


    // P2L_ij
    // MatrixXd P2L_fillin_ij = - on->P2L_operator[iNgbr]*on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iSelf])- beta_is*alpha_self.partialPivLu().solve(nu_sj);
    // MatrixXd P2L_fillin_ij = - on->P2L_operator[iNgbr]*P2P_operator_self_LU.solve(on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iSelf]) - beta_is*alpha_self_LU.solve(nu_sj);
    // MatrixXd P2L_fillin_ij = - on->P2L_operator[iNgbr]*omicron_sj - beta_is*alpha_self_LU.solve(nu_sj);
    // MatrixXd P2L_fillin_ij = - on->P2L_operator[iNgbr]*omicron_si[jNgbr] - beta_is[iNgbr]*zeta_si[jNgbr];
    // TEST MatrixXd P2L_fillin_ij = - on->P2L_operator[iNgbr]*eta_si[jNgbr];


    // if (on->ngbr(jNgbr)->id==177 || on->ngbr(jNgbr)->id==174 || on->ngbr(jNgbr)->id==272){
    //   cout << "M2P_fillin_ji" << M2P_fillin_ji << endl;
    //   cout << "P2L_fillin_ij" << P2L_fillin_ij << endl;

    //   cout << "on->P2P_operator[jNgbr]" << endl << on->P2P_operator[jNgbr] << endl;
    //   cout << "tau_si[iNgbr]" << endl <<tau_si[iNgbr] << endl;
    //   cout << "mu_is[jNgbr]" << endl <<mu_is[jNgbr]<< endl;
    //   cout << "xi_si[iNgbr]" << endl <<xi_si[iNgbr]<< endl;

    //   cout << "alpha_self" << endl << alpha_self << endl;
    //   cout << "gamma_si[iNgbr]" << endl << gamma_si[iNgbr] << endl;

    //   cout << "on->P2M_operator.rows()" << on->P2M_operator.rows() << endl;
    //   cout << "on->P2M_operator.cols()" << on->P2M_operator.cols() << endl;

    //   cout << "on->L2P_operator.rows()" << on->L2P_operator.rows() << endl;
    //   cout << "on->L2P_operator.cols()" << on->L2P_operator.cols() << endl;

    //   cout << "on->L2P_operator" << on->L2P_operator << endl;
    //   cout << "on->P2M_operator" << on->P2M_operator << endl;
    // }


    // if (isnan(M2P_fillin_ji.norm())){
    //   cout << "on->P2P_operator[jNgbr]" << endl << on->P2P_operator[jNgbr] << endl; 
    //   cout << "tau_si[iNgbr]" << endl << tau_si[iNgbr] << endl;
    //   cout << "mu_is[jNgbr]" << endl << mu_is[jNgbr] << endl;
    //   cout << "xi_si[iNgbr]" << endl << xi_si[iNgbr] << endl;
    // }



    if (WellSeparated==false){

      //20200108	      start_operations=clock();
      // M2P_ji
      // on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_jNgbr]+=M2P_fillin_ji;
      // // on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_jNgbr].noalias()+=- on->P2P_operator[jNgbr]*psi_si[iNgbr];

      // P2L_ij
      // on->ngbr(jNgbr)->P2L_operator[index_jNgbr_to_iNgbr]+=P2L_fillin_ij;
      // // on->ngbr(jNgbr)->P2L_operator[index_jNgbr_to_iNgbr].noalias()+=- on->P2L_operator[iNgbr]*eta_si[jNgbr];

      //20200108      MatrixXd temp1 = - on->P2P_operator[jNgbr]*psi_si[iNgbr];
      //20200108      // MatrixXd temp2 = - on->P2L_operator[iNgbr]*eta_si[jNgbr];
      //20200108
      //20200108      on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_jNgbr].noalias()+= temp1; 
      //20200108      if (curr_level==l){
      //20200108	on->ngbr(jNgbr)->P2L_operator[index_jNgbr_to_iNgbr].noalias()+= temp1.transpose(); // USE SYMMETRY
      //20200108      }
      //20200108      else{
      //20200108	on->ngbr(jNgbr)->P2L_operator[index_jNgbr_to_iNgbr].noalias()+= - on->P2L_operator[iNgbr]*eta_si[jNgbr];
      //20200108      }

      on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_jNgbr].noalias()+= - on->P2P_operator[jNgbr]*psi_si[iNgbr];
      on->ngbr(jNgbr)->P2L_operator[index_jNgbr_to_iNgbr].noalias()+= - on->P2L_operator[iNgbr]*eta_si[jNgbr];

      //20200108	      stop_operations=clock();
      //20200108	      Time_operations+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);

      //20200108	      number_operations_g_1++;
      //20200108	      Time_operations_g_1+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);


    }
    else{
      /////////////////////////////////////////////////////////////////
      // iNGBR and jNGBR ARE WELL SEPARATED -> COMPRESS AND REDIRECT //
      /////////////////////////////////////////////////////////////////
      // timeType t0, t00, t1, t2, t3,t4,t5,t6,t7,t8; // Time variables

      // cout << "Compress and redirect" << endl;

      // MatrixXd M2P_fillin_ji = - on->P2P_operator[jNgbr]*psi_si[iNgbr];
      // MatrixXd P2L_fillin_ij = - on->P2L_operator[iNgbr]*eta_si[jNgbr];

      //20200108	      start_operations=clock();

      MatrixXd M2P_fillin_ji;
      MatrixXd P2L_fillin_ij;

      if (LR_mode !=3){
	M2P_fillin_ji = - on->P2P_operator[jNgbr]*psi_si[iNgbr];
	P2L_fillin_ij = - on->P2L_operator[iNgbr]*eta_si[jNgbr];
      }


      //20200108	      stop_operations=clock();
      //20200108	      Time_operations+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);

      //20200108	      number_operations_g_2++;
      //20200108	      Time_operations_g_2+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);



      int jSV_M2P=0;
      int jSV_P2L=0;
      int SV_U_comb_j=0;
      int SV_V_comb_j=0;

      MatrixXd U_j_prime, K_ji_prime;

      MatrixXd K_ij_prime, V_j_prime;
      MatrixXd U_comb_j, V_comb_j;

      MatrixXd Sigma_ji_prime, Sigma_ij_prime;
                    
      MatrixXd U_recomp_j, R_recomp_U_j;
      MatrixXd V_recomp_j, R_recomp_V_j;

      MatrixXd Sigma_diag_U_comb_j, Sigma_diag_V_comb_j;

      bool NeedUpdate=true;


      if (LR_mode==0){// SVD
#include "LR0-1.cpp"
      }
      if(LR_mode==1){ // ACA
#include "LR1-1.cpp"
      }

      if (LR_mode==2){ // RSVD
#include "LR2-1.cpp"
      }
      if (LR_mode==3){ // RSVD - optimized
#include "LR3-1.cpp"
      }
                

      ///////////////////////////  
      // UPDATE EDGES OF jNgbr //
      ///////////////////////////  

      if (NeedUpdate==true){

	//20200108		start_updates =clock();

	MatrixXd Sigma_L2P_j_inverse=MatrixXd::Zero(on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->rank);
	MatrixXd Sigma_P2M_j_inverse=MatrixXd::Zero(on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->rank);
	for (int jRank=0; jRank< on->ngbr(jNgbr)->rank; jRank++){
	  Sigma_L2P_j_inverse(jRank,jRank)=1.0/on->ngbr(jNgbr)->Sigma_L2P(jRank,jRank);
	  Sigma_P2M_j_inverse(jRank,jRank)=1.0/on->ngbr(jNgbr)->Sigma_P2M(jRank,jRank);
	}


	MatrixXd Sigma_ij_inverse=MatrixXd::Zero(jSV_M2P,jSV_M2P);
	MatrixXd Sigma_ji_inverse=MatrixXd::Zero(jSV_M2P,jSV_M2P);
	for (int index_rank=0; index_rank<jSV_M2P; index_rank++){
	  Sigma_ij_inverse(index_rank,index_rank)=1.0/Sigma_ij_prime(index_rank,index_rank);
	  Sigma_ji_inverse(index_rank,index_rank)=1.0/Sigma_ji_prime(index_rank,index_rank);
	}



	MatrixXd R_U_j = R_recomp_U_j.block(0,0,SV_U_comb_j,on->ngbr(jNgbr)->rank)*Sigma_L2P_j_inverse;
	MatrixXd R_U_j_prime = R_recomp_U_j.block(0,on->ngbr(jNgbr)->rank,SV_U_comb_j,jSV_M2P)*Sigma_ji_inverse;
                
	MatrixXd R_V_j_T = Sigma_P2M_j_inverse*R_recomp_V_j.block(0,0,SV_V_comb_j,on->ngbr(jNgbr)->rank).transpose();
	MatrixXd R_V_j_prime_T = Sigma_ij_inverse*R_recomp_V_j.block(0,on->ngbr(jNgbr)->rank,SV_V_comb_j,jSV_P2L).transpose();



	on->ngbr(jNgbr)->P2M_operator=V_recomp_j.transpose();
	on->ngbr(jNgbr)->L2P_operator=U_recomp_j;

	on->ngbr(jNgbr)->Sigma_P2M = Sigma_diag_V_comb_j;
	on->ngbr(jNgbr)->Sigma_L2P = Sigma_diag_U_comb_j;

                

#if defined(_OPENMP) && defined(IFMM_PARALLELIZE)
#pragma omp critical
	{
#endif
	  
	  // NEIGHBOUR LIST
	  for (int lNgbr=0; lNgbr < on->ngbr(jNgbr)->ngbr.m; lNgbr++){
	    // on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+lNgbr] = on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+lNgbr]*R_recomp_V_j.block(0,0,SV_V_comb_j,on->ngbr(jNgbr)->rank).transpose();
	    on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+lNgbr] = on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+lNgbr]*R_V_j_T;
	    
	    
	    // also reverse
	    int index_lNgbr_to_jNgbr=0;
	    for (int i=0; i < on->ngbr(jNgbr)->ngbr(lNgbr)->ngbr.m;i++){
	      if (on->ngbr(jNgbr)->ngbr(lNgbr)->ngbr(i)->id == on->ngbr(jNgbr)->id)    
		{
		  index_lNgbr_to_jNgbr=i; // index l to i
		  break;
		}
	    }
	    // on->ngbr(jNgbr)->ngbr(lNgbr)->M2L_operator[on->ngbr(jNgbr)->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_jNgbr] = R_recomp_U_j.block(0,0,SV_U_comb_j,on->ngbr(jNgbr)->rank)*on->ngbr(jNgbr)->ngbr(lNgbr)->M2L_operator[on->ngbr(jNgbr)->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_jNgbr];                  // cout << "on->ngbr(jNgbr)->ngbr(lNgbr)->M2L_operator[on->ngbr(jNgbr)->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_jNgbr] " << endl << on->ngbr(jNgbr)->ngbr(lNgbr)->M2L_operator[on->ngbr(jNgbr)->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_jNgbr]  << endl;
	    on->ngbr(jNgbr)->ngbr(lNgbr)->M2L_operator[on->ngbr(jNgbr)->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_jNgbr] = R_U_j*on->ngbr(jNgbr)->ngbr(lNgbr)->M2L_operator[on->ngbr(jNgbr)->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_jNgbr];       
	    
	  }
	  
	  // on->ngbr(iNgbr)->M2L_operator[index_iNgbr_to_jNgbr] = R_recomp_U_j.block(0,0,SV_U_comb_j,on->ngbr(jNgbr)->rank)*on->ngbr(iNgbr)->M2L_operator[index_iNgbr_to_jNgbr]+R_recomp_U_j.block(0,on->ngbr(jNgbr)->rank,SV_U_comb_j,jSV_M2P)*K_ji_prime; 
	  on->ngbr(iNgbr)->M2L_operator[index_iNgbr_to_jNgbr] = R_U_j*on->ngbr(iNgbr)->M2L_operator[index_iNgbr_to_jNgbr]+R_U_j_prime*K_ji_prime; 
	  
	  
          
	  // INTERACTION LIST
	  for (int lInteract=0; lInteract < on->ngbr(jNgbr)->Ilist.m; lInteract++){
	    if(lInteract == index_jNgbr_to_iNgbr){
	      // on->ngbr(jNgbr)->M2L_operator[index_jNgbr_to_iNgbr] = on->ngbr(jNgbr)->M2L_operator[index_jNgbr_to_iNgbr]*R_recomp_V_j.block(0,0,SV_V_comb_j,on->ngbr(jNgbr)->rank).transpose()+K_ij_prime*R_recomp_V_j.block(0,on->ngbr(jNgbr)->rank,SV_V_comb_j,jSV_P2L).transpose(); 
	      on->ngbr(jNgbr)->M2L_operator[index_jNgbr_to_iNgbr] = on->ngbr(jNgbr)->M2L_operator[index_jNgbr_to_iNgbr]*R_V_j_T+K_ij_prime*R_V_j_prime_T; 
	      
	    }
	    else{
	      // on->ngbr(jNgbr)->M2L_operator[lInteract] =  on->ngbr(jNgbr)->M2L_operator[lInteract]*R_recomp_V_j.block(0,0,SV_V_comb_j,on->ngbr(jNgbr)->rank).transpose();
	      on->ngbr(jNgbr)->M2L_operator[lInteract] =  on->ngbr(jNgbr)->M2L_operator[lInteract]*R_V_j_T;
	      
	      
	      // also reverse
	      int index_lNgbr_to_jNgbr=0;
	      for (int i=0; i < on->ngbr(jNgbr)->Ilist(lInteract)->Ilist.m;i++){
		if (on->ngbr(jNgbr)->Ilist(lInteract)->Ilist(i)->id == on->ngbr(jNgbr)->id)    
		  {
		    index_lNgbr_to_jNgbr=i; // index l to i
		    break;
		  }
	      }
	      // on->ngbr(jNgbr)->Ilist(lInteract)->M2L_operator[index_lNgbr_to_jNgbr] = R_recomp_U_j.block(0,0,SV_U_comb_j,on->ngbr(jNgbr)->rank)*on->ngbr(jNgbr)->Ilist(lInteract)->M2L_operator[index_lNgbr_to_jNgbr];
	      on->ngbr(jNgbr)->Ilist(lInteract)->M2L_operator[index_lNgbr_to_jNgbr] = R_U_j*on->ngbr(jNgbr)->Ilist(lInteract)->M2L_operator[index_lNgbr_to_jNgbr];
	      
	    }
	  }


#if defined(_OPENMP) && defined(IFMM_PARALLELIZE)
	} // omp critical
#endif	

	if(curr_level>2){
	  // UPDATE PARENT OPERATORS
	  int id_jNgbr_local=0;
	  for (int iChild=0;iChild<8;iChild++){
	    if (on->ngbr(jNgbr)->parent->child[iChild] != NULL){
	      if (on->ngbr(jNgbr)->parent->child[iChild]->id == on->ngbr(jNgbr)->id){
		id_jNgbr_local=iChild;
		break;
	      }
	    }
	  }
	  // on->ngbr(jNgbr)->parent->M2M_operator[id_jNgbr_local] = on->ngbr(jNgbr)->parent->M2M_operator[id_jNgbr_local]*R_recomp_V_j.block(0,0,SV_V_comb_j,on->ngbr(jNgbr)->rank).transpose();
	  // on->ngbr(jNgbr)->parent->L2L_operator[id_jNgbr_local] = R_recomp_U_j.block(0,0,SV_U_comb_j,on->ngbr(jNgbr)->rank)*on->ngbr(jNgbr)->parent->L2L_operator[id_jNgbr_local];

	  on->ngbr(jNgbr)->parent->M2M_operator[id_jNgbr_local] = on->ngbr(jNgbr)->parent->M2M_operator[id_jNgbr_local]*R_V_j_T;
	  on->ngbr(jNgbr)->parent->L2L_operator[id_jNgbr_local] = R_U_j*on->ngbr(jNgbr)->parent->L2L_operator[id_jNgbr_local];

	}

	// UPDATE RANK OF jNGBR
	// printf("Before - on->ngbr(jNgbr)->rank: %d\n",on->ngbr(jNgbr)->rank);
	on->ngbr(jNgbr)->rank=SV_U_comb_j;
	// printf("After - on->ngbr(jNgbr)->rank: %d\n",on->ngbr(jNgbr)->rank);

	//20200108		stop_updates =clock();
	//20200108		Time_updates += double(stop_updates-start_updates)/double(CLOCKS_PER_SEC);

      }


    }
  }

}


#if defined(IFMM_PARALLELIZE) 
#if !defined(IFMM_PARALLELIZE_CONCURRENCY_SEP)
#define IFMM_PARALLELIZE_CONCURRENCY_SEP 3
#endif
#endif


void IFMM_Matrix::eliminate_level_node(double & sigma0_Aextended, int curr_level, ofstream& outfile, const double machine_eps, const int rank_RSVD_min, MatrixXd &P, const int k)
{
  // This function was a part of eliminate_level().

  
  // CURRENT NODE
  oct_node * on = lvl_p(k);
  // printf("Eliminating node... %d\n", on->id);
  
  // timeType t0, t1;
  // t0=Timer();
  
  // ELIMINATE CURRENT NODE
  on->IsEliminated=1;
  
  int iSelf=0; // index to current node (in the neighbour list of the current node)
  for (int iNgbr=0; iNgbr < on-> ngbr.m; iNgbr++){
    if (on->ngbr(iNgbr)->id==on->id){
      iSelf=iNgbr;
      break;
    }
  }
    
    
  //  // SVD OF P2M
  // JacobiSVD<MatrixXd> svd_P2M_operator(on->P2M_operator.transpose(), ComputeThinU | ComputeThinV);
  // VectorXd Sigma_i_P2M = svd_P2M_operator.singularValues();
  // cout << "Sigma_i_P2M" << endl << Sigma_i_P2M << endl;
    

  // printf("on->rank:%d\n",on->rank);
  // printf("on->nDOF:%d\n",on->nDOF);

  // if (on->id==46 || on->id==48 || on->id==97 || on->id==100){
  // cout << "Start pre-computation..." << endl;
  // }


  //20200108  start_precomp=clock();
  // PRECOMPUTE CERTAIN MATRICES HERE
  // LU DECOMPOSITION OF P2P_ss
  // Eigen::PartialPivLU<MatrixXd> P2P_operator_self_LU = on->P2P_operator[iSelf].partialPivLu();
  on->P2P_operator_self_LU = on->P2P_operator[iSelf].partialPivLu();
  MatrixXd rho_self=on->P2P_operator_self_LU.solve(on->L2P_operator);

  // for (int iChild=0;iChild<8;iChild++){
  //    if (on->child[iChild] != NULL){
  //     printf("on->child[iChild]->id: %d\n",on->child[iChild]->id);
  //     }
  // }


  MatrixXd alpha_self=MatrixXd::Zero(on->rank,on->rank);
  // alpha_self = - on->P2M_operator*(on->P2P_operator[iSelf].partialPivLu().solve(on->L2P_operator));
  // alpha_self = - on->P2M_operator*(P2P_operator_self_LU.solve(on->L2P_operator)); // EFF LU
  alpha_self = - on->P2M_operator*rho_self; // EFF LU

  // LU DECOMPOSITION OF alpha_ii
  // Eigen::PartialPivLU<MatrixXd> alpha_self_LU = alpha_self.partialPivLu(); // EFF LU
  on->alpha_self_LU = alpha_self.partialPivLu(); // EFF LU
  MatrixXd alpha_self_inv=on->alpha_self_LU.solve(MatrixXd::Identity(on->rank,on->rank));

#warning Insert patch.cpp here.

  VectorXd kappa_s = on->P2P_operator_self_LU.solve(on->RHS_leaf);
  VectorXd iota_s = on->alpha_self_LU.solve(on->P2M_operator*(kappa_s));

  // vector<MatrixXd> tau_si(on->ngbr.m);
  MatrixXd tau_si;
  vector<MatrixXd> beta_is(on->ngbr.m);
  // vector<MatrixXd> gamma_si(on->ngbr.m);
  MatrixXd gamma_si;
  // vector<MatrixXd> delta_ii(on->ngbr.m);
  vector<MatrixXd> xi_si(on->ngbr.m);
  vector<MatrixXd> psi_si(on->ngbr.m);

  // vector<MatrixXd> omicron_si(on->ngbr.m);
  MatrixXd omicron_si;
  // vector<MatrixXd> nu_si(on->ngbr.m);
  MatrixXd nu_si;
  vector<MatrixXd> mu_is(on->ngbr.m);
  // vector<MatrixXd> chi_ii(on->ngbr.m);
  vector<MatrixXd> zeta_si(on->ngbr.m);
  vector<MatrixXd> eta_si(on->ngbr.m);

  vector<MatrixXd> p2lto(on->ngbr.m);
  vector<MatrixXd> p2pto(on->ngbr.m);


  //////////////////////////////////////////////////////////////////////////
  // START LOOP OVER ALL POSSIBLE PAIRS OF NEIGHBOURS OF THE CURRENT NODE //
  //////////////////////////////////////////////////////////////////////////  
  // LOOP OVER ALL NEIGHBOURS iNgbr
  for (int iNgbr=0; iNgbr < on->ngbr.m; iNgbr++){

    // GET THE RIGHT INDICES
    int index_iNgbr_to_iSelf=0;
    for (int i=0; i < on->ngbr(iNgbr)->ngbr.m;i++){
      if (on->ngbr(iNgbr)->ngbr(i)->id == on-> ngbr(iSelf)->id)
	{
	  index_iNgbr_to_iSelf=i; // index to current node (in the neighbour list of i)
	  break;
	}
    }
    // int index_iNgbr_to_iNgbr=on->ngbr(iNgbr)->ngbr.m;
    // for (int i=0; i < on->ngbr(iNgbr)->ngbr.m;i++){
    //   if (on->ngbr(iNgbr)->ngbr(i)->id == on-> ngbr(iNgbr)->id)
    //   {
    //     index_iNgbr_to_iNgbr=i; // index to i (in the neighbour list of i)
    //     break;
    //   }
    // }

    // PRECOMPUTE CERTAIN MATRICES HERE
    // MatrixXd beta_is = - (on->P2L_operator[iNgbr]*on->P2P_operator[iSelf].partialPivLu().solve(on->L2P_operator));
    // MatrixXd gamma_si = - on->P2M_operator*on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_iSelf]);
    // MatrixXd delta_ii = - on->P2L_operator[iNgbr]*on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_iSelf]);

    // MatrixXd nu_si = - on->P2M_operator*(on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf]));
    // MatrixXd mu_is = - on->P2P_operator[iNgbr]*(on->P2P_operator[iSelf].partialPivLu().solve(on->L2P_operator));
    // MatrixXd chi_ii = - on->P2P_operator[iNgbr]*(on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf]));
      
    // MatrixXd beta_is = - (on->P2L_operator[iNgbr]*P2P_operator_self_LU.solve(on->L2P_operator));
    // MatrixXd beta_is = - (on->P2L_operator[iNgbr]*rho_self);
    // MatrixXd gamma_si = - on->P2M_operator*P2P_operator_self_LU.solve(on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_iSelf]);
    // MatrixXd gamma_si = - on->P2M_operator*tau_si;
    // MatrixXd delta_ii = - on->P2L_operator[iNgbr]*P2P_operator_self_LU.solve(on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_iSelf]);
    // MatrixXd delta_ii = - on->P2L_operator[iNgbr]*tau_si;

    // MatrixXd nu_si = - on->P2M_operator*(P2P_operator_self_LU.solve(on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf]));
    // // MatrixXd mu_is = - on->P2P_operator[iNgbr]*(P2P_operator_self_LU.solve(on->L2P_operator));
    // MatrixXd mu_is = - on->P2P_operator[iNgbr]*rho_self;
    // MatrixXd chi_ii = - on->P2P_operator[iNgbr]*(P2P_operator_self_LU.solve(on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf]));


    if (on->ngbr(iNgbr)->IsEliminated){
      if(iNgbr!=iSelf){
	// tau_si[iNgbr] = P2P_operator_self_LU.solve(on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_iSelf]);
	tau_si = on->P2P_operator_self_LU.solve(on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_iSelf]);
	beta_is[iNgbr] = - on->P2L_operator[iNgbr]*rho_self;
	gamma_si = - on->P2M_operator*tau_si;
	// delta_ii[iNgbr] = - on->P2L_operator[iNgbr]*tau_si;
	xi_si[iNgbr] = on->alpha_self_LU.solve(gamma_si);
	psi_si[iNgbr] = tau_si - rho_self*xi_si[iNgbr];

	MatrixXd O(on->ngbr(iNgbr)->rank,rank_RSVD_min);
	sample_gaussian(O);
	p2lto[iNgbr] = on->P2L_operator[iNgbr].transpose()*O;
      }
    }
    else{
      // omicron_si[iNgbr] = P2P_operator_self_LU.solve(on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf]);
      omicron_si = on->P2P_operator_self_LU.solve(on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf]);
      nu_si = - on->P2M_operator*omicron_si;
      mu_is[iNgbr] = - on->P2P_operator[iNgbr]*rho_self;
      // chi_ii[iNgbr] = - on->P2P_operator[iNgbr]*omicron_si;
      zeta_si[iNgbr] = on->alpha_self_LU.solve(nu_si);
      eta_si[iNgbr] = omicron_si - rho_self*zeta_si[iNgbr];

      MatrixXd O(on->ngbr(iNgbr)->nDOF,rank_RSVD_min);
      sample_gaussian(O);
      p2pto[iNgbr] = on->P2P_operator[iNgbr].transpose()*O;
    }
    
  }

  //20200108  stop_precomp=clock();
  //20200108  Time_precomp+= double(stop_precomp-start_precomp)/double(CLOCKS_PER_SEC);
      

  // MatrixXd tau_si, beta_is, gamma_si, delta_ii, xi_si;
  // MatrixXd omicron_si, nu_si, mu_is, chi_ii, zeta_si;
      
  // if (on->ngbr(iNgbr)->IsEliminated){
  //   tau_si = P2P_operator_self_LU.solve(on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_iSelf]);
        
  //   beta_is = - on->P2L_operator[iNgbr]*rho_self;
  //   gamma_si = - on->P2M_operator*tau_si;
  //   delta_ii = - on->P2L_operator[iNgbr]*tau_si;
        
  //   xi_si=alpha_self_LU.solve(gamma_si);
  // }
  // else{
  //   omicron_si = P2P_operator_self_LU.solve(on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf]);

  //   nu_si = - on->P2M_operator*omicron_si;
  //   mu_is = - on->P2P_operator[iNgbr]*rho_self;
  //   chi_ii = - on->P2P_operator[iNgbr]*omicron_si;

  //   zeta_si=alpha_self_LU.solve(nu_si);
  // }
      
  // LOOP OVER ALL NEIGHBOURS iNgbr
  for (int iNgbr=0; iNgbr < on-> ngbr.m; iNgbr++){
    // printf("on->ngbr(iNgbr)->id: %d \n",on->ngbr(iNgbr)->id);

    // GET THE RIGHT INDICES
    // int index_iNgbr_to_iSelf=0;
    //   for (int i=0; i < on->ngbr(iNgbr)->ngbr.m;i++){
    //     if (on->ngbr(iNgbr)->ngbr(i)->id == on-> ngbr(iSelf)->id)
    //     {
    //       index_iNgbr_to_iSelf=i; // index to current node (in the neighbour list of i)
    //       break;
    //     }
    //   }
    int index_iNgbr_to_iNgbr=on->ngbr(iNgbr)->ngbr.m;
    for (int i=0; i < on->ngbr(iNgbr)->ngbr.m;i++){
      if (on->ngbr(iNgbr)->ngbr(i)->id == on-> ngbr(iNgbr)->id)
	{
	  index_iNgbr_to_iNgbr=i; // index to i (in the neighbour list of i)
	  break;
	}
    }

    // LOOP OVER ALL NEIGHBOURS jNgbr
    for (int jNgbr=iNgbr; jNgbr < on-> ngbr.m; jNgbr++){
      // printf("on->ngbr(iNgbr)->id,on->ngbr(jNgbr)->id %d %d \n",on->ngbr(iNgbr)->id,on->ngbr(jNgbr)->id);
      // printf("on->ngbr(iNgbr)->rank,on->ngbr(jNgbr)->rank %d %d \n",on->ngbr(iNgbr)->rank,on->ngbr(jNgbr)->rank);

      //20200108      start_ij=clock();

      // Get the right indices
      int index_iNgbr_to_jNgbr=on->ngbr(iNgbr)->ngbr.m;
      for (int i=0; i < on->ngbr(iNgbr)->ngbr.m;i++){
	if (on->ngbr(iNgbr)->ngbr(i)->id == on-> ngbr(jNgbr)->id)
	  {
	    index_iNgbr_to_jNgbr=i; // index to j (in the neighbour list of i)
	    break;
	  }
      }
      // int index_jNgbr_to_iSelf=0;
      // for (int i=0; i < on->ngbr(jNgbr)->ngbr.m;i++){
      //   if (on->ngbr(jNgbr)->ngbr(i)->id == on-> ngbr(iSelf)->id)
      //   {
      //     index_jNgbr_to_iSelf=i; // index to current node (in the neighbour list of j)
      //     break;
      //   }
      // }
      int index_jNgbr_to_iNgbr=on->ngbr(jNgbr)->ngbr.m;
      for (int i=0; i < on->ngbr(jNgbr)->ngbr.m;i++){
	if (on->ngbr(jNgbr)->ngbr(i)->id == on-> ngbr(iNgbr)->id)
	  {
	    index_jNgbr_to_iNgbr=i; // index to i (in the neighbour list of j)
	    break;
	  }
      }
      // int index_jNgbr_to_jNgbr=on->ngbr(jNgbr)->ngbr.m;
      // for (int i=0; i < on->ngbr(jNgbr)->ngbr.m;i++){
      //   if (on->ngbr(jNgbr)->ngbr(i)->id == on-> ngbr(jNgbr)->id)
      //   {
      //     index_jNgbr_to_jNgbr=i; // index to j (in the neighbour list of j)
      //     break;
      //   }
      // }

      /////////////////////////////////////////////////////////////////
      // Determine whether iNgbr and jNgbr are well separated or not //
      /////////////////////////////////////////////////////////////////
      bool WellSeparated=false;
      if (index_iNgbr_to_jNgbr==on->ngbr(iNgbr)->ngbr.m && index_jNgbr_to_iNgbr==on->ngbr(jNgbr)->ngbr.m){
	WellSeparated=true; // iNgbr and jNgbr are neighbours of the current node, but not of each other (so they are in each others interaction list)

	for (int i=0; i < on->ngbr(iNgbr)->Ilist.m;i++){
	  if (on->ngbr(iNgbr)->Ilist(i)->id == on-> ngbr(jNgbr)->id)
	    {
	      index_iNgbr_to_jNgbr=i; // index to j (in the interaction list of i)
	      break;
	    }
	}
	for (int i=0; i < on->ngbr(jNgbr)->Ilist.m;i++){
	  if (on->ngbr(jNgbr)->Ilist(i)->id == on-> ngbr(iNgbr)->id)
	    {
	      index_jNgbr_to_iNgbr=i; // index to i (in the interaction list of j)
	      break;
	    }
	}
      }

      //////////////////////
      // iNgbr eliminated //
      //////////////////////
      if (on->ngbr(iNgbr)->IsEliminated){
	//////////////////////
	// jNgbr eliminated //
	//////////////////////
	if (on->ngbr(jNgbr)->IsEliminated){
	  // cout << "ij: 1 1"<< endl;

	  if (iNgbr==jNgbr){
	    if (iNgbr==iSelf){

	      //20200108	      start_operations=clock();

	      // M2L_ii
	      // on->ngbr(iNgbr)->M2L_operator[on->Ilist.m+index_iNgbr_to_iNgbr] += - alpha_self.partialPivLu().solve(MatrixXd::Identity(on->rank,on->rank));
	      // on->ngbr(iNgbr)->M2L_operator[on->Ilist.m+index_iNgbr_to_iNgbr] += - alpha_self_LU.solve(MatrixXd::Identity(on->rank,on->rank));
	      on->ngbr(iNgbr)->M2L_operator[on->Ilist.m+index_iNgbr_to_iNgbr] += - alpha_self_inv;

	      //20200108	      stop_operations=clock();
	      //20200108	      Time_operations+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);

	      //20200108	      number_operations_a++;
	      //20200108	      Time_operations_a+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);

                
	    }
	    else{

	      //20200108	      start_operations=clock();
	      // M2L_ii
	      // on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+index_iNgbr_to_iNgbr] += delta_ii - beta_is*alpha_self.partialPivLu().solve(gamma_si);
	      // on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+index_iNgbr_to_iNgbr] += delta_ii - beta_is*alpha_self_LU.solve(gamma_si);
	      // on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+index_iNgbr_to_iNgbr] += delta_ii - beta_is*xi_si;
	      // on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+index_iNgbr_to_iNgbr] += delta_ii[iNgbr] - beta_is[iNgbr]*xi_si[iNgbr];
	      on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+index_iNgbr_to_iNgbr].noalias() += -on->P2L_operator[iNgbr]*psi_si[iNgbr];

	      //20200108	      stop_operations=clock();
	      //20200108	      Time_operations+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);

	      //20200108	      number_operations_b++;
	      //20200108	      Time_operations_b+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);

	    }
	  }
	  else{
	    if(jNgbr==iSelf){

	      //20200108	      start_operations=clock();
	      // M2L_ji
	      // on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+index_iNgbr_to_jNgbr] += alpha_self.partialPivLu().solve(gamma_si);
	      // on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+index_iNgbr_to_jNgbr] += alpha_self_LU.solve(gamma_si);
	      // on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+index_iNgbr_to_jNgbr] += xi_si;
	      on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+index_iNgbr_to_jNgbr] += xi_si[iNgbr];

	      // M2L_ij
	      // on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+index_jNgbr_to_iNgbr] += beta_is*alpha_self.partialPivLu().solve(MatrixXd::Identity(on->rank,on->rank));
	      // on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+index_jNgbr_to_iNgbr] += beta_is*alpha_self_LU.solve(MatrixXd::Identity(on->rank,on->rank));
	      // on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+index_jNgbr_to_iNgbr] += beta_is*alpha_self_inv;
	      on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+index_jNgbr_to_iNgbr].noalias() += beta_is[iNgbr]*alpha_self_inv;

	      //20200108	      stop_operations=clock();
	      //20200108	      Time_operations+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);

	      //20200108	      number_operations_c++;
	      //20200108	      Time_operations_c+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);

#if defined(IFMM_PARALLELIZE)
	    } else if (iNgbr == iSelf) {
	      
	      // M2L_ij
	      on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+index_jNgbr_to_iNgbr] += xi_si[jNgbr];

	      // M2L_ji
	      on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+index_iNgbr_to_jNgbr].noalias() += beta_is[jNgbr]*alpha_self_inv;

#endif // IFMM_PARALLELIZE

	    }
	    else{
	      // MatrixXd beta_js = - (on->P2L_operator[jNgbr]*on->P2P_operator[iSelf].partialPivLu().solve(on->L2P_operator));
	      // MatrixXd gamma_sj = - on->P2M_operator*on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(jNgbr)->M2P_operator[index_jNgbr_to_iSelf]);
                
	      // // MatrixXd tau_sj=P2P_operator_self_LU.solve(on->ngbr(jNgbr)->M2P_operator[index_jNgbr_to_iSelf]);

	      // MatrixXd beta_js = - (on->P2L_operator[jNgbr]*P2P_operator_self_LU.solve(on->L2P_operator));
	      // //  MatrixXd beta_js = - (on->P2L_operator[jNgbr]*rho_self);
	      // MatrixXd gamma_sj = - on->P2M_operator*P2P_operator_self_LU.solve(on->ngbr(jNgbr)->M2P_operator[index_jNgbr_to_iSelf]);
	      // // MatrixXd gamma_sj = - on->P2M_operator*tau_sj;


	      if (WellSeparated==false){

		//20200108		start_operations=clock();
		// M2L_ji
		// on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+index_iNgbr_to_jNgbr] += -on->P2L_operator[jNgbr]*on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_iSelf])-beta_js*alpha_self.partialPivLu().solve(gamma_si);
		// on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+index_iNgbr_to_jNgbr] += -on->P2L_operator[jNgbr]*P2P_operator_self_LU.solve(on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_iSelf])-beta_js*alpha_self_LU.solve(gamma_si);
		// on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+index_iNgbr_to_jNgbr] += -on->P2L_operator[jNgbr]*tau_si-beta_js*alpha_self_LU.solve(gamma_si);
		// on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+index_iNgbr_to_jNgbr] += -on->P2L_operator[jNgbr]*tau_si-beta_js*xi_si;
		// on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+index_iNgbr_to_jNgbr] += -on->P2L_operator[jNgbr]*tau_si[iNgbr]-beta_is[jNgbr]*xi_si[iNgbr];
		on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+index_iNgbr_to_jNgbr].noalias() += -on->P2L_operator[jNgbr]*psi_si[iNgbr];



		// M2L_ij
		// on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+index_jNgbr_to_iNgbr] += -on->P2L_operator[iNgbr]*on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(jNgbr)->M2P_operator[index_jNgbr_to_iSelf])-beta_is*alpha_self.partialPivLu().solve(gamma_sj);
		// on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+index_jNgbr_to_iNgbr] += -on->P2L_operator[iNgbr]*P2P_operator_self_LU.solve(on->ngbr(jNgbr)->M2P_operator[index_jNgbr_to_iSelf])-beta_is*alpha_self_LU.solve(gamma_sj);
		// on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+index_jNgbr_to_iNgbr] += -on->P2L_operator[iNgbr]*tau_sj-beta_is*alpha_self_LU.solve(gamma_sj);
		// on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+index_jNgbr_to_iNgbr] += -on->P2L_operator[iNgbr]*tau_si[jNgbr]-beta_is[iNgbr]*xi_si[jNgbr];
		on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+index_jNgbr_to_iNgbr].noalias() += -on->P2L_operator[iNgbr]*psi_si[jNgbr];

		//20200108		stop_operations=clock();
		//20200108		Time_operations+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);

		//20200108		number_operations_d++;
		//20200108		Time_operations_d+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);


	      }
	      else{

		//20200108		start_operations=clock();
		// M2L_ji
		// on->ngbr(iNgbr)->M2L_operator[index_iNgbr_to_jNgbr] += -on->P2L_operator[jNgbr]*on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_iSelf])-beta_js*alpha_self.partialPivLu().solve(gamma_si);
		// on->ngbr(iNgbr)->M2L_operator[index_iNgbr_to_jNgbr] += -on->P2L_operator[jNgbr]*P2P_operator_self_LU.solve(on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_iSelf])-beta_js*alpha_self_LU.solve(gamma_si);
		// on->ngbr(iNgbr)->M2L_operator[index_iNgbr_to_jNgbr] += -on->P2L_operator[jNgbr]*tau_si-beta_js*alpha_self_LU.solve(gamma_si);
		// on->ngbr(iNgbr)->M2L_operator[index_iNgbr_to_jNgbr] += -on->P2L_operator[jNgbr]*tau_si-beta_js*xi_si;
		// on->ngbr(iNgbr)->M2L_operator[index_iNgbr_to_jNgbr] += -on->P2L_operator[jNgbr]*tau_si[iNgbr]-beta_is[jNgbr]*xi_si[iNgbr];
		on->ngbr(iNgbr)->M2L_operator[index_iNgbr_to_jNgbr].noalias() += -on->P2L_operator[jNgbr]*psi_si[iNgbr];



		// M2L_ij
		// on->ngbr(jNgbr)->M2L_operator[index_jNgbr_to_iNgbr] += -on->P2L_operator[iNgbr]*on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(jNgbr)->M2P_operator[index_jNgbr_to_iSelf])-beta_is*alpha_self.partialPivLu().solve(gamma_sj);
		// on->ngbr(jNgbr)->M2L_operator[index_jNgbr_to_iNgbr] += -on->P2L_operator[iNgbr]*P2P_operator_self_LU.solve(on->ngbr(jNgbr)->M2P_operator[index_jNgbr_to_iSelf])-beta_is*alpha_self_LU.solve(gamma_sj);
		// on->ngbr(jNgbr)->M2L_operator[index_jNgbr_to_iNgbr] += -on->P2L_operator[iNgbr]*tau_sj-beta_is*alpha_self_LU.solve(gamma_sj);
		// on->ngbr(jNgbr)->M2L_operator[index_jNgbr_to_iNgbr] += -on->P2L_operator[iNgbr]*tau_si[jNgbr]-beta_is[iNgbr]*xi_si[jNgbr];
		on->ngbr(jNgbr)->M2L_operator[index_jNgbr_to_iNgbr].noalias() += -on->P2L_operator[iNgbr]*psi_si[jNgbr];

		//20200108		stop_operations=clock();
		//20200108		Time_operations+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);

		//20200108		number_operations_e++;
		//20200108		Time_operations_e+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);

	      }
	    }
	  }
	}
	//////////////////////////////
	// jNgbr not yet eliminated //
	//////////////////////////////
	else{

	  didIdoJ(on, iSelf, iNgbr, jNgbr, index_iNgbr_to_jNgbr, index_jNgbr_to_iNgbr, rho_self, epsilon_rel_fillin, epsilon_rel_basis, sigma0_Aextended, machine_eps, rank_RSVD_min, P, LR_mode, Delay_update, outfile, curr_level, WellSeparated, alpha_self, alpha_self_inv, kappa_s, iota_s, beta_is, xi_si, psi_si, mu_is, zeta_si, eta_si, p2lto, p2pto);

	} // jNgbr not yet eliminated
      }
      //////////////////////////////
      // iNgbr not yet eliminated //
      //////////////////////////////
      else{
	//////////////////////
	// jNgbr eliminated //
	//////////////////////
	if (on->ngbr(jNgbr)->IsEliminated){

	  //20200108	  // cout << "ij: 0 1"<< endl;
	  //20200108	  cout << "Complete code..." << endl;
	  //20200108
	  //20200108	  // P2L_operator;
	  //20200108	  // M2P_operator;
          
#if defined(IFMM_PARALLELIZE)
	  didIdoJ(on, iSelf, jNgbr, iNgbr, index_jNgbr_to_iNgbr, index_iNgbr_to_jNgbr, rho_self, epsilon_rel_fillin, epsilon_rel_basis, sigma0_Aextended, machine_eps, rank_RSVD_min, P, LR_mode, Delay_update, outfile, curr_level, WellSeparated, alpha_self, alpha_self_inv, kappa_s, iota_s, beta_is, xi_si, psi_si, mu_is, zeta_si, eta_si, p2lto, p2pto);
#endif


	}
	//////////////////////////////
	// jNgbr not yet eliminated //
	//////////////////////////////
	else{

	  // cout << "ij: 0 0"<< endl;

	  if (iNgbr==jNgbr){
              

	    //20200108	    start_operations=clock();

	    // on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iNgbr]+=chi_ii - mu_is*(alpha_self.partialPivLu().solve(nu_si));
	    // on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iNgbr]+=chi_ii - mu_is*(alpha_self_LU.solve(nu_si));
	    // on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iNgbr]+=chi_ii - mu_is*zeta_si;
	    // on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iNgbr]+=chi_ii[iNgbr] - mu_is[iNgbr]*zeta_si[iNgbr];
	    on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iNgbr].noalias() += -on->P2P_operator[iNgbr]*eta_si[iNgbr];;

	    //20200108	    stop_operations=clock();
	    //20200108	    Time_operations+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);

	    //20200108	    number_operations_h++;
	    //20200108	    Time_operations_h+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);

	  }
	  else{

	    // // MatrixXd omicron_sj = P2P_operator_self_LU.solve(on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iSelf]);

	    // MatrixXd nu_sj = - on->P2M_operator*(on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iSelf]));
	    // MatrixXd mu_js = - on->P2P_operator[jNgbr]*(on->P2P_operator[iSelf].partialPivLu().solve(on->L2P_operator));
              
	    // MatrixXd nu_sj = - on->P2M_operator*(P2P_operator_self_LU.solve(on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iSelf]));
	    // // MatrixXd nu_sj = - on->P2M_operator*omicron_sj;
	    // MatrixXd mu_js = - on->P2P_operator[jNgbr]*(P2P_operator_self_LU.solve(on->L2P_operator));
	    // // MatrixXd mu_js = - on->P2P_operator[jNgbr]*rho_self;


	    // MatrixXd P2P_fillin_ji= - on->P2P_operator[jNgbr]*(on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf]));
	    // P2P_fillin_ji += - mu_js*(alpha_self.partialPivLu().solve(nu_si));
	    // MatrixXd P2P_fillin_ji= - on->P2P_operator[jNgbr]*(P2P_operator_self_LU.solve(on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf]));
	    // MatrixXd P2P_fillin_ji= - on->P2P_operator[jNgbr]*(omicron_si);
	    // MatrixXd P2P_fillin_ji= - on->P2P_operator[jNgbr]*(omicron_si[iNgbr]);
	    // P2P_fillin_ji += - mu_js*(alpha_self_LU.solve(nu_si));
	    // P2P_fillin_ji += - mu_js*zeta_si;
	    // P2P_fillin_ji += - mu_is[jNgbr]*zeta_si[iNgbr];

	    // MatrixXd P2P_fillin_ji= - on->P2P_operator[jNgbr]*omicron_si[iNgbr] - mu_is[jNgbr]*zeta_si[iNgbr];
	    //TEST MatrixXd P2P_fillin_ji= - on->P2P_operator[jNgbr]*eta_si[iNgbr];


	    // MatrixXd P2P_fillin_ij= - on->P2P_operator[iNgbr]*(on->P2P_operator[iSelf].partialPivLu().solve(on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iSelf]));
	    // P2P_fillin_ij += - mu_is*(alpha_self.partialPivLu().solve(nu_sj));
	    // MatrixXd P2P_fillin_ij= - on->P2P_operator[iNgbr]*(P2P_operator_self_LU.solve(on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iSelf]));
	    // MatrixXd P2P_fillin_ij= - on->P2P_operator[iNgbr]*omicron_sj;
	    // MatrixXd P2P_fillin_ij= - on->P2P_operator[iNgbr]*omicron_si[jNgbr];
	    // P2P_fillin_ij += - mu_is*(alpha_self_LU.solve(nu_sj));
	    // P2P_fillin_ij += - mu_is[iNgbr]*zeta_si[jNgbr];

	    // MatrixXd P2P_fillin_ij= - on->P2P_operator[iNgbr]*omicron_si[jNgbr] - mu_is[iNgbr]*zeta_si[jNgbr];
	    //TEST MatrixXd P2P_fillin_ij= - on->P2P_operator[iNgbr]*eta_si[jNgbr];




	    if (WellSeparated==false){

	      //20200108	      start_operations=clock();

	      // on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_jNgbr].noalias()+=P2P_fillin_ji;
	      // on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iNgbr].noalias()+=P2P_fillin_ij;

	      // // on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_jNgbr].noalias()+=- on->P2P_operator[jNgbr]*eta_si[iNgbr];
	      // // on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iNgbr].noalias()+=- on->P2P_operator[iNgbr]*eta_si[jNgbr];

	      //20200108	      MatrixXd temp1 = -on->P2P_operator[jNgbr]*eta_si[iNgbr];
	      //20200108
	      //20200108
	      //20200108	      on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_jNgbr].noalias() += temp1;
	      //20200108	      if (curr_level==l){
	      //20200108		on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iNgbr] += temp1.transpose(); // USE SYMMETRY
	      //20200108	      }
	      //20200108	      else{
	      //20200108		on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iNgbr].noalias() += - on->P2P_operator[iNgbr]*eta_si[jNgbr];
	      //20200108	      }

	      MatrixXd temp_ji = on->P2P_operator[jNgbr]*eta_si[iNgbr];
	      on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_jNgbr] += - temp_ji;
	      MatrixXd temp_ij = on->P2P_operator[iNgbr]*eta_si[jNgbr];
	      on->ngbr(jNgbr)->P2P_operator[index_jNgbr_to_iNgbr] += - temp_ij;

	      //20200108	      stop_operations=clock();
	      //20200108	      Time_operations+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);
                
	      //20200108	      number_operations_i_1++;
	      //20200108	      Time_operations_i_1+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);



	    }
	    else{
	      // timeType t0, t00, t1, t2, t3,t4,t5,t6,t7,t8,t9,t10,t11,t12; // Time variables

	      /////////////////////////////////////////////////////////////////
	      // iNGBR and jNGBR ARE WELL SEPARATED -> COMPRESS AND REDIRECT //
	      /////////////////////////////////////////////////////////////////

	      // cout << "Compress and redirect" << endl;

	      //20200108	      start_operations=clock();

	      // // MatrixXd P2P_fillin_ji= - on->P2P_operator[jNgbr]*eta_si[iNgbr];
	      // // MatrixXd P2P_fillin_ij= - on->P2P_operator[iNgbr]*eta_si[jNgbr];


	      MatrixXd P2P_fillin_ji;
	      MatrixXd P2P_fillin_ij;

	      if (LR_mode!=3){
		P2P_fillin_ji = - on->P2P_operator[jNgbr]*eta_si[iNgbr];
		P2P_fillin_ij = - on->P2P_operator[iNgbr]*eta_si[jNgbr];
	      }

	      //20200108	      stop_operations=clock();
	      //20200108	      Time_operations+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);
                
	      //20200108	      number_operations_i_2++;
	      //20200108	      Time_operations_i_2+= double(stop_operations-start_operations)/double(CLOCKS_PER_SEC);


	      int iSV=0;
	      int SV_U_comb_i=0;
	      int SV_V_comb_j=0;
	      int jSV=0;
	      int SV_U_comb_j=0;
	      int SV_V_comb_i=0;


	      MatrixXd U_i_prime, V_j_prime;
	      MatrixXd U_comb_i, V_comb_j;

	      MatrixXd U_j_prime, V_i_prime;
	      MatrixXd K_ij_prime, K_ji_prime;

	      MatrixXd U_comb_j, V_comb_i;

	      MatrixXd U_recomp_i, R_recomp_U_i;
	      MatrixXd U_recomp_j, R_recomp_U_j;
	      MatrixXd V_recomp_i, R_recomp_V_i;
	      MatrixXd V_recomp_j, R_recomp_V_j;

	      bool NeedUpdate=true;

	      MatrixXd Sigma_diag_U_comb_i, Sigma_diag_V_comb_i;    
	      MatrixXd Sigma_diag_U_comb_j, Sigma_diag_V_comb_j;



	      if (LR_mode==0){//SVD
#include "LR0-2.cpp"
	      }
	      if(LR_mode==1){// ACA
#include "LR1-2.cpp"
	      }
	      if (LR_mode==2){//RSVD
#include "LR2-2.cpp"
	      }
	      if (LR_mode==3){//RSVD - optimized
#include "LR3-2.cpp"
	      }


	      ///////////////////////////  
	      // UPDATE EDGES OF iNgbr //
	      ///////////////////////////
	      //20200108	      start_updates =clock();


	      if (NeedUpdate==true){

		MatrixXd Sigma_L2P_i_inverse=MatrixXd::Zero(on->ngbr(iNgbr)->rank,on->ngbr(iNgbr)->rank);
		MatrixXd Sigma_P2M_i_inverse=MatrixXd::Zero(on->ngbr(iNgbr)->rank,on->ngbr(iNgbr)->rank);
		for (int iRank=0; iRank< on->ngbr(iNgbr)->rank; iRank++){
		  Sigma_L2P_i_inverse(iRank,iRank)=1.0/on->ngbr(iNgbr)->Sigma_L2P(iRank,iRank);
		  Sigma_P2M_i_inverse(iRank,iRank)=1.0/on->ngbr(iNgbr)->Sigma_P2M(iRank,iRank);
		}
                
                
		MatrixXd Sigma_L2P_j_inverse=MatrixXd::Zero(on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->rank);
		MatrixXd Sigma_P2M_j_inverse=MatrixXd::Zero(on->ngbr(jNgbr)->rank,on->ngbr(jNgbr)->rank);
		for (int jRank=0; jRank< on->ngbr(jNgbr)->rank; jRank++){
		  Sigma_L2P_j_inverse(jRank,jRank)=1.0/on->ngbr(jNgbr)->Sigma_L2P(jRank,jRank);
		  Sigma_P2M_j_inverse(jRank,jRank)=1.0/on->ngbr(jNgbr)->Sigma_P2M(jRank,jRank);
		}

		MatrixXd Sigma_ij_inverse=MatrixXd::Zero(iSV,iSV);
		MatrixXd Sigma_ji_inverse=MatrixXd::Zero(jSV,jSV);
		for (int index_rank=0; index_rank<iSV; index_rank++){
		  Sigma_ij_inverse(index_rank,index_rank)=1.0/K_ij_prime(index_rank,index_rank);
		  Sigma_ji_inverse(index_rank,index_rank)=1.0/K_ji_prime(index_rank,index_rank);
		}

		// UPDATE P2M AND L2P
		MatrixXd R_V_i_T = Sigma_P2M_i_inverse*R_recomp_V_i.block(0,0,SV_V_comb_i,on->ngbr(iNgbr)->rank).transpose();                
		MatrixXd R_V_j_T = Sigma_P2M_j_inverse*R_recomp_V_j.block(0,0,SV_V_comb_j,on->ngbr(jNgbr)->rank).transpose();

		MatrixXd R_U_i = R_recomp_U_i.block(0,0,SV_U_comb_i,on->ngbr(iNgbr)->rank)*Sigma_L2P_i_inverse;
		MatrixXd R_U_j = R_recomp_U_j.block(0,0,SV_U_comb_j,on->ngbr(jNgbr)->rank)*Sigma_L2P_j_inverse;
                
		MatrixXd R_U_i_prime = R_recomp_U_i.block(0,on->ngbr(iNgbr)->rank,SV_U_comb_i,iSV)*Sigma_ij_inverse;
		MatrixXd R_U_j_prime = R_recomp_U_j.block(0,on->ngbr(jNgbr)->rank,SV_U_comb_j,jSV)*Sigma_ji_inverse;
                
		MatrixXd R_V_i_prime_T = Sigma_ji_inverse*R_recomp_V_i.block(0,on->ngbr(iNgbr)->rank,SV_V_comb_i,jSV).transpose();
		MatrixXd R_V_j_prime_T = Sigma_ij_inverse*R_recomp_V_j.block(0,on->ngbr(jNgbr)->rank,SV_V_comb_j,iSV).transpose();

		// UPDATE P2M AND L2P
		on->ngbr(iNgbr)->P2M_operator=V_recomp_i.transpose();
		on->ngbr(iNgbr)->L2P_operator=U_recomp_i;

		on->ngbr(iNgbr)->Sigma_P2M = Sigma_diag_V_comb_i;
		on->ngbr(iNgbr)->Sigma_L2P = Sigma_diag_U_comb_i;


#if defined(_OPENMP) && defined(IFMM_PARALLELIZE)
#pragma omp critical
		{
#endif

		  // NEIGHBOUR LIST
		  for (int lNgbr=0; lNgbr < on->ngbr(iNgbr)->ngbr.m; lNgbr++){
		    // M2L_li
		    // on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+lNgbr] = on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+lNgbr]*R_recomp_V_i.block(0,0,SV_V_comb_i,on->ngbr(iNgbr)->rank).transpose();
		    on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+lNgbr] = on->ngbr(iNgbr)->M2L_operator[on->ngbr(iNgbr)->Ilist.m+lNgbr]*R_V_i_T;
		    
		    
		    // also reverse
		    int index_lNgbr_to_iNgbr=0;
		    for (int i=0; i < on->ngbr(iNgbr)->ngbr(lNgbr)->ngbr.m;i++){
		      if (on->ngbr(iNgbr)->ngbr(lNgbr)->ngbr(i)->id == on->ngbr(iNgbr)->id)    
			{
			  index_lNgbr_to_iNgbr=i; // index l to i
			  break;
			}
		    }
		    // M2L_il
		    // on->ngbr(iNgbr)->ngbr(lNgbr)->M2L_operator[on->ngbr(iNgbr)->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_iNgbr] = R_recomp_U_i.block(0,0,SV_U_comb_i,on->ngbr(iNgbr)->rank)*on->ngbr(iNgbr)->ngbr(lNgbr)->M2L_operator[on->ngbr(iNgbr)->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_iNgbr];
		    on->ngbr(iNgbr)->ngbr(lNgbr)->M2L_operator[on->ngbr(iNgbr)->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_iNgbr] = R_U_i*on->ngbr(iNgbr)->ngbr(lNgbr)->M2L_operator[on->ngbr(iNgbr)->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_iNgbr];
		    
		  }


		  // INTERACTION LIST
		  for (int lInteract=0; lInteract < on->ngbr(iNgbr)->Ilist.m; lInteract++){
		    if(lInteract == index_iNgbr_to_jNgbr){
		      
		      // R_recomp_U_j.block(0,on->ngbr(jNgbr)->rank,SV_U_comb_j,jSV)*K_ji_prime* R_recomp_V_i.block(0,on->ngbr(iNgbr)->rank,SV_V_comb_i,jSV).transpose();
		      // M2L_ji
		      // on->ngbr(iNgbr)->M2L_operator[index_iNgbr_to_jNgbr] = R_recomp_U_j.block(0,0,SV_U_comb_j,on->ngbr(jNgbr)->rank)*on->ngbr(iNgbr)->M2L_operator[index_iNgbr_to_jNgbr]*R_recomp_V_i.block(0,0,SV_V_comb_i,on->ngbr(iNgbr)->rank).transpose()+R_recomp_U_j.block(0,on->ngbr(jNgbr)->rank,SV_U_comb_j,jSV)*K_ji_prime* R_recomp_V_i.block(0,on->ngbr(iNgbr)->rank,SV_V_comb_i,jSV).transpose();
		      on->ngbr(iNgbr)->M2L_operator[index_iNgbr_to_jNgbr] = R_U_j*on->ngbr(iNgbr)->M2L_operator[index_iNgbr_to_jNgbr]*R_V_i_T+R_U_j_prime*K_ji_prime*R_V_i_prime_T;
		      
		    }
		    else{
		      // M2L_li
		      // on->ngbr(iNgbr)->M2L_operator[lInteract] =  on->ngbr(iNgbr)->M2L_operator[lInteract]*R_recomp_V_i.block(0,0,SV_V_comb_i,on->ngbr(iNgbr)->rank).transpose();
		      on->ngbr(iNgbr)->M2L_operator[lInteract] =  on->ngbr(iNgbr)->M2L_operator[lInteract]*R_V_i_T;
		      
		      
		      // also reverse
		      int index_lNgbr_to_iNgbr=0;
		      for (int i=0; i < on->ngbr(iNgbr)->Ilist(lInteract)->Ilist.m;i++){
			if (on->ngbr(iNgbr)->Ilist(lInteract)->Ilist(i)->id == on->ngbr(iNgbr)->id)    
			  {
			    index_lNgbr_to_iNgbr=i; // index l to i
			    break;
			  }
		      }

		      // M2L_il
		      // on->ngbr(iNgbr)->Ilist(lInteract)->M2L_operator[index_lNgbr_to_iNgbr] = R_recomp_U_i.block(0,0,SV_U_comb_i,on->ngbr(iNgbr)->rank)*on->ngbr(iNgbr)->Ilist(lInteract)->M2L_operator[index_lNgbr_to_iNgbr];
		      on->ngbr(iNgbr)->Ilist(lInteract)->M2L_operator[index_lNgbr_to_iNgbr] = R_U_i*on->ngbr(iNgbr)->Ilist(lInteract)->M2L_operator[index_lNgbr_to_iNgbr];
		      
		      
		    }
		    
		  }

#if defined(_OPENMP) && defined(IFMM_PARALLELIZE)
		} // omp critical
#endif	

		if(curr_level>2){
		  // UPDATE PARENT OPERATORS
		  int id_iNgbr_local=0;
		  for (int iChild=0;iChild<8;iChild++){
		    if (on->ngbr(iNgbr)->parent->child[iChild] != NULL){
		      if (on->ngbr(iNgbr)->parent->child[iChild]->id == on->ngbr(iNgbr)->id){
			id_iNgbr_local=iChild;
			break;
		      }
		    }
		  }

		  // on->ngbr(iNgbr)->parent->M2M_operator[id_iNgbr_local] = on->ngbr(iNgbr)->parent->M2M_operator[id_iNgbr_local]*R_recomp_V_i.block(0,0,SV_V_comb_i,on->ngbr(iNgbr)->rank).transpose();
		  // on->ngbr(iNgbr)->parent->L2L_operator[id_iNgbr_local] = R_recomp_U_i.block(0,0,SV_U_comb_i,on->ngbr(iNgbr)->rank)*on->ngbr(iNgbr)->parent->L2L_operator[id_iNgbr_local];

		  on->ngbr(iNgbr)->parent->M2M_operator[id_iNgbr_local] = on->ngbr(iNgbr)->parent->M2M_operator[id_iNgbr_local]*R_V_i_T;
		  on->ngbr(iNgbr)->parent->L2L_operator[id_iNgbr_local] = R_U_i*on->ngbr(iNgbr)->parent->L2L_operator[id_iNgbr_local];
		}


		///////////////////////////  
		// UPDATE EDGES OF jNgbr //
		///////////////////////////  

		// printf("SV_U_comb_j,SV_V_comb_j: %d %d \n",SV_U_comb_j,SV_V_comb_j);
                
		// UPDATE P2M AND L2P
		on->ngbr(jNgbr)->P2M_operator=V_recomp_j.transpose();
		on->ngbr(jNgbr)->L2P_operator=U_recomp_j;

		on->ngbr(jNgbr)->Sigma_P2M = Sigma_diag_V_comb_j;
		on->ngbr(jNgbr)->Sigma_L2P = Sigma_diag_U_comb_j;

		// cout << "on->ngbr(jNgbr)->P2M_operator.norm()" << endl << on->ngbr(jNgbr)->P2M_operator.norm() << endl;


#if defined(_OPENMP) && defined(IFMM_PARALLELIZE)
#pragma omp critical
		{
#endif
		  
		  // NEIGHBOUR LIST
		  for (int lNgbr=0; lNgbr < on->ngbr(jNgbr)->ngbr.m; lNgbr++){
		    
		    // M2L_lj
		    // on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+lNgbr] = on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+lNgbr]*R_recomp_V_j.block(0,0,SV_V_comb_j,on->ngbr(jNgbr)->rank).transpose();
		    on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+lNgbr] = on->ngbr(jNgbr)->M2L_operator[on->ngbr(jNgbr)->Ilist.m+lNgbr]*R_V_j_T;
		    
		    
		    // also reverse
		    int index_lNgbr_to_jNgbr=0;
		    for (int i=0; i < on->ngbr(jNgbr)->ngbr(lNgbr)->ngbr.m;i++){
		      if (on->ngbr(jNgbr)->ngbr(lNgbr)->ngbr(i)->id == on->ngbr(jNgbr)->id)    
			{
			  index_lNgbr_to_jNgbr=i; // index l to i
			  break;
			}
		    }
		    // M2L_jl
		    // on->ngbr(jNgbr)->ngbr(lNgbr)->M2L_operator[on->ngbr(jNgbr)->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_jNgbr] = R_recomp_U_j.block(0,0,SV_U_comb_j,on->ngbr(jNgbr)->rank)*on->ngbr(jNgbr)->ngbr(lNgbr)->M2L_operator[on->ngbr(jNgbr)->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_jNgbr];
		    on->ngbr(jNgbr)->ngbr(lNgbr)->M2L_operator[on->ngbr(jNgbr)->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_jNgbr] = R_U_j*on->ngbr(jNgbr)->ngbr(lNgbr)->M2L_operator[on->ngbr(jNgbr)->ngbr(lNgbr)->Ilist.m+index_lNgbr_to_jNgbr];
		    
		    
		  }
		  
		  
		  // INTERACTION LIST
		  for (int lInteract=0; lInteract < on->ngbr(jNgbr)->Ilist.m; lInteract++){
		    if(lInteract == index_jNgbr_to_iNgbr){
		      // M2L_ij
		    // on->ngbr(jNgbr)->M2L_operator[index_jNgbr_to_iNgbr] = R_recomp_U_i.block(0,0,SV_U_comb_i,on->ngbr(iNgbr)->rank)*on->ngbr(jNgbr)->M2L_operator[index_jNgbr_to_iNgbr]*R_recomp_V_j.block(0,0,SV_V_comb_j,on->ngbr(jNgbr)->rank).transpose()+R_recomp_U_i.block(0,on->ngbr(iNgbr)->rank,SV_U_comb_i,iSV)*K_ij_prime*R_recomp_V_j.block(0,on->ngbr(jNgbr)->rank,SV_V_comb_j,iSV).transpose(); 
		      on->ngbr(jNgbr)->M2L_operator[index_jNgbr_to_iNgbr] = R_U_i*on->ngbr(jNgbr)->M2L_operator[index_jNgbr_to_iNgbr]*R_V_j_T+R_U_i_prime*K_ij_prime*R_V_j_prime_T; 
		      
		    }
		    else{
		      // M2L_lj
		      // on->ngbr(jNgbr)->M2L_operator[lInteract] = on->ngbr(jNgbr)->M2L_operator[lInteract]*R_recomp_V_j.block(0,0,SV_V_comb_j,on->ngbr(jNgbr)->rank).transpose();
		      on->ngbr(jNgbr)->M2L_operator[lInteract] = on->ngbr(jNgbr)->M2L_operator[lInteract]*R_V_j_T;
		      
		      
		      // also reverse
		      int index_lNgbr_to_jNgbr=0;
		      for (int i=0; i < on->ngbr(jNgbr)->Ilist(lInteract)->Ilist.m;i++){
			if (on->ngbr(jNgbr)->Ilist(lInteract)->Ilist(i)->id == on->ngbr(jNgbr)->id)    
			  {
			    index_lNgbr_to_jNgbr=i; // index l to i
			    break;
			  }
		      }
		      // M2L_jl
		      // on->ngbr(jNgbr)->Ilist(lInteract)->M2L_operator[index_lNgbr_to_jNgbr] = R_recomp_U_j.block(0,0,SV_U_comb_j,on->ngbr(jNgbr)->rank)*on->ngbr(jNgbr)->Ilist(lInteract)->M2L_operator[index_lNgbr_to_jNgbr];
		      on->ngbr(jNgbr)->Ilist(lInteract)->M2L_operator[index_lNgbr_to_jNgbr] = R_U_j*on->ngbr(jNgbr)->Ilist(lInteract)->M2L_operator[index_lNgbr_to_jNgbr];
		      
		    }
		  }

#if defined(_OPENMP) && defined(IFMM_PARALLELIZE)
		} // omp critical
#endif	

		if (curr_level>2){
		  // UPDATE PARENT OPERATORS
		  int id_jNgbr_local=0;
		  for (int iChild=0;iChild<8;iChild++){
		    if (on->ngbr(jNgbr)->parent->child[iChild] != NULL){
		      if (on->ngbr(jNgbr)->parent->child[iChild]->id == on->ngbr(jNgbr)->id){
			id_jNgbr_local=iChild;
			break;
		      }
		    }
		  }
		  // on->ngbr(jNgbr)->parent->M2M_operator[id_jNgbr_local] = on->ngbr(jNgbr)->parent->M2M_operator[id_jNgbr_local]*R_recomp_V_j.block(0,0,SV_V_comb_j,on->ngbr(jNgbr)->rank).transpose();
		  // on->ngbr(jNgbr)->parent->L2L_operator[id_jNgbr_local] = R_recomp_U_j.block(0,0,SV_U_comb_j,on->ngbr(jNgbr)->rank)*on->ngbr(jNgbr)->parent->L2L_operator[id_jNgbr_local];

		  on->ngbr(jNgbr)->parent->M2M_operator[id_jNgbr_local] = on->ngbr(jNgbr)->parent->M2M_operator[id_jNgbr_local]*R_V_j_T;
		  on->ngbr(jNgbr)->parent->L2L_operator[id_jNgbr_local] = R_U_j*on->ngbr(jNgbr)->parent->L2L_operator[id_jNgbr_local];
                  
                  
		}

		// UPDATE RANK OF iNGBR AND jNGBR
		// printf("Before - on->ngbr(iNgbr)->rank: %d\n",on->ngbr(iNgbr)->rank);
		// printf("Before - on->ngbr(jNgbr)->rank: %d\n",on->ngbr(jNgbr)->rank);

		on->ngbr(iNgbr)->rank=SV_U_comb_i;
		on->ngbr(jNgbr)->rank=SV_U_comb_j;

		// printf("After - on->ngbr(iNgbr)->rank: %d\n",on->ngbr(iNgbr)->rank);
		// printf("After - on->ngbr(jNgbr)->rank: %d\n",on->ngbr(jNgbr)->rank);

	      }

	      //20200108	      stop_updates =clock();
	      //20200108	      Time_updates += double(stop_updates-start_updates)/double(CLOCKS_PER_SEC);


	    }
	  } 
	}
      }

      //20200108      stop_ij=clock();
      //20200108      Time_ij+= double(stop_ij-start_ij)/double(CLOCKS_PER_SEC);

    } // jNgbr



    ////////////////
    // UPDATE RHS //
    ////////////////
    //20200108    start_RHS =clock();
    if (on->ngbr(iNgbr)->IsEliminated){
      if (iNgbr==iSelf){
	// DO NOTHING YET, AS THIS WOULD AFFECT THE UPDATES OF THE OTHER NODES
      }
      else{
	// on->ngbr(iNgbr)->RHS += -on->P2L_operator[iNgbr]*(on->P2P_operator[iSelf].partialPivLu().solve(on->RHS_leaf))+beta_is*alpha_self.partialPivLu().solve(on->P2M_operator*(on->P2P_operator[iSelf].partialPivLu().solve(on->RHS_leaf)));
	// on->ngbr(iNgbr)->RHS += -on->P2L_operator[iNgbr]*(P2P_operator_self_LU.solve(on->RHS_leaf))+beta_is*alpha_self_LU.solve(on->P2M_operator*(P2P_operator_self_LU.solve(on->RHS_leaf)));
	// on->ngbr(iNgbr)->RHS += -on->P2L_operator[iNgbr]*(kappa_s)+beta_is*alpha_self_LU.solve(on->P2M_operator*(kappa_s));
	// on->ngbr(iNgbr)->RHS += -on->P2L_operator[iNgbr]*(kappa_s)+beta_is[iNgbr]*alpha_self_LU.solve(on->P2M_operator*(kappa_s));
	on->ngbr(iNgbr)->RHS.noalias() += -on->P2L_operator[iNgbr]*kappa_s+beta_is[iNgbr]*iota_s;
      }
    }
    else{
      // on->ngbr(iNgbr)->RHS_leaf += - on->P2P_operator[iNgbr]*on->P2P_operator[iSelf].partialPivLu().solve(on->RHS_leaf)+mu_is*(alpha_self.partialPivLu().solve(on->P2M_operator*(on->P2P_operator[iSelf].partialPivLu().solve(on->RHS_leaf))));
      // on->ngbr(iNgbr)->RHS_leaf += - on->P2P_operator[iNgbr]*P2P_operator_self_LU.solve(on->RHS_leaf)+mu_is*(alpha_self_LU.solve(on->P2M_operator*(P2P_operator_self_LU.solve(on->RHS_leaf))));
      // on->ngbr(iNgbr)->RHS_leaf += - on->P2P_operator[iNgbr]*kappa_s+mu_is*(alpha_self_LU.solve(on->P2M_operator*(kappa_s)));
      // on->ngbr(iNgbr)->RHS_leaf += - on->P2P_operator[iNgbr]*kappa_s+mu_is[iNgbr]*(alpha_self_LU.solve(on->P2M_operator*(kappa_s)));
      on->ngbr(iNgbr)->RHS_leaf.noalias() += - on->P2P_operator[iNgbr]*kappa_s+mu_is[iNgbr]*iota_s;
    }
  
  }

  // PERFORM UPDATE OF iSELF NOW
  // if (IsLeaf){
  if (on->RHS.rows() != on->rank){
    on->RHS = VectorXd::Zero(on->rank);
  }
  // on->RHS += - alpha_self.partialPivLu().solve(on->P2M_operator*(on->P2P_operator[iSelf].partialPivLu().solve(on->RHS_leaf)));
  // on->RHS += - alpha_self_LU.solve(on->P2M_operator*(P2P_operator_self_LU.solve(on->RHS_leaf)));
  // on->RHS += - alpha_self_LU.solve(on->P2M_operator*(kappa_s));
  on->RHS += - iota_s;

  // printf("on->rank:  %d\n",on->rank);

  // CURRENT NODE IS ELIMINATED FROM THE GRAPH
  // t1=Timer();

  // cout << "Time for eliminating node:" << endl << t1-t0 << endl;

  //20200108  stop_RHS =clock();
  //20200108  Time_RHS += double(stop_RHS-start_RHS)/double(CLOCKS_PER_SEC);

}

void IFMM_Matrix::eliminate_level(double & sigma0_Aextended, int curr_level, bool IsLeaf, ofstream& outfile){

  // int LR_mode=0; // SVD
  // int LR_mode=1; // ACA
  const double machine_eps = std::numeric_limits<double>::epsilon();
  // cout << "machine_eps:" << endl << machine_eps << endl;

  // (LR_mode==0) ? printf("SVD\n"):printf("ACA\n");
  if (LR_mode==0){printf("SVD\n");}
  if (LR_mode==1){printf("ACA\n");}
  if (LR_mode==2){printf("RSVD\n");}
  if (LR_mode==3){printf("RSVD - optimized\n");}

  const int rank_RSVD_min=15;

  MatrixXd P(rank_RSVD_min, rank_RSVD_min);
  sample_gaussian(P);

  //20200108  clock_t start_SVD, stop_SVD; double TimeSVD=0;
  //20200108  clock_t start_ACA, stop_ACA; double Time_ACA=0;
  //20200108  clock_t start_precomp, stop_precomp; double Time_precomp=0;
  //20200108  clock_t start_updates, stop_updates; double Time_updates=0;
  //20200108  clock_t start_RHS, stop_RHS; double Time_RHS=0;
  //20200108  clock_t start_ij, stop_ij; double Time_ij=0;
  //20200108  clock_t start_operations, stop_operations; double Time_operations=0;
  //20200108
  //20200108  double number_operations_a=0;  double number_operations_b=0; double number_operations_c=0;
  //20200108  double number_operations_d=0;  double number_operations_e=0; double number_operations_f=0;
  //20200108  double number_operations_g_1=0; double number_operations_g_2=0; double number_operations_h=0;
  //20200108  double number_operations_i_1=0; double number_operations_i_2=0;
  //20200108
  //20200108  double Time_operations_a=0; double Time_operations_b=0; double Time_operations_c=0;
  //20200108  double Time_operations_d=0; double Time_operations_e=0; double Time_operations_f=0;
  //20200108  double Time_operations_g_1=0; double Time_operations_g_2=0; double Time_operations_h=0;
  //20200108  double Time_operations_i_1=0; double Time_operations_i_2=0;

#if defined(IFMM_PARALLELIZE)

  seek2_concurrent_nodes(curr_level, IFMM_PARALLELIZE_CONCURRENCY_SEP);
  
  // Loop over all nodes of the current level
  for (int k = lvl_idx(curr_level); k < lvl_idx(curr_level + 1); ++ k) {
    
    oct_node *on = lvl_p(k);
    
#if defined(_OPENMP)
#pragma omp parallel for //num_threads(IFMM_PARALLELIZE)
#endif
    for (int i = 0; i < on->concurrent.m; i ++) {
      eliminate_level_node(sigma0_Aextended, curr_level, outfile, machine_eps, rank_RSVD_min, P, on->concurrent(i)->id);
    }

  }

#else // ! IFMM_PARALLELIZE
  
  // Loop over all nodes of the current level
  for (int k = lvl_idx(curr_level); k < lvl_idx(curr_level + 1); ++k) {
    eliminate_level_node(sigma0_Aextended, curr_level, outfile, machine_eps, rank_RSVD_min, P, k);
  }
  
#endif // ! IFMM_PARALLELIZE

  //20200108  outfile << "Time_RSVD: " << TimeSVD << endl;
  //20200108  outfile << "Time_ACA: " << Time_ACA << endl;
  //20200108  outfile << "Time_precomp: " << Time_precomp << endl;
  //20200108  outfile << "Time_updates: " << Time_updates << endl;
  //20200108  outfile << "Time_RHS: " << Time_RHS << endl;
  //20200108  outfile << "Time_ij: " << Time_ij << endl;
  //20200108  outfile << "Time_operations: " << Time_operations << endl;
  //20200108
  //20200108  outfile <<"number_operations_a,Time_operations_a:" <<  number_operations_a <<  " , " << Time_operations_a << endl;
  //20200108  outfile <<"number_operations_b,Time_operations_b:" << number_operations_b <<  " , "<< Time_operations_b << endl;
  //20200108  outfile <<"number_operations_c,Time_operations_c:" << number_operations_c <<  " , "<< Time_operations_c << endl;
  //20200108  outfile <<"number_operations_d,Time_operations_d:" << number_operations_d <<  " , "<< Time_operations_d << endl;
  //20200108  outfile <<"number_operations_e,Time_operations_e:" << number_operations_e <<  " , "<< Time_operations_e<< endl;
  //20200108  outfile <<"number_operations_f,Time_operations_f:" << number_operations_f <<  " , "<< Time_operations_f<< endl;
  //20200108  outfile <<"number_operations_g_1,Time_operations_g_1:" << number_operations_g_1 <<  " , "<< Time_operations_g_1<< endl;
  //20200108  outfile <<"number_operations_g_2,Time_operations_g_2:" << number_operations_g_2 <<  " , "<< Time_operations_g_2<< endl;  
  //20200108  outfile <<"number_operations_h,Time_operations_h: " << number_operations_h <<  " , "<< Time_operations_h<< endl;
  //20200108  outfile <<"number_operations_i_1,Time_operations_i_1:" << number_operations_i_1 <<  " , "<< Time_operations_i_1<< endl;
  //20200108  outfile <<"number_operations_i_2,Time_operations_i_2:" << number_operations_i_2 <<  " , "<< Time_operations_i_2<< endl;


}

 
void IFMM_Matrix::solve_top_level(){

    // cout << "Solve top level ... " << endl;

    int curr_level = 2;
    int nNodes_before=lvl_idx(curr_level); // NUMBER OF  NODES AT LOWER LEVELS
    int nNodes=lvl_idx(curr_level + 1)-lvl_idx(curr_level); // NUMBER OF NODES AT THE CURRENT LEVEL

    // CUMULATIVE RANK
    ivec rank_cum(nNodes+1);
    rank_cum(0)=0;

   for (int k = lvl_idx(curr_level); k < lvl_idx(curr_level + 1); ++k) {
      oct_node * on = lvl_p(k);
      rank_cum(k-nNodes_before+1)=rank_cum(k-nNodes_before)+on->rank;
    }

    // ASSEMBLE DENSE MATRIX (CONSISTING OF M2L INTERACTIONS)
    MatrixXd DenseMatrix_top_level = MatrixXd::Zero(rank_cum(nNodes),rank_cum(nNodes));
    VectorXd RHS_top_level = VectorXd::Zero(rank_cum(nNodes));

#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (int k = lvl_idx(curr_level); k < lvl_idx(curr_level + 1); ++k) {

      // CURRENT NODE
      oct_node * on = lvl_p(k);

      if (reuse==false){
        // INTERACTION LIST
        for (int iInteract = 0; iInteract < on->Ilist.m; ++iInteract){
          DenseMatrix_top_level.block(rank_cum(on->Ilist(iInteract)->id-nNodes_before),rank_cum(on->id-nNodes_before),on->Ilist(iInteract)->rank,on->rank)=on->M2L_operator[iInteract];
        }      

        // NEIGHBOUR LIST
        for (int iNgbr=0; iNgbr < on-> ngbr.m; iNgbr++){
          DenseMatrix_top_level.block(rank_cum(on->ngbr(iNgbr)->id-nNodes_before),rank_cum(on->id-nNodes_before),on->ngbr(iNgbr)->rank,on->rank)=on->M2L_operator[on->Ilist.m+iNgbr];
        }
      }

      // ASSEMBLE RHS
      RHS_top_level.block(rank_cum(on->id-nNodes_before),0,on->rank,1) = on->RHS;

    }

    // SOLVE TOP SYSTEM OF EQUATIONS AT LEVEL 2 (DENSE MATRIX) 
    // VectorXd Moments_toplevel = DenseMatrix_top_level.partialPivLu().solve(RHS_top_level);
    if (reuse==false){
      DenseMatrix_top_level_LU = DenseMatrix_top_level.partialPivLu();
      // VectorXd Moments_toplevel = DenseMatrix_top_level.partialPivLu().solve(RHS_top_level);
    }

    VectorXd Moments_toplevel = DenseMatrix_top_level_LU.solve(RHS_top_level);
    // cout << "Here is the Matrix DenseMatrix_top_level:\n" << DenseMatrix_top_level << endl;
    // cout << "Here is the Vector Moments_toplevel:\n" << Moments_toplevel << endl;


    // TRANSFER TO ON->Y
#if defined(_OPENMP)
#pragma omp parallel for
#endif
   for (int k = lvl_idx(curr_level); k < lvl_idx(curr_level + 1); ++k) {
      oct_node * on = lvl_p(k);
      for (int i=0; i<on->rank; i++){
        on->y(i) = Moments_toplevel(rank_cum(on->id-nNodes_before)+i);
      }
    }

    if (reuse==false){
      reuse=true;
    }
    
}


void IFMM_Matrix::downward_pass_node(dvec &x_iFMM, ofstream &outfile, const int curr_level, const int k)
{
  // CURRENT NODE
  oct_node * on = lvl_p(k);
  // printf("on->id, on->rank: %d %d \n",on->id, on->rank);
  // printf("on->id, on->nDOF: %d %d \n",on->id, on->nDOF);
  
  outfile << "Node: " << on->id << ", Rank (final): " << on->rank << endl;

  int iSelf=0; // index to current node (in the neighbour list of the current node)
  for (int iNgbr=0; iNgbr < on-> ngbr.m; iNgbr++){
    if (on->ngbr(iNgbr)->id==on->id){
      iSelf=iNgbr;
      break;
    }
  }


  // Eigen::PartialPivLU<MatrixXd> P2P_operator_self_LU = on->P2P_operator[iSelf].partialPivLu();

  // MatrixXd alpha_self=MatrixXd::Zero(on->rank,on->rank);
  // alpha_self = - on->P2M_operator*(on->P2P_operator[iSelf].partialPivLu().solve(on->L2P_operator));
  // alpha_self = - on->P2M_operator*(P2P_operator_self_LU.solve(on->L2P_operator));
  VectorXd temp=VectorXd::Zero(on->nDOF);


  // VectorXd temptest1=VectorXd::Zero(on->nDOF);
  // VectorXd temptest2=VectorXd::Zero(on->nDOF);


  // NEIGHBOUR LIST
  for (int iNgbr = 0; iNgbr < on->ngbr.m; ++iNgbr) {
    int index_iNgbr_to_iSelf=0;
    for (int i=0; i < on->ngbr(iNgbr)->ngbr.m;i++){
      if (on->ngbr(iNgbr)->ngbr(i)->id == on->id)
	{
	  index_iNgbr_to_iSelf=i; // index to current node (in the neighbour list of i)
	  break;
	}
    }

    if (iNgbr!=iSelf){// not self!
      if(on->ngbr(iNgbr)->IsEliminated){
	temp.noalias()+=on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_iSelf]*on->ngbr(iNgbr)->y;
	// temptest1+=on->ngbr(iNgbr)->M2P_operator[index_iNgbr_to_iSelf]*on->ngbr(iNgbr)->y;
      }
      else{
	temp.noalias()+=on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf]*on->ngbr(iNgbr)->x;   
	// temptest2+=on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf]*on->ngbr(iNgbr)->x; 
	// cout <<"on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf]" << on->ngbr(iNgbr)->P2P_operator[index_iNgbr_to_iSelf] << endl;           
      }
    }
  }

  // cout << "temptest1 " << endl << temptest1 << endl;
  // cout << "temptest2 " << endl << temptest2 << endl;
  // cout << "temp " << endl << temp << endl;

  // if (on->id==5){
  // cout << "P2P_operator_self_LU.solve(on->RHS_leaf-temp)" << P2P_operator_self_LU.solve(on->RHS_leaf-temp) << endl;
  // cout << "on->RHS_leaf-temp" << on->RHS_leaf-temp << endl;
  // }


  // cout << "alpha_self" << endl << alpha_self << endl;
  // cout << "alpha_self_inv" << endl << alpha_self.partialPivLu().solve(MatrixXd::Identity(on->rank,on->rank)) << endl;
  // cout << "on->P2M_operator" << endl << on->P2M_operator << endl;
      

  // CALCULATE LOCAL
  // on->z = alpha_self.partialPivLu().solve(on->y - on->P2M_operator*(on->P2P_operator[iSelf].partialPivLu().solve(on->RHS_leaf-temp)));
  on->z = on->alpha_self_LU.solve(on->y - on->P2M_operator*(on->P2P_operator_self_LU.solve(on->RHS_leaf-temp)));

  // CALCULATE CHARGES
  // on->x = on->P2P_operator[iSelf].partialPivLu().solve(on->RHS_leaf - temp - on->L2P_operator*on->z);
  on->x = on->P2P_operator_self_LU.solve(on->RHS_leaf - temp - on->L2P_operator*on->z);


  // CURRENT NODE IS BACK-SUBSTITUTED
  on->IsEliminated=0;

  // VectorXd test_ones=VectorXd::Ones(on->nDOF);
  // double test= (test_ones - on->x).norm() / test_ones.norm(); 
  // cout << "The relative error is:\n" << test << endl;

  // cout << "on->x " << on->x  << endl;
  // cout << "on->y " << on->y  << endl;
  // cout << "on->z " << on->z  << endl;

  if (curr_level!=l){
    for (int iChild=0; iChild<8; iChild++){
      if (on->child[iChild] != NULL){
	on->child[iChild]->y = on->x.block(on->rank_children_cum(iChild),0,on->child[iChild]->rank,1);
      }
    }
  }
  else{
    list<int>::iterator it = on->pidx.begin();
    for (int i = 0; i < on->nDOF; ++it, ++i){
      x_iFMM(*it) = on->x(i);  
    }
    // cout << "on->x " << on->x  << endl;
    // cout << "on->y " << on->y  << endl;
    // cout << "on->z " << on->z  << endl;

    // cout << "on->RHS_leaf - temp - on->L2P_operator*on->z" << on->RHS_leaf - temp - on->L2P_operator*on->z << endl;
    // cout << "on->P2P_operator[iSelf]" << endl << on->P2P_operator[iSelf] << endl;

    // VectorXd corrsol = VectorXd::Ones(on->nDOF);
    // VectorXd diff = corrsol - on->x;
        
    // double aaa=abs(on->x-corrsol)/abs(corrsol);
    // printf("relerror node: %f \n",diff.norm()/corrsol.norm());

    // if (diff.norm()/corrsol.norm() > 1e-1 || k==lvl_idx(curr_level + 1)-1){
    //   cout << "on->x " << on->x  << endl;
    //   cout << "on->y " << on->y  << endl;
    //   cout << "on->z " << on->z  << endl;
    //   cout << "on->RHS_leaf - temp - on->L2P_operator*on->z" << on->RHS_leaf - temp - on->L2P_operator*on->z << endl;

    //   cout << "on->P2P_operator[iSelf]" << endl << on->P2P_operator[iSelf] << endl;
    //   cout << "alpha_self" << endl << alpha_self << endl;
    //   cout << "on->RHS_leaf" << endl << on->RHS_leaf << endl;
    //   cout << "temp" << endl << temp << endl;

    //   cout << "on->P2M_operator" << endl << on->P2M_operator << endl;
    //   cout << "on->L2P_operator" << endl << on->L2P_operator << endl;




    //    for(int j=0;j<on->pnt.m;++j){
    //   //   printf("point(i).x, point(i).y, point(i).z : %f %f %f  ; \n",on->pnt.p[j].x,on->pnt.p[j].y,on->pnt.p[j].z);
    //     printf("on->sigma(j): %f\n",on->sigma(j));
    //   }

    //   list<int>::iterator it = on->pidx.begin();
    //    for (int i = 0; i < on->nDOF; ++it, ++i){
    //     printf("*it: %d\n",*it);
    //   }

          
    // }


  }

}


void IFMM_Matrix::downward_pass(dvec & x_iFMM,ofstream& outfile){
  
  ///////////////////
  // DOWNWARD PASS //
  ///////////////////
  
  outfile << endl;
  
  for (int curr_level=2; curr_level <= l;curr_level++){
    // printf("curr_level, l: %d %d\n",curr_level, l);
    
    // cout << "Downward_pass level " << curr_level << " ... " << endl;
    outfile << "Level " << curr_level << endl;
    
    for (int k=lvl_idx(curr_level + 1)-1; k >= lvl_idx(curr_level); --k) {

#if defined(IFMM_PARALLELIZE)

      oct_node *on = lvl_p(k);

#if defined(_OPENMP)
#pragma omp parallel for //num_threads(IFMM_PARALLELIZE)
#endif
      for (int i = 0; i < on->concurrent.m; i ++) {
	downward_pass_node(x_iFMM, outfile, curr_level, on->concurrent(i)->id);
      }

#else // ! IFMM_PARALLELIZE
      
      downward_pass_node(x_iFMM, outfile, curr_level, k);
      
#endif // ! IFMM_PARALLELIZE

    } // k

  }

}


void IFMM_Matrix::substitution(dvec & x_iFMM,ofstream& outfile){
  
  for (int curr_level=2; curr_level <= l;curr_level++){
    for (int k=lvl_idx(curr_level + 1)-1; k >= lvl_idx(curr_level); --k) {
      oct_node * on = lvl_p(k);
      set_xyz(on);
    }
  }
  
  //20200116    clock_t start_solve_top_level, stop_solve_top_level; double Time_solve_top_level=0;
  //20200116    clock_t start_downward_pass, stop_downward_pass; double Time_downward_pass=0;
  
  
  //20200116  start_solve_top_level =clock();
  TICK(timer_solve_top_level);
  solve_top_level();
  //20200116  stop_solve_top_level =clock();
  //20200116  Time_solve_top_level += double(stop_solve_top_level-start_solve_top_level)/double(CLOCKS_PER_SEC);
  TACK(timer_solve_top_level, outfile);
  
  //20200116  start_downward_pass=clock();
  TICK(timer_downward_pass);
  downward_pass(x_iFMM,outfile);
  //20200116  stop_downward_pass=clock();
  //20200116  Time_downward_pass += double(stop_downward_pass-start_downward_pass)/double(CLOCKS_PER_SEC);
  TACK(timer_downward_pass, outfile);
  
  //20200116  outfile << "Time_solve_top_level: " << Time_solve_top_level << endl;
  //20200116  outfile << "Time_downward_pass: " << Time_downward_pass << endl;
  
}


void IFMM_Matrix::eliminate_level_node_reuse(const int k)
{

  // CURRENT NODE
  oct_node * on = lvl_p(k);
  // printf("\n Eliminating node... %d \n", on->id);
  // printf("on->rank: %d\n", on->rank);

  // ELIMINATE CURRENT NODE
  on->IsEliminated=1;

  int iSelf=0; // index to current node (in the neighbour list of the current node)
  for (int iNgbr=0; iNgbr < on-> ngbr.m; iNgbr++){
    if (on->ngbr(iNgbr)->id==on->id){
      iSelf=iNgbr;
      break;
    }
  }

  // PRECOMPUTE CERTAIN MATRICES HERE
  // LU DECOMPOSITION OF P2P_ss
  // Eigen::PartialPivLU<MatrixXd> P2P_operator_self_LU = on->P2P_operator[iSelf].partialPivLu();
  MatrixXd rho_self=on->P2P_operator_self_LU.solve(on->L2P_operator);

  // MatrixXd alpha_self=MatrixXd::Zero(on->rank,on->rank);
  // alpha_self = - on->P2M_operator*rho_self; // EFF LU

  // LU DECOMPOSITION OF alpha_ii
  // Eigen::PartialPivLU<MatrixXd> alpha_self_LU = alpha_self.partialPivLu(); // EFF LU
  // MatrixXd alpha_self_inv=alpha_self_LU.solve(MatrixXd::Identity(on->rank,on->rank));

  VectorXd kappa_s = on->P2P_operator_self_LU.solve(on->RHS_leaf);
  VectorXd iota_s = on->alpha_self_LU.solve(on->P2M_operator*(kappa_s));

  vector<MatrixXd> beta_is(on->ngbr.m);
  vector<MatrixXd> mu_is(on->ngbr.m);

#if(0) // Tuning, N35, p10

  VectorXcd b_s = kappa_s + rho_self * iota_s;

  for (int iNgbr=0; iNgbr < on->ngbr.m; iNgbr ++) {
    if (iNgbr != iSelf) {
      if (on->ngbr(iNgbr)->IsEliminated) {
	on->ngbr(iNgbr)->RHS.noalias() += - on->P2L_operator[iNgbr] * b_s;
      } else {
	on->ngbr(iNgbr)->RHS_leaf.noalias() += - on->P2P_operator[iNgbr] * b_s;
      }
    }
  }

#else // original

  // LOOP OVER ALL NEIGHBOURS iNgbr
  for (int iNgbr=0; iNgbr < on-> ngbr.m; iNgbr++){
    if (on->ngbr(iNgbr)->IsEliminated){
      beta_is[iNgbr] = - on->P2L_operator[iNgbr]*rho_self;
    }
    else{
      mu_is[iNgbr] = - on->P2P_operator[iNgbr]*rho_self;
    }
  }

  for (int iNgbr=0; iNgbr < on-> ngbr.m; iNgbr++){
    if (on->ngbr(iNgbr)->IsEliminated){
      if (iNgbr==iSelf){
        // DO NOTHING YET, AS THIS WOULD AFFECT THE UPDATES OF THE OTHER NODES
      }
      else{
        on->ngbr(iNgbr)->RHS.noalias() += -on->P2L_operator[iNgbr]*kappa_s+beta_is[iNgbr]*iota_s;
      }
    }
    else{
      on->ngbr(iNgbr)->RHS_leaf.noalias() += - on->P2P_operator[iNgbr]*kappa_s+mu_is[iNgbr]*iota_s;
    }
  }

#endif // original

  on->RHS += - iota_s;

}

// void IFMM_Matrix::eliminate_level_reuse(int curr_level, bool IsLeaf, ofstream& outfile){
void IFMM_Matrix::eliminate_level_reuse(int curr_level, bool IsLeaf){

  // Loop over all nodes of the current level
  for (int k = lvl_idx(curr_level); k < lvl_idx(curr_level + 1); ++k) {
  
#if defined(IFMM_PARALLELIZE)
    
    oct_node *on = lvl_p(k);
    
#if defined(_OPENMP)
#pragma omp parallel for //num_threads(IFMM_PARALLELIZE)
#endif
    for (int i = 0; i < on->concurrent.m; i ++) {
      eliminate_level_node_reuse(on->concurrent(i)->id);
    }
    
#else // ! IFMM_PARALLELIZE
    
    eliminate_level_node_reuse(k);
    
#endif // ! IFMM_PARALLELIZE
    
  } // k

}


// void IFMM_Matrix::elimination_reuse(ofstream& outfile){
void IFMM_Matrix::elimination_reuse(){


  // // LEAF LEVEL
  // int iLevel=l;
  // bool IsLeaf=true;
  // for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
  //   oct_node * on = lvl_p(k);
  //   on->IsEliminated=0;
  //   // set RHS
  //   set_RHS(on,IsLeaf);
  // }
  // // NON-LEAF LEVELS (UP TO LEVEL 2)
  // for (int iLevel=l-1;iLevel>=2;iLevel--){
  //   IsLeaf=false;

  //   for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
  //     oct_node * on = lvl_p(k);
  //     on->IsEliminated=0;
  //     // set RHS
  //     set_RHS(on,IsLeaf);
  //   }
  // }


  for (int iLevel=l;iLevel>=2;iLevel--){
    // printf("Level %d\n",iLevel );
    // outfile << "Level: " << iLevel << endl;

    bool IsLeaf = (iLevel==l) ? true:false;

    // Loop over all levels
    eliminate_level_reuse(iLevel,IsLeaf);

    // PREPARE PARENT LEVEL
    if (iLevel>2){
      // bool isl_p=false;
#if defined(_OPENMP)
#pragma omp parallel for
#endif
      for (int k = lvl_idx(iLevel-1); k < lvl_idx(iLevel ); ++k) {
        oct_node * on = lvl_p(k);

        // transfer RHS
        transfer_RHS(on);
      }
    }
  }
}


// void set_P2M_operator(MatrixXd & P2M_operator,dmat & Sr_x, dmat & Sr_y, dmat & Sr_z)

// void set_M2M_operator(oct_node * node, MatrixXd & P2M_operator, int &n,dvec & Sr);

// void set_L2P_operator(MatrixXd & L2P_operator,dmat & Sr_x, dmat & Sr_y, dmat & Sr_z);

// void set_L2L_operator(oct_node * node, MatrixXd & L2P_operator, int &n,dvec & Sr);

// void set_M2L_operator(oct_node * node, int &n);

// void set_nDOF(oct_node * node, int &n, bool &IsLeaf);

// void set_xyz(oct_node * node);

// void set_P2L_operator(oct_node * node, int &n, bool &IsLeaf);

// void set_M2P_operator(oct_node * node, int &n, bool &IsLeaf);

// void set_RHS(oct_node * node, bool &IsLeaf);

// void set_P2P_operator(oct_node * node,kernel * kfun, bool &IsLeaf);

void IFMM_Matrix::create_lvl_indexing() { // ???
  // Number of nodes at each level
  lvl_idx.resize(l + 2);
  zero(lvl_idx);
  int curr_l = 0;
  nnode_at_level(false, curr_l, &root, lvl_idx, lvl_p);

  // Cumulative sum
  cum_sum(lvl_idx);

  int n_nodes = lvl_idx(lvl_idx.m - 1); // Total number of nodes
  lvl_p.resize(n_nodes);

  // Node IDs
  ivec nnode_(l + 1);
  for (int i = 0; i < lvl_idx.m - 1; ++i)
    nnode_(i) = lvl_idx(i);
  curr_l = 0;
  nnode_at_level(true, curr_l, &root, nnode_, lvl_p);

}


void allocate_M_L(int n3, Vector<oct_node*> & lvl_p) {
  // Allocate memory for M and L data
  for (int i = 0; i < lvl_p.m; ++i) {
    oct_node * on = lvl_p(i);
    on->M.resize(n3);
    on->L.resize(n3);
  }
}

// void IFMM_Matrix::setup_leaf_data(Vector<vec3> & point, dvec & b) {
void IFMM_Matrix::setup_leaf_data(Vector<vec3> & point) {

  // Initialize the vector with RHS for leaf nodes

  for (int i = lvl_idx(lvl_idx.m - 2); i < lvl_idx(lvl_idx.m - 1); ++i) {
    oct_node * on = lvl_p(i);


    list<int> & l = on->pidx;

    on->pnt.resize(l.size());
    on->sigma.resize(l.size());
    list<int>::iterator it = l.begin();

    int j = 0;
    for (; it != l.end(); ++it, ++j) {
      on->pnt(j) = point(*it);
      // on->sigma(j) = b(*it);
    }
  }
}

void IFMM_Matrix::setup_x_leaf_check(Vector<vec3> & point, dvec & x_leaf_check) {
  // Initialize the vector with RHS for leaf nodes

#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for (int i = lvl_idx(lvl_idx.m - 2); i < lvl_idx(lvl_idx.m - 1); ++i) {
    oct_node * on = lvl_p(i);

    list<int> & l = on->pidx;    

    on->pnt.resize(l.size());
    on->x_leaf_check.resize(l.size());
    list<int>::iterator it = l.begin();

    int j = 0;
    for (; it != l.end(); ++it, ++j) {
      on->pnt(j) = point(*it);
      on->x_leaf_check(j) = x_leaf_check(*it);
    }
  }
}

void IFMM_Matrix::init_M2L_p(ivec & Kidx, dvec & M2LOp) {
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

void IFMM_Matrix::setRHS(dvec & b){
#if defined(_OPENMP)
#pragma omp for
#endif
  for (int i = lvl_idx(lvl_idx.m - 2); i < lvl_idx(lvl_idx.m - 1); ++i) {
    oct_node * on = lvl_p(i);

    list<int> & l = on->pidx;
    
    // on->pnt.resize(l.size());
    // on->sigma.resize(l.size()*NDIM);
    on->sigma.resize(l.size());
    list<int>::iterator it = l.begin();

    int j = 0;
    for (; it != l.end(); ++it, ++j) {
      // on->pnt(j) = point(*it);
      // for (int iDim=0; iDim<NDIM; iDim++){
        // on->sigma(j*NDIM+iDim) = b((*it)*NDIM+iDim);
      // }
      on->sigma(j) = b(*it);
    }
  }

    // LEAF LEVEL
  int iLevel=l;
  bool IsLeaf=true;
#if defined(_OPENMP)
#pragma omp for
#endif
  for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
    oct_node * on = lvl_p(k);
    on->IsEliminated=0;
    // set RHS
    set_RHS(on,IsLeaf);
  }
  // NON-LEAF LEVELS (UP TO LEVEL 2)
  for (int iLevel=l-1;iLevel>=2;iLevel--){
    IsLeaf=false;

#if defined(_OPENMP)
#pragma omp for
#endif
    for (int k = lvl_idx(iLevel); k < lvl_idx(iLevel + 1); ++k) {
      oct_node * on = lvl_p(k);
      on->IsEliminated=0;
      // set RHS
      set_RHS(on,IsLeaf);
    }
  }

  
}

void IFMM_Matrix::initialization(Vector<vec3> & point, kernel * kfun,ofstream& outfile){


  //20200116  clock_t start_initialization, stop_initialization; double Time_initialization=0;
  //20200116  start_initialization =clock();    
  TICK(timer_initialization);

  // Calculate the center and size of the box containing all the points
  vec3 xmin(1e32, 1e32, 1e32);
  vec3 xmax(-1e32, -1e32, -1e32);
  for (int i = 0; i < point.m; ++i) {
    xmin.x = min(xmin.x, point(i).x);
    xmin.y = min(xmin.y, point(i).y);
    xmin.z = min(xmin.z, point(i).z);

    xmax.x = max(xmax.x, point(i).x);
    xmax.y = max(xmax.y, point(i).y);
    xmax.z = max(xmax.z, point(i).z);
  }

  double S;
  vec3 ctr;
  // if (Cubic_lattice==false){
    ctr = 0.5 * (xmin + xmax);
    S = norm_max(xmax - xmin);
  // }
  // else{
    // double Nx = pow(numSpheres,1./3.);
    // S =  length*pow(2,ceil(log2(Nx))); // for 3D regular lattice only!
    // ctr = 0.5 * (xmin + xmax);
    // ctr.x += (S-Nx*length)/2; // for 3D regular lattice only!
    // ctr.y += (S-Nx*length)/2; // for 3D regular lattice only!
    // ctr.z += (S-Nx*length)/2; // for 3D regular lattice only!
  // }



  // Distribute the sources and field points and set up the interaction and neighbor lists
  init_tree_Ilist(ctr, S, point);

  // Initialize the vectors of nodes at each level
  create_lvl_indexing();

  ivec Kidx(343);
  dvec M2LOp;

  // Initialize the M2L operator
  init_M2L_op(kfun, Kidx, S, M2LOp);

  // Initialize M2L pointers
  init_M2L_p(Kidx, M2LOp);

  // Copy data into leaf nodes
  setup_leaf_data(point);

  // Initiliaze IFMM operators
  initialize_ifmm_operators(kfun,outfile);

  //20200116  stop_initialization =clock();
  //20200116  Time_initialization += double(stop_initialization-start_initialization)/double(CLOCKS_PER_SEC);  
  //20200116  outfile << "Time_initialization: " << Time_initialization << endl;
  TACK(timer_initialization, outfile);

}


void ifmm(int l, int n, double epsilon, double epsilon_rel, int LR_mode, Vector<vec3> & point, dvec & b, kernel * kfun,dvec & x_iFMM, bool & Delay_update) {
// void ifmm(int l, int n, double epsilon, int LR_mode, Vector<vec3> & point, dvec & b, diag_plus_lr * dlr_fun,dvec & x_iFMM, dvec & xdirect_check) {

  // assert(l>0);
  // assert(n>0);
  // assert(point.m > 0);
  // assert(point.m == b.m);

  // const int N = point.m;

  // timeType t0, t1, t2, t3, t4; // Time variables

  // // // // int n3 = n * n * n; // n3 = n^3


  // // cout << "test" << endl;

  // // int Ntest=200;
  // // MatrixXd Atest = MatrixXd::Zero(Ntest,Ntest);
  // // for (int i=0;i<Ntest;i++){
  // //   for (int j=0;j<Ntest;j++){
  // //     Atest(i,j)=frand(0,10);
  // //   }
  // // }

  // // MatrixXd U, V;
  // // int ranktest=0;

  // // clock_t start = clock();
  // // ACA_FullyPivoted(Atest,U,V,epsilon,ranktest,Ntest);
  // // clock_t end = clock();    
  // // double Time  = double(end-start)/double(CLOCKS_PER_SEC);
  // // cout << "Time ACA: " << Time << endl;
  // // printf("ranktest: %d\n",ranktest);


  // // start = clock();
  // // JacobiSVD<MatrixXd> svd_Atest(Atest, ComputeThinU | ComputeThinV);
  // // end = clock();    
  // // Time  = double(end-start)/double(CLOCKS_PER_SEC);
  // // cout << "Time SVD: " << Time << endl;


  // ofstream outfile;
  // char filename [100];
  // // sprintf(filename,"output/benchmark_ii/2D_circle_rank_N_ %d _n_ %d _l_ %d _epsilon_%1.1e .txt", N,n,l,epsilon);
  // // sprintf(filename,"output/benchmark_ii/2D_circle_rank_test_SVD_N_ %d _n_ %d _l_ %d _epsilon_%1.1e .txt", N,n,l,epsilon);

  // // sprintf(filename,"output/benchmark_ii/2D_circle_test_rank_SVD_N_ %d _n_ %d _l_ %d _epsilon_%1.1e .txt", N,n,l,epsilon);
  // // sprintf(filename,"output/benchmark_ii/2D_circle_test_rank_ACA_N_ %d _n_ %d _l_ %d _epsilon_%1.1e .txt", N,n,l,epsilon);

  // // sprintf(filename,"output/benchmark_ii/2D_randompoints_test_rank_ACA_N_ %d _n_ %d _l_ %d _epsilon_%1.1e .txt", N,n,l,epsilon);
  // // sprintf(filename,"output/benchmark_ii/2D_randompoints_test_rank_SVD_N_ %d _n_ %d _l_ %d _epsilon_%1.1e .txt", N,n,l,epsilon);

  // // sprintf(filename,"output/benchmark_ii/3D_sphere_test_rank_SVD_N_ %d _n_ %d _l_ %d _epsilon_%1.1e .txt", N,n,l,epsilon);
  // // sprintf(filename,"output/benchmark_ii/3D_randompoints_test_rank_SVD_N_ %d _n_ %d _l_ %d _epsilon_%1.1e .txt", N,n,l,epsilon);

  // // sprintf(filename,"output/benchmark_ii/2D_circle_final_rank_SVD_N_ %d _n_ %d _l_ %d _epsilon_%1.1e .txt", N,n,l,epsilon);
  // // sprintf(filename,"output/benchmark_ii/2D_randompoints_final_rank_SVD_N_ %d _n_ %d _l_ %d _epsilon_%1.1e .txt", N,n,l,epsilon);

  // // sprintf(filename,"output/benchmark_ii/3D_randompoints_final_rank_SVD_N_ %d _n_ %d _l_ %d _epsilon_%1.1e .txt", N,n,l,epsilon);
  // // sprintf(filename,"output/benchmark_ii/3D_sphere_final_rank_SVD_N_%d_n_%d_l_%d_epsilon_%1.1e.txt", N,n,l,epsilon);

  // // sprintf(filename,"output/benchmark_ii/2D_randompoints_ThinKing_rank_RSVD_N_%d_n_%d_l_%d_epsilon_%1.1e_a_1e-6.txt", N,n,l,epsilon);
  // // sprintf(filename,"output/benchmark_ii/3D_randompoints_ThinKing_rank_RSVD_N_%d_n_%d_l_%d_epsilon_%1.1e_a_1e-4.txt", N,n,l,epsilon);

  // sprintf(filename,"output/diag_plus_poly/2D_randompoints_rank_RSVD_N_%d_n_%d_l_%d_epsilon_%1.1e_d_1e-3.txt", N,n,l,epsilon);
  // // sprintf(filename,"output/diag_plus_poly/3D_randompoints_rank_RSVD_N_%d_n_%d_l_%d_epsilon_%1.1e_d_1e3.txt", N,n,l,epsilon);



  





  
  // outfile.open(filename);

  // outfile << "N " << N << endl;
  // outfile << "n " << n << endl;
  // outfile << "l " << l << endl;
  // outfile << "epsilon " << epsilon << endl;
  // outfile << endl;

  // // Begin pre-computation timing
  // t0 = Timer();

  // // Calculate the center and size of the box containing all the points
  // vec3 xmin(1e32, 1e32, 1e32);
  // vec3 xmax(-1e32, -1e32, -1e32);
  // for (int i = 0; i < N; ++i) {
  //   xmin.x = min(xmin.x, point(i).x);
  //   xmin.y = min(xmin.y, point(i).y);
  //   xmin.z = min(xmin.z, point(i).z);

  //   xmax.x = max(xmax.x, point(i).x);
  //   xmax.y = max(xmax.y, point(i).y);
  //   xmax.z = max(xmax.z, point(i).z);
  // }

  // vec3 ctr = 0.5 * (xmin + xmax);
  // double S = norm_max(xmax - xmin);
  
  // IFMM_Matrix fmmd(l, n, epsilon, epsilon_rel, LR_mode);

  // // Distribute the sources and field points and set up the interaction and neighbor lists
  // fmmd.init_tree_Ilist(ctr, S, point);

  // // Initialize the vectors of nodes at each level
  // fmmd.create_lvl_indexing();

  // // Copy data into leaf nodes
  // // // // fmmd.setup_leaf_data(point, xdirect_check);
  

  // // Allocate memory for M and L
  // // // // allocate_M_L(n3, fmmd.lvl_p);

  // // Initialize the M2L operator
  // ivec Kidx(343);
  // dvec M2LOp;

  // // cout << "start" << endl;
  // fmmd.init_M2L_op(kfun, Kidx, S, M2LOp);
  // // fmmd.init_M2L_op(dlr_fun, Kidx, S, M2LOp);
  // // cout << "end" << endl;



  // // Initialize M2L pointers
  // fmmd.init_M2L_p(Kidx, M2LOp);

  // // End pre-computation timing
  // t1 = Timer();

  // // Compute the field using FMM
  // // // // fmmd.compute(kfun, x_iFMM);


  // // Copy data into leaf nodes
  // fmmd.setup_leaf_data(point, b);
  // // // // fmmd.setup_x_leaf_check(point, xdirect_check);


  // fmmd.initialize_ifmm_operators(kfun,outfile);
  // // fmmd.initialize_ifmm_operators(dlr_fun,outfile);
  // t2 = Timer();
  // fmmd.elimination(outfile);
  // t3 = Timer();
  // fmmd.substitution(x_iFMM,outfile);


  // // End FMM timing
  // t4 = Timer();

  // // printf("Pre-compution time: %.3fs,\n" "fmm time: %.3fs.\n" "ifmm time: %.3fs.\n", t1-t0, t2-t1,t3-t2);
  // printf("ifmm setup time: %.3fs.\n", t2-t0);
  // printf("ifmm elimination time: %.3fs.\n", t3-t2);
  // printf("ifmm substitution time: %.3fs.\n", t4-t3);
  // printf("ifmm total time: %.3fs.\n", t4-t0);

  // outfile << endl;
  
  // outfile << "ifmm setup time: " << t2-t0 << endl;
  // outfile << "ifmm elimination time: " << t3-t2 << endl;
  // outfile << "ifmm substitution time: " << t4-t3 << endl;
  // outfile << "ifmm total time: " << t4-t0 << endl;

  // outfile.close();

}
