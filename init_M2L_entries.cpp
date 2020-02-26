#include "init_M2L_entries.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;
class kernel;


void init_M2L_entries(int n, int l, double S, ivec & Kidx, kernel * kfun, dvec & M2LOp) {
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
          }
        }
      }
    }
    S *= 2.0;
  }
}

inline int M2L_index_map(vec3 & i) {
  return int(49 * (i.x + 3) + 7 * (i.y + 3) + i.z + 3);
}