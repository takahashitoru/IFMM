#include "init_M2L_block.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

void init_M2L_block(double b_size, dvec & nodes, vec3 & c_f, vec3 & c_s, dmat & Kmat, kernel * kfun) {
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