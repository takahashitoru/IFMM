#include "create_FMM_tree.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

void create_FMM_tree(int levels, int idx, vec3 & point, oct_node * node) {
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
}