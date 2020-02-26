#if defined(IFMM_PARALLELIZE)

#include "ifmm.h"

#define SEEK2_EMPTY ( - 1 )
#define SEEK2_EXIST (   0 )
#define SEEK2_CHECK (   1 )

void IFMM_Matrix::seek2_concurrent_nodes(const int curr_level, const int m)
{
  /* Create the groups of semi-concurrent nodes for a given
     degree of concurrency separation m */
  
  ///////////////
  ASSERT(m >= 3); // semi-concurrent degree
  ///////////////
  
  const int n = POW2(curr_level); // number of possible nodes at leaf level per dimension
  
  /* Create the map from the index of a node to its pointer and the status of the nodes */
  
  const oct_node *root = lvl_p(lvl_idx(0)); // root node

  vector<oct_node *> map(n * n * n, NULL); // NULL stands for empty
  vector<int> status(n * n * n, SEEK2_EMPTY);
  
  for (int k = lvl_idx(curr_level); k < lvl_idx(curr_level + 1); k ++) { // (non-empty) nodes in the leaf level
    oct_node *on = lvl_p(k);
    int i0 = (int)((on->ctr.x - (root->ctr.x - root->S * 0.5)) / on->S);
    int j0 = (int)((on->ctr.y - (root->ctr.y - root->S * 0.5)) / on->S);
    int k0 = (int)((on->ctr.z - (root->ctr.z - root->S * 0.5)) / on->S);
    //////////////////////////
    ASSERT(0 <= i0 && i0 < n);
    ASSERT(0 <= j0 && j0 < n);
    ASSERT(0 <= k0 && k0 < n);
    //////////////////////////
    int idx0 = i0 + n * (j0 + n * k0);
    map[idx0] = on;
    status[idx0] = SEEK2_EXIST;
    on->concurrent.resize(0); // necessary?
  }

  /* Create the groups */

  vector<vector<oct_node *>> group;

  for (int k0 = 0; k0 < n; k0 ++) {
    for (int j0 = 0; j0 < n; j0 ++) {
      for (int i0 = 0; i0 < n; i0 ++) {
	
	int idx0 = i0 + n * (j0 + n * k0);

	if (status[idx0] == SEEK2_EXIST) { // then, create a new group

	  vector<oct_node *> list;
	  list.push_back(map[idx0]);
	  status[idx0] = SEEK2_CHECK; // "eliminated" in the actual situation
	  bool jskip = true;
	  bool kskip = true;

	  int k = k0;

	  while (k < n) {

	    int j;
	    if (k == k0) {
	      j = j0;
	    } else {
	      j = 0;
	    }

	    while (j < n) {

	      int i;
	      if (j == j0 && k == k0) {
		i = i0 + m;
	      } else {
		i = 0;
	      }
	      
	      while (i < n) {

		int idx = i + n * (j + n * k);
		
		if (status[idx] == SEEK2_EXIST) {
		  list.push_back(map[idx]);
		  status[idx] = SEEK2_CHECK;
		  i += m;
		  jskip = true;
		  kskip = true;
		} else {
		  i ++;
		}

	      } // while(i<n)

	      if (jskip == true) {
		j += m;
		jskip = false;
	      } else {
		j ++;
	      }

	    } // while(j<n)

	    if (kskip == true) {
	      k += m;
	      kskip = false;
	    } else {
	      k ++;
	    }

	  } // while(k<n)

	  group.push_back(list);
	
	} 

      } // i0
    } // j0
  } // k0
  
#if defined(MYDEBUG)
  for (int k = 0; k < n; k ++) {
    for (int j = 0; j < n; j ++) {
      for (int i = 0; i < n; i ++) {
	int idx = i + n * (j + n * k);
	if (map[idx] != NULL) { // not empty but exist
	  ASSERT(status[idx] == SEEK2_CHECK); // all the non-empty nodes must be considered.
	}
      }
    }
  }
#endif
  
#if(0) // check
  INFO("group.size()=%zu\n", group.size());
  for (size_t igroup = 0; igroup < group.size(); igroup ++) {
    INFO("igroup=%zu: size=%zu:", igroup, group[igroup].size());
    for (size_t k = 0; k < group[igroup].size(); k ++) {
      fprintf(stderr, " %d,", group[igroup][k]->icell);
    }
    fprintf(stderr, "\n");
  }
#endif
  
  /* Copy the lists to the nodes */
  
  for (size_t igroup = 0; igroup < group.size(); igroup ++) {
    oct_node *on = group[igroup][0]; // the first node in this group is chosen as the representative.
    ///////////////////
    ASSERT(on != NULL);
    ///////////////////
    on->concurrent.resize(group[igroup].size());
    for (size_t k = 0; k < group[igroup].size(); k ++) {
      on->concurrent(k) = group[igroup][k];
    }
  }
  
}

void IFMM_Matrix::statistics_concurrent_nodes()
{
  for (int curr_level = 2; curr_level <= l; curr_level ++) {

    const int nnodes = lvl_idx(curr_level + 1) - lvl_idx(curr_level); // number of nodes in this level
    
    int *count = new int [nnodes + 1];

    for (int i = 0; i <= nnodes; i ++) {
      count[i] = 0; // number of nodes that can be processed concurrently with other i-1 nodes, i.e., nodes of concurrency of i.
    }

    for (int k = lvl_idx(curr_level); k < lvl_idx(curr_level + 1); k ++) {
      oct_node *on = lvl_p(k);
      count[on->concurrent.m] ++;
      ///////////////////////////////////////////////////////////////////////////////
      DBG("curr_level=%d k=%d: concurrent.m=%d\n", curr_level, k, on->concurrent.m);
      ///////////////////////////////////////////////////////////////////////////////
    }

    double average = 0;
    int minimum = nnodes + 1;
    int maximum = 0;
    int sum = 0;
    double nelim = 0; // number of eliminations at this level in parallel algorithm
    //int num_threads = omp_get_max_threads();
    int num_threads = IFMM_PARALLELIZE;
    for (int i = 1; i <= nnodes; i ++) { // exclude the nodes of concurrency of zero, as they are not actually computed.
      average += count[i] * i;
      sum += count[i];
      if (count[i] > 0) {
	if (i < minimum) minimum = i;
	if (i > maximum) maximum = i;
      }
      nelim += ceil((double)i / num_threads) * count[i];
    }
    average /= sum;

    double sp = nnodes / nelim; // theoretical speedup

    //170902    INFO("curr_level=%d nnodes=%d sum=%d average=%f minimum=%d maximum=%d\n", curr_level, nnodes, sum, average, minimum, maximum);
    INFO("level=%d nnodes=%d sum=%d ave=%f min=%d max=%d nelim=%f sp=%f\n", curr_level, nnodes, sum, average, minimum, maximum, nelim, sp);

    char str[256];
    sprintf(str, "concurrency_histogram_level%d.txt", curr_level);
    FILE *fp = fopen(str, "w");
    for (int i = 1; i <= lvl_idx(curr_level + 1) - lvl_idx(curr_level); i ++) {
      fprintf(fp, "%d %d\n", i, count[i]);
    }
    INFO("Created %s.\n", str);
    fclose(fp);

    delete [] count;

  } // curr_level

}

#endif // IFMM_PARALLELIZE
