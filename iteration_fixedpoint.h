#ifndef _ITERATION_FIXEDPOINT_
#define _ITERATION_FIXEDPOINT_

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


#include "container.h"
#include "ifmm.h"
#include "../BBFMM3D_pieter/include/bbfmm3d.hpp"
// #include "../BBFMM3D_pieter/include/kernel_Types.hpp"


using namespace std;
using namespace Eigen;


class myKernel: public H2_3D_Tree {
public:
    myKernel(double L, int level, int n,  double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n, epsilon, use_chebyshev){};
    virtual void setHomogen(string& kernelType,doft *dof) {
        homogen = 0;
        symmetry = 1;
        kernelType = "myKernel";
        dof->s = 1;
        dof->f = 1;
    }
    virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos, double *K, doft *dof) {
        vector3 diff;        
        // Compute 1/r
        diff.x = sourcepos.x - fieldpos.x;
        diff.y = sourcepos.y - fieldpos.y;
        diff.z = sourcepos.z - fieldpos.z;
        double r = sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
        double a=1e-2;


        if (r==0){
          *K = 1.0;
        }
        else{
          if(r<a){
            *K = r/a;
          }
          else{
            *K = a/r;
          }
        }
    }
};


// void iteration_fixedpoint(double &L, int &l, int &n_FMM, double &eps, int &m, int & use_chebyshev, vector3 * point_fmm, int &N, dvec &x, dvec &x0, dvec & b, int & max_iter, double & tol, IFMM_Matrix & ifmm, bool &precond);
void iteration_fixedpoint(H2_3D_Compute<myKernel> &compute, myKernel &Atree, vector3 * point_fmm, int &N, dvec &x, dvec &x0, dvec & b, int & max_iter, double & tol, IFMM_Matrix & ifmm, bool &precond);


#endif