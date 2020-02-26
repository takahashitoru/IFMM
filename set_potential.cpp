#include "init_M2L_entries.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

void set_potential(double L, dvec & b, Vector<vec3> & point) {

  // cout << "sphere " << endl;
  // cout << "circle " << endl;

  // double pi = M_PI;
  // double phi, theta;


  for (int i = 0; i < point.m; i++) {
    // SET RHS
    // b(i) = frand(-1, 1);
    b(i) = 1;
    // b(i) = (double) (i+1)/point.m;

    // printf("b(i)%f\n",b(i));


    // 3D SPHERE
    // theta = frand(0, 1)*2*pi;  // RAND
    // phi = acos(2*frand(0, 1)-1);
    
    // 2D CIRCLE
    // phi = pi/2; 
    // // // theta = frand(0, 1)*2*pi;  // RAND
    // theta = ((double) i/point.m)*2*pi; // UNIFORM

    
    // point(i).x = 0.95*L * cos(theta) * sin(phi);
    // point(i).y = 0.95*L * sin(theta) * sin(phi);
    // point(i).z = 0.95*L * cos(phi);

    // printf("point(i).x, point(i).y, point(i).z : %f %f %f  ; \n",point(i).x, point(i).y, point(i).z);

    // RANDOM POINTS IN A CUBE
    // point(i).x = frand(-0.5, 0.5) * L;
    // point(i).y = frand(-0.5, 0.5) * L;
    // point(i).z = frand(-0.5, 0.5) * L;

    point(i).x = frand(-1, 1) * L;
    point(i).y = frand(-1, 1) * L;
    point(i).z = frand(-1, 1) * L;

    // 2D/1D
    // point(i).x = 0.95*L;
    // point(i).y = 0.95*L;
    point(i).z = 0;

    // 1D
    // point(i).x = (((double) i/point.m)-0.5)*L; // UNIFORM
    // point(i).y = 0;
    // point(i).z = 0;

    // printf("point(i).x, point(i).y, point(i).z : %f %f %f  ; \n",point(i).x, point(i).y, point(i).z);


  }
}