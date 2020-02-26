#include "init_M2L_entries.h"
#include "ifmm.h"
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

void set_points(double L, Vector<vec3> & point) {

  // cout << "sphere " << endl;
  // cout << "circle " << endl;

  // double pi = M_PI;
  // double phi, theta;


  for (int i = 0; i < point.m; i++) {
    // 3D SPHERE
    // theta = frand(0, 1)*2*pi;  // RAND
    // phi = acos(2*frand(0, 1)-1);
    

    // point(i).x = L * cos(theta) * sin(phi);
    // point(i).y = L * sin(theta) * sin(phi);
    // point(i).z = L * cos(phi);


    // RANDOM POINTS IN A CUBE
    point(i).x = frand(-1, 1) * L;
    point(i).y = frand(-1, 1) * L;
    point(i).z = frand(-1, 1) * L;

    

  }
}