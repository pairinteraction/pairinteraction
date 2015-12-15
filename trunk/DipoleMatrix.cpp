#include "Numerov.hpp"
#include "DipoleMatrix.hpp"
#include "SQLite.hpp"

#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <gsl/gsl_sf_coupling.h>

real_t wigner3j(real_t ja, real_t jb, real_t jc,
                real_t ma, real_t mb, real_t mc) {
  return gsl_sf_coupling_3j(2*ja, 2*jb, 2*jc,
                            2*ma, 2*mb, 2*mc);
}

real_t wigner6j(real_t ja, real_t jb, real_t jc,
                real_t ma, real_t mb, real_t mc) {
  return gsl_sf_coupling_6j(2*ja, 2*jb, 2*jc,
                            2*ma, 2*mb, 2*mc);
}

size_t find(std::vector<real_t> x, real_t d) {
  size_t i;
  for (i = 0; i < x.size(); ++i) {
    if (x[i] == d)
      break;
  }
  return i;
}


real_t radial_element_cached(std::string species, int n1, int l1, real_t j1,
                             int power,
                             int n2, int l2, real_t j2) {
  real_t mu = 0;
  std::stringstream ss;

  SQLite3 db("radial_elements.db");
  ss << "select mu "
     << "from radial_elements where ("
     << "(species = '" << species << "') "
     << "and (n1 = " << n1 << ") "
     << "and (l1 = " << l1 << ") "
     << "and (j1 = " << j1 << ") "
     << "and (power = " << power << ") "
     << "and (n2 = " << n2 << ") "
     << "and (l2 = " << l2 << ") "
     << "and (j2 = " << j2 << ") "
     << ");";

  SQLite3Result res = db.query(ss.str().c_str());
  ss.str(std::string());
  if (res.nRow > 0) {
    printf("Element found!\n");
    ss << res.azResult[res.nRow*res.nColumn + res.nColumn-1];
    ss >> mu;
  }
  printf("Element not found!\n");

  mu = radial_element(species, n1, l1, j1, power, n2, l2, j2);

  ss.str(std::string());
  ss << "insert or ignore into radial_elements "
     << "(species,n1,l1,j1,power,n2,l2,j2,mu)"
     << "values("
     << "'" << species << "'," << n1 << "," << l1 << "," << j1 << ","
     << power << ","
     << n2 << "," << l2 << "," << j2 << ","
     << mu << ");";
  db.query(ss.str().c_str());

  return mu;
}


real_t radial_element(std::string species, int n1, int l1, real_t j1, int power,
                      int n2, int l2, real_t j2) {
  Numerov N1(species, n1, l1, j1);
  Numerov N2(species, n2, l2, j2);

  std::vector<real_t> x1 = N1.axis();
  std::vector<real_t> y1 = N1.integrate();
  std::vector<real_t> x2 = N2.axis();
  std::vector<real_t> y2 = N2.integrate();

  real_t xmin = N1.xmin >= N2.xmin ? N1.xmin : N2.xmin;
  real_t xmax = N1.xmax <= N2.xmax ? N1.xmax : N2.xmax;

  real_t mu = 0;
  // If there is an overlap, calculate the matrix element
  if (xmin <= xmax) {
    int start1 = find(x1, xmin);
    int end1   = find(x1, xmax);
    int start2 = find(x2, xmin);
    int end2   = find(x2, xmax);

    int i1, i2;
    for (i1 = start1, i2 = start2; i1 < end1 && i2 < end2; i1++, i2++) {
      mu += y1[i1]*y2[i2] * pow(x1[i1], 2*power+2) * N1.dx;
    }
    mu = fabs(2*mu);
  }
  
  return mu;
}

real_t angular_L_basis(int l1, int l2) {
  real_t A = pow( -1, l1 );
  real_t B = sqrt((2*l1+1)*(2*l2+1));
  real_t C = wigner3j(l1,1,l2 , 0,0,0);
  return A*B*C;
}

real_t angular_J_basis(int l1, real_t j1, real_t m1, int l2, real_t j2, real_t m2, real_t q) {
  real_t s = 0.5;
  real_t A = pow( -1, roundf(j1-m1 + l1+s+j2+1) );
  real_t B = sqrt((2*j1+1)*(2*j2+1));
  real_t C = wigner3j(j1,1,j2 , -m1,q,m2);
  real_t D = wigner6j(l1,j1,s , j2,l2,1);
  return A*B*C*D;
}

real_t angular_element(int l1, real_t j1, real_t m1, int l2, real_t j2, real_t m2, real_t q) {
  return angular_L_basis(l1, l2) * angular_J_basis(l1, j1, m1, l2, j2, m2, q);
}

bool selection_rules(int l1, real_t j1, real_t m1, int l2, real_t j2, real_t m2) {
  int delta_l = abs(l1-l2);
  real_t delta_j = fabs(j1-j2);
  real_t delta_m = fabs(m1-m2);
  bool select_l = ( delta_l == 1 );
  bool select_j = ( delta_j == 0 || delta_j == 1 );
  bool select_m = ( delta_m == 0 || delta_m == 1 );
  return ( select_l && select_j && select_m );
}
