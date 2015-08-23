#include "Numerov.hpp"
#include "DipoleMatrix.hpp"

#include <string>
#include <vector>
#include <cmath>

//#include <boost/math/special_functions/spherical_harmonic.hpp>

size_t find(std::vector<real_t> x, real_t d) {
  size_t i;
  for (i = 0; i < x.size(); ++i) {
    if (x[i] == d)
      break;
  }
  return i;
}

real_t radial_element(std::string species1, int n1, int l1, real_t j1, int power,
                      std::string species2, int n2, int l2, real_t j2) {
  Numerov N1(species1, n1, l1, j1);
  Numerov N2(species2, n2, l2, j2);

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
