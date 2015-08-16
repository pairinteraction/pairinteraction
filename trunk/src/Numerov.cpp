#include "Numerov.hpp"
#include "utils.hpp"
#include "QuantumDefect.hpp"


#include <iostream>
void vprint(std::vector<real_t> x) {
  for (std::vector<real_t>::iterator it = x.begin(); it != x.end(); it++)
    std::cout << *it << ' ';
  std::cout << std::endl;
}

#include <string>
#include <complex>
#include <vector>
#include <cmath>

//#include <boost/math/special_functions/spherical_harmonic.hpp>

Numerov::Numerov(std::string species, int n, int l, double j)
  : species(species), n(n), l(l), j(j), qd(species,n,l,j) {
  nsteps = 1000;
  limit_inner = 2.0e-5;
  limit_outer = 2.0e+2;
  dx = (limit_outer - limit_inner)/(nsteps-1);
  x = utils::linspace(limit_outer, limit_inner, nsteps);
  y.assign(nsteps,0.0);
}


std::vector<real_t> Numerov::axis() {
  return x;
}


std::vector<real_t> Numerov::integrate() {
  // Set the initial condition
  if ( (n-l) % 2 == 0 )
    y[1] = -1e-10;
  else
    y[1] = 1e-10;

  // Perform the integration using Numerov's scheme
  for (int i = 2; i < nsteps; i++)
    y[i] = step(i);

  // Normalization
  for (int i = 0; i < nsteps; i++)
    y[i] /= x[i];

  real_t norm = 0;
  for (int i = 0; i < nsteps; i++)
    norm += y[i]*y[i] * x[i]*x[i] * dx;
  norm = sqrt(norm);

  for (int i = 0; i < nsteps; i++)
    y[i] /= norm;
  
  return y;
}


real_t Numerov::V(real_t x) {
  real_t Z_l = 1 + (qd.Z - 1)*exp(-qd.a1*x)
    - x*(qd.a3 + qd.a4*x)*exp(-qd.a2*x);
  real_t V_c = -Z_l / x;
  real_t V_p = -qd.ac/(2.*x*x*x*x) * (1-exp(-pow(x/qd.rc,6)));
  real_t V_so = 0.0;
  if ( l < 4 ) {
    real_t alpha = 7.29735257e-3;// ~1/137 fine-structure constant
    V_so = alpha*alpha/(4 * x*x*x) * (j*(j+1) - l*(l+1) - 0.5*(1+0.5));
  }
  return V_c + V_p + V_so;
}


real_t Numerov::g(real_t x) {
  return l*(l+1)/(x*x) + 2.*(V(x) - qd.energy);
  //return (2.*l+.5)*(2.*l+1.5)/(x*x) + 8.*(x*x)*(V(x*x) - qd.energy);
}


real_t Numerov::step(int i) {
  real_t A = (2. + 5./6. * dx*dx * g(x[i-1])) * y[i-1];
  real_t B = (1. - 1./12.* dx*dx * g(x[i-2])) * y[i-2];
  real_t C =  1. - 1./12.* dx*dx * g(x[i]);
  return (A - B)/C  ;
}

