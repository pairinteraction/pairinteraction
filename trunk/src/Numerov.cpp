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

Numerov::Numerov(std::string species, int n, int l, int m, double j, int Z)
  : species(species), n(n), l(l), m(m), Z(Z), qd(species,n,l,j) {
  this->nsteps = 10000;
  this->limit_inner = 2.0e-5;
  this->limit_outer = 2.0e+2;
  this->dx = (this->limit_outer - this->limit_inner)/(this->nsteps-1);
  this->x = utils::linspace(this->limit_outer, this->limit_inner, this->nsteps);
  this->y.assign(this->nsteps,0.0);
}

std::vector<real_t> Numerov::axis() {
  return this->x;
}

std::vector<real_t> Numerov::integrate() {
  // Set the initial condition
  if ( (this->n-this->l) % 2 == 0 )
    this->y[1] = -1e-10;
  else
    this->y[1] = 1e-10;

  // Perform the integration using Numerov's scheme
  for (int i = 2; i < this->nsteps; i++)
    this->y[i] = step(i);

  // Normalisation
  for (int i = 0; i < this->nsteps; i++)
    this->y[i] /= this->x[i];

  real_t norm = 0;
  for (int i = 0; i < this->nsteps; i++)
    norm += this->y[i]*this->y[i] * this->x[i]*this->x[i] * this->dx;
  norm = sqrt(norm);

  for (int i = 0; i < this->nsteps; i++)
    this->y[i] /= norm;
  
  return y;
}

real_t Numerov::E() {
  return -.5/(this->n*this->n);
}

real_t Numerov::V(real_t x) {
  return -this->Z / x;
}

real_t Numerov::g(real_t x) {
  return 2.*(this->E() - this->V(x)) - this->l*(this->l+1)/(x*x);
}

real_t Numerov::step(int i) {
  real_t A = (2. - 5./6. * this->dx*this->dx * this->g(this->x[i-1])) * this->y[i-1];
  real_t B = (1. + 1./12.* this->dx*this->dx * this->g(this->x[i-2])) * this->y[i-2];
  real_t C =  1. + 1./12.* this->dx*this->dx * this->g(this->x[i]);
  return (A - B)/C  ;
}
