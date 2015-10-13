#ifndef NUMEROV_HPP
#define NUMEROV_HPP

#include <string>
#include <vector>
#include "dtypes.h"
#include "QuantumDefect.hpp"

//typedef double real_t;

class Numerov {
private:
  std::string species;
  int n, l;
  real_t j;
  QuantumDefect qd;
  int nsteps_;
  real_t xmin_;
  real_t xmax_;
  real_t dx_;
  std::vector<real_t> x;
  std::vector<real_t> y;

  real_t V(real_t x);
  real_t g(real_t x);
  real_t step(int i);
  
public:
  const int &nsteps;
  const real_t &xmin;
  const real_t &xmax;
  const real_t &dx;
  Numerov(std::string species, int n, int l, real_t j);
  std::vector<real_t> axis();
  std::vector<real_t> integrate();
};

#endif // NUMEROV_HPP
