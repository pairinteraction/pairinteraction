#ifndef NUMEROV_HPP
#define NUMEROV_HPP

#include <string>
#include <vector>
#include "QuantumDefect.hpp"

typedef double real_t;

class Numerov {
private:
  std::string species;
  int n, l, m, Z;
  QuantumDefect qd;
  int nsteps;
  real_t limit_inner;
  real_t limit_outer;
  real_t dx;
  std::vector<real_t> x;
  std::vector<real_t> y;

  real_t E();
  real_t V(real_t x);
  real_t g(real_t x);
  real_t step(int i);
  
public:
  Numerov(std::string species, int n, int l, int m, double j, int Z);
  std::vector<real_t> axis();
  std::vector<real_t> integrate();
};

#endif // NUMEROV_HPP
