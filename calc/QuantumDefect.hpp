#ifndef QUANTUM_DEFECT_HPP
#define QUANTUM_DEFECT_HPP

#include <string>
#include "dtypes.h"

//typedef double real_t;

class QuantumDefect {
private:
  real_t ac_;
  int Z_;
  real_t a1_, a2_, a3_, a4_;
  real_t rc_;
  real_t energy_;
public:
  QuantumDefect(std::string species, int n, int l, real_t j);
  const real_t &ac;
  const int &Z;
  const real_t &a1, &a2, &a3, &a4;
  const real_t &rc;
  const real_t &energy;
};

real_t energy_level(std::string species, int n, int l, real_t j);

#endif // QUANTUM_DEFECT_HPP
