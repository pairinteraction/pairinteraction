#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

#include <string>
#include <vector>
#include "dtypes.h"
#include "QuantumDefect.hpp"

// --- Common interface ---

class Wavefunction {
protected:
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
  void initialize(double xmin);
public:
  const int &nsteps;
  const real_t &xmin;
  const real_t &xmax;
  const real_t &dx;
  Wavefunction(std::string species, int n, int l, real_t j)
    : species(species), n(n), l(l), j(j), qd(species,n,l,j),
      nsteps(nsteps_), xmin(xmin_), xmax(xmax_), dx(dx_) {};
  std::vector<real_t> axis();
  std::vector<real_t> integrate();
};

// --- Numerov's method ---

class Numerov : public Wavefunction {
private:
  real_t V(real_t x);
  real_t g(real_t x);
  real_t step(int i);
public:
  Numerov(std::string species, int n, int l, real_t j);
  std::vector<real_t> axis();
  std::vector<real_t> integrate();
};

// --- Whittaker method ---

class Whittaker : public Wavefunction {
public:
  Whittaker(std::string species, int n, int l, real_t j);
  std::vector<real_t> axis();
  std::vector<real_t> integrate();
};

// --- Matrix element calculation ---

size_t findidx(std::vector<real_t> x, real_t d);


template<typename T>
real_t IntegrateRadialElement(T N1, int power, T N2) {
    std::vector<real_t> x1 = N1.axis();
    std::vector<real_t> y1 = N1.integrate();
    std::vector<real_t> x2 = N2.axis();
    std::vector<real_t> y2 = N2.integrate();

    real_t xmin = N1.xmin >= N2.xmin ? N1.xmin : N2.xmin;
    real_t xmax = N1.xmax <= N2.xmax ? N1.xmax : N2.xmax;

    real_t mu = 0;
    // If there is an overlap, calculate the matrix element
    if (xmin <= xmax) {
        int start1 = findidx(x1, xmin);
        int end1   = findidx(x1, xmax);
        int start2 = findidx(x2, xmin);
        int end2   = findidx(x2, xmax);

        int i1, i2;
        for (i1 = start1, i2 = start2; i1 < end1 && i2 < end2; i1++, i2++) {
            mu += y1[i1]*y2[i2] * pow(x1[i1], 2*power+2) * N1.dx;
        }
        mu = fabs(2*mu);
    }

    return mu;
}

#endif // WAVEFUNCTION_HPP
