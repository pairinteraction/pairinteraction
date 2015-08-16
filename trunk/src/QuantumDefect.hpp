#ifndef QUANTUM_DEFECT_HPP
#define QUANTUM_DEFECT_HPP

#include <string>

class QuantumDefect {
private:
  double ac_;
  int Z_;
  double a1_, a2_, a3_, a4_;
  double rc_;
  double energy_;
  void H(int n);
  void Rb87(int n, int l, double j);
public:
  QuantumDefect(std::string species, int n, int l, double j);
  const double &ac;
  const int &Z;
  const double &a1, &a2, &a3, &a4;
  const double &rc;
  const double &energy;
};

#endif // QUANTUM_DEFECT_HPP
