#ifndef QUANTUM_DEFECT_HPP
#define QUANTUM_DEFECT_HPP

#include <string>

class QuantumDefect {
private:
  double ac_;
  int Z_;
  double a1_, a2_, a3_, a4_;
  double rc_;
  double d0, d2, d4, d6, d8;
  double delta_;
  void Rb87(int n, int l, double j);
  void Rb85(int n, int l, double j);
  void Cs(int n, int l, double j);
  void calc_delta(int n);
public:
  QuantumDefect(std::string species, int n, int l, double j);
  const double& delta;
  const double& ac;
  const int& Z;
  const double& a1;
  const double& a2;
  const double& a3;
  const double& a4;
  const double& rc;
};

#endif // QUANTUM_DEFECT_HPP
