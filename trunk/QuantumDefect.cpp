#include "QuantumDefect.hpp"
#include "SQLite.hpp"

#include <string>
#include <iostream>
#include <sstream>
#include <exception>
#include <cmath>


struct no_defect : public std::exception {
  const char* what () const throw ()
  {
    return "There is no defect available";
  }
};


void QuantumDefect::H(int n) {
  ac_ = 0.0;
  Z_  = 1;
  a1_ = 0.0;
  a2_ = 0.0;
  a3_ = 0.0;
  a4_ = 0.0;
  rc_ = 1.0;
  energy_ = -.5/(n*n);
}


void QuantumDefect::Rb(int l /*, real_t j*/) {
  //***************
  //* Rubidium 87 *
  //***************
  ac_ = 9.076;
  Z_ = 37;
  switch (l) {
  case 0: // S-state
    a1_ = 3.69628474;
    a2_ = 1.64915255;
    a3_ = -9.86069196;
    a4_ = 0.19579987;
    rc_ = 1.66242117;
    break;
  case 1: // P-state
    a1_ = 4.44088978;
    a2_ = 1.92828831;
    a3_ = -16.79597770;
    a4_ = -0.81633314;
    rc_ = 1.50195124;
    break;
  case 2: // D-state
    a1_ = 3.78717363;
    a2_ = 1.57027864;
    a3_ = -11.6558897;
    a4_ = 0.52942835;
    rc_ = 4.86851938;
    break;
  default: // L>= 3
    a1_ = 2.39848933;
    a2_ = 1.76810544;
    a3_ = -12.0710678;
    a4_ = 0.77256589;
    rc_ = 4.79831327;
    break;
  }
}


QuantumDefect::QuantumDefect(std::string species, int n, int l, real_t j)
  : ac(ac_), Z(Z_), a1(a1_), a2(a2_), a3(a3_), a4(a4_), rc(rc_), energy(energy_)
{
  if (species == std::string("Rb"))
    Rb(l/*, j*/);
  else if (species == std::string("H"))
    H(n);
  else throw no_defect();

  std::stringstream ss;
  ss << "select d0,d2,d4,d6,d8,Ry from quantum_defects where ("
     << "(element = '" << species << "') "
     << "and (L = " << l << ") "
     << "and (J = " << j << ") "
     << ");";
  SQLite3 db("quantum_defects.db");
  SQLite3Result res = db.query(ss.str().c_str());

  ss.str(std::string());
  real_t nstar = n;
  real_t Ry_inf = 109737.31568525;
  real_t d0, d2, d4, d6, d8, Ry = Ry_inf;
  if (res.size() > 0) {
    res.first() >> d0 >> d2 >> d4 >> d6 >> d8 >> Ry;
    nstar -= d0 + d2/pow(n-d0,2) + d4/pow(n-d0,4)
      + d6/pow(n-d0,6) + d8/pow(n-d0,8);
  }
  energy_ = -.5*(Ry/Ry_inf)/(nstar*nstar);
}

real_t energy_level(std::string species, int n, int l, real_t j) {
    QuantumDefect qd(species, n, l, j);
    return qd.energy;
}
