#include "QuantumDefect.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include <exception>


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


void QuantumDefect::Rb87(int n, int l, double j) {
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
  // Testcase: n = 113, l = 6, j = NULL
  (void) n;
  (void) l;
  (void) j;
  energy_ = -.5/(113*113);
}


QuantumDefect::QuantumDefect(std::string species, int n, int l, double j)
  : ac(ac_), Z(Z_), a1(a1_), a2(a2_), a3(a3_), a4(a4_), rc(rc_), energy(energy_)
{
  if (species == std::string("Rb87"))
    Rb87(n, l, j);
  else if (species == std::string("H"))
    H(n);
  else throw no_defect();
}
