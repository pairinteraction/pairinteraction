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


QuantumDefect::QuantumDefect(std::string species, int n, int l, real_t j)
  : ac(ac_), Z(Z_), a1(a1_), a2(a2_), a3(a3_), a4(a4_), rc(rc_), energy(energy_)
{
  std::stringstream ss;
  SQLite3 db("quantum_defects.db");
  int pot_max_l, ryd_max_l;
  double ryd_max_j;

  // Determine maximal L for model potentials
  ss.str(std::string());
  ss << "select MAX(L) from model_potential where (element = '" << species << "');";
  SQLite3Result res1 = db.query(ss.str().c_str());
  ss.str(std::string());
  if (res1.size() > 0)
    *res1.first() >> pot_max_l;
  else throw no_defect();

  // Determine maximal L for Rydberg-Ritz coefficients
  ss.str(std::string());
  ss << "select MAX(L) from rydberg_ritz where (element = '" << species << "');";
  SQLite3Result res2 = db.query(ss.str().c_str());
  ss.str(std::string());
  if (res2.size() > 0)
    *res2.first() >> ryd_max_l;
  else throw no_defect();

  // The l to be used is the minimum of the three below
  l = std::min({l, pot_max_l, ryd_max_l});


  // Determine maximal J for Rydberg-Ritz coefficients
  ss.str(std::string());
  ss << "select MAX(J) from rydberg_ritz where  ("
     << "(element = '" << species << "') "
     << "and (L = " << l << ") "
     << ");";
  SQLite3Result res3 = db.query(ss.str().c_str());
  ss.str(std::string());
  if (res3.size() > 0)
    *res3.first() >> ryd_max_j;
  else throw no_defect();

  // The j to be used is the minimum of the two below
  j = std::min({j, ryd_max_j});


  // Load model potentials from database
  ss << "select ac,Z,a1,a2,a3,a4,rc from model_potential where ("
     << "(element = '" << species << "') "
     << "and (L = " << l << ") "
     << ");";
  SQLite3Result res4 = db.query(ss.str().c_str());
  ss.str(std::string());
  if (res4.size() > 0)
    *res4.first() >> ac_ >> Z_ >> a1_ >> a2_ >> a3_ >> a4_ >> rc_;
  else throw no_defect();


  // Load Rydberg-Ritz coefficients from database
  ss.str(std::string());
  ss << "select d0,d2,d4,d6,d8,Ry from rydberg_ritz where ("
     << "(element = '" << species << "') "
     << "and (L = " << l << ") "
     << "and (J = " << j << ") "
     << ");";
  SQLite3Result res = db.query(ss.str().c_str());
  ss.str(std::string());
  real_t nstar = n;
  real_t Ry_inf = 109737.31568525;
  real_t d0, d2, d4, d6, d8, Ry = Ry_inf;
  if (res.size() > 0)
  {
    *res.first() >> d0 >> d2 >> d4 >> d6 >> d8 >> Ry;
    nstar -= d0 + d2/pow(n-d0,2) + d4/pow(n-d0,4)
      + d6/pow(n-d0,6) + d8/pow(n-d0,8);
  }
  else throw no_defect();

  energy_ = -.5*(Ry/Ry_inf)/(nstar*nstar);
}


real_t energy_level(std::string species, int n, int l, real_t j) {
    QuantumDefect qd(species, n, l, j);
    return qd.energy;
}
