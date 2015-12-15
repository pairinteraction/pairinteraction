#ifndef CONF_PARSER_HPP
#define CONF_PARSER_HPP

#include <string>
#include "dtypes.h"

class Configuration {
public:
  int load_from_json(std::string filename);
  std::string element;
  int n1;
  int n2;
  int l1;
  int l2;
  real_t j1;
  real_t j2;
  real_t m1;
  real_t m2;
  std::string symmetry;
  real_t angle;
  bool enable_dipdip;
  bool enable_dipquad;
  bool enable_quadquad;
  real_t efield_strength;
  bool efield_increasing;
  real_t bfield_strength;
  bool bfield_increasing;
  int delta_n;
  int delta_l;
  int delta_m;
  real_t delta_energy;
  bool preserve_M;
  bool preserve_submatrix;
  bool preserve_parityL;
  real_t distance_min;
  real_t distance_max;
  int distance_steps;
  real_t energyrange;
/*
 diagonalizer;
 cutoff_eigenVec;
 cutoff_interactionMat;
 slepc_numeigenpairs;
 slepc_solvertollerance;
 slepc_maxiteration;
 slp_blockdimension;
 slp_enablesquaregrid;
*/
};


#endif // CONF_PARSER_HPP
