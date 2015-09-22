#ifndef DIPOLE_MATRIX_HPP
#define DIPOLE_MATRIX_HPP

#include <string>
#include <vector>
#include "DipoleMatrix.hpp"
#include "Numerov.hpp"

typedef double real_t;

real_t radial_element(std::string species1, int n1, int l1, real_t j1, int power,
                      std::string species2, int n2, int l2, real_t j2);

real_t angular_element(int l1, real_t j1, real_t m1, int l2, real_t j2, real_t m2, real_t q);
#endif // DIPOLE_MATRIX_HPP
