/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * This file is part of the pairinteraction library.
 *
 * The pairinteraction library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The pairinteraction library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

constexpr const double au2GHz = 6579683.920757349;
constexpr const double au2Vcm = 5142206707.0;
constexpr const double au2G = 2350517550.0;
constexpr const double au2um = 5.2917721067e-05;
constexpr const double inverse_electric_constant = au2Vcm * au2Vcm * au2um * au2um * au2um / au2GHz;
constexpr const double sqrt_inverse_electric_constant =
    0.7717003798774048; // equals sqrt(inverse_electric_constant)
constexpr const double inverse_electron_rest_mass = au2Vcm * au2Vcm / au2G / au2G / au2GHz;

constexpr const double coulombs_constant = au2Vcm * au2Vcm * au2um * au2um * au2um / au2GHz;
constexpr const double electron_rest_mass = au2G * au2G * au2GHz / au2Vcm / au2Vcm;
constexpr const double elementary_charge = au2GHz / au2um / au2Vcm;
constexpr const double bohr_magneton = 0.5 * au2GHz / au2G;
constexpr const double reduced_planck_constant = au2GHz * au2um * au2G / au2Vcm;
constexpr const double speed_of_light = 137.035999139 * au2Vcm / au2G;

constexpr const double muB = 0.5; // in atomic units
constexpr const double gS = 2.0023192;
constexpr const double gL = 1;

constexpr const int ARB = 32767;

#ifdef WITH_INTEL_MKL
constexpr const bool mkl_enabled = true;
#else  // WITH_INTEL_MKL
constexpr const bool mkl_enabled = false;
#endif // WITH_INTEL_MKL

#ifdef WITH_GSL
constexpr const bool gsl_enabled = true;
#else  // WITH_GSL
constexpr const bool gsl_enabled = false;
#endif // WITH_GSL

#endif
