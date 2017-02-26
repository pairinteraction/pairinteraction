/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef QUANTUM_DEFECT_H
#define QUANTUM_DEFECT_H

#include <string>
#include "dtypes.h"

//typedef double real_t;

/** \brief Quantum defect storage
 *
 * The quantum defect queries a database for Rydberg-Ritz coefficients and
 * model potential parameters. The Rydberg-Ritz coefficients are not exposed to
 * the outside but included in the energy data member.
 *
 * The slight deviation from the hydrogen atom for Rydberg atoms are captured
 * in the so-called quantum defect. The hydrogen energy is given by \f$ E =
 * \mathrm{Ry}/n^2 \f$. For Rydberg atoms this is slightly modified to \f$ E^*
 * = \mathrm{Ry}/n^{*2} \f$ where \f$ n^* = n - \delta_{nlj} \f$. This
 * correction can be expanded in a series
 * \f[
 *      \delta_{nlj} = \delta_0
 *      + \frac{\delta_2}{(n-\delta_0)^2}
 *      + \frac{\delta_4}{(n-\delta_0)^4}
 *      + \frac{\delta_6}{(n-\delta_0)^6}
 *      + \cdots \;
 * \f]
 * The coefficients of the series \f$ \delta_i \f$ (Rydberg-Ritz coefficients)
 * are loaded from the database.
 *
 * The object also loads the parameters for the model potentials which are used
 * in the Numerov object.
 */
class QuantumDefect {
private:
    real_t ac_;
    int Z_;
    real_t a1_, a2_, a3_, a4_;
    real_t rc_;
    real_t nstar_;
    real_t energy_;
public:
    /** \brief Constructor
     *
     * Save the input and query the database.
     *
     * \param[in] species   atomic species of the atom
     * \param[in] n         principal quantum number
     * \param[in] l         angular quantum number
     * \param[in] j         magnetic quantum number
     */
    QuantumDefect(std::string const& species, int n, int l, real_t j);

    /** \brief Atomic species */
    const std::string species;

    /** \brief Principal quantum number */
    const int n;

    /** \brief Angular quantum number */
    const int l;

    /** \brief Magnetic quantum number */
    const real_t j;

    /** \brief Polarizability */
    const real_t &ac;

    /** \brief Core charge */
    const int &Z;

    /** \brief Parameter of effective Coulomb potential */
    const real_t &a1;

    /** \brief Parameter of effective Coulomb potential */
    const real_t &a2;

    /** \brief Parameter of effective Coulomb potential */
    const real_t &a3;

    /** \brief Parameter of effective Coulomb potential */
    const real_t &a4;

    /** \brief Effective core size */
    const real_t &rc;

    /** \brief Effective principal quantum number */
    const real_t &nstar;

    /** \brief State energy */
    const real_t &energy;
};


/** \brief Returns only the energy of a state
 *
 * \warning This function has a considerable overhead because it allocates a new
 * quantum defect which triggers a database query. Use this function only for
 * debugging because many successive calls to this will cause a massive slow
 * down.
 *
 * \param[in] species   atomic species of the atom
 * \param[in] n         principal quantum number
 * \param[in] l         angular quantum number
 * \param[in] j         magnetic quantum number
 */
real_t energy_level(std::string const& species, int n, int l, real_t j);

#endif // QUANTUM_DEFECT_H
