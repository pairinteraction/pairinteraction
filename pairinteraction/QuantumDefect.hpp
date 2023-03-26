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

#ifndef QUANTUM_DEFECT_H
#define QUANTUM_DEFECT_H

#include "Cache.hpp"
#include "SQLite.hpp"
#include "dtypes.hpp"
#include "utils.hpp"

#include <mutex>
#include <string>
#include <tuple>
#include <unordered_map>

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
    /** \brief Key type for quantum defect cache */
    struct Key {
        std::string species;
        int n;
        int l;
        double j;

        /** \brief Comparison for keys */
        bool operator==(Key const &o) const {
            return (species == o.species && n == o.n && l == o.l && j == o.j);
        }
    };

    //** \brief Hash for Key */
    struct Hash {
        /** \brief Hash the key using utils::hash_combine */
        std::size_t operator()(Key const &key) const {
            std::size_t seed = 0;
            utils::hash_combine(seed, key.species);
            utils::hash_combine(seed, key.n);
            utils::hash_combine(seed, key.l);
            utils::hash_combine(seed, key.j);
            return seed;
        }
    };

    /** \brief Element in the quantum defect cache */
    struct Element {
        double ac;
        int Z;
        double a1, a2, a3, a4;
        double rc;
        double nstar;
        double energy;
    };

    Element e;

    void setup(sqlite3 * /*db*/, std::string const &db_name);

    QuantumDefect(std::string _species, int _n, int _l, double _j, std::nullptr_t);

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
    QuantumDefect(std::string const &species, int n, int l, double j);

    /** \overload QuantumDefect(std::string const& species, int n, int l, double
     * j) */
    QuantumDefect(std::string const &species, int n, int l, double j, std::string const &database);

    /** \brief Atomic species */
    const std::string species;

    /** \brief Principal quantum number */
    const int n;

    /** \brief Angular quantum number */
    const int l;

    /** \brief Magnetic quantum number */
    const double j;

    /** \brief Polarizability */
    const double &ac;

    /** \brief Core charge */
    const int &Z;

    /** \brief Parameter of effective Coulomb potential */
    const double &a1;

    /** \brief Parameter of effective Coulomb potential */
    const double &a2;

    /** \brief Parameter of effective Coulomb potential */
    const double &a3;

    /** \brief Parameter of effective Coulomb potential */
    const double &a4;

    /** \brief Effective core size */
    const double &rc;

    /** \brief Effective principal quantum number */
    const double &nstar;

    /** \brief State energy */
    const double &energy;
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
double energy_level(std::string const &species, int n, int l, double j,
                    std::string const &database = "");

double nstar(std::string const &species, int n, int l, double j, std::string const &database = "");

#endif // QUANTUM_DEFECT_H
