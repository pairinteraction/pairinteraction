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

#ifndef HAMILTONIAN_TWO_H
#define HAMILTONIAN_TWO_H

#include "Basisnames.hpp"
#include "ConfParser.hpp"
#include "Hamiltonian.hpp"
#include "HamiltonianOne.hpp"
#include "MatrixElements.hpp"
#include "SQLite.hpp"
#include "dtypes.hpp"
#include "filesystem.hpp"

#include <boost/algorithm/hex.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include <cmath>
#include <iostream>
#include <memory>

template <typename Scalar>
class HamiltonianTwo : public Hamiltonian<Scalar,BasisnamesTwo> {
public:
    HamiltonianTwo(const Configuration &config, fs::path &path_cache,
                   const std::shared_ptr<HamiltonianOne<Scalar>> &hamiltonian_one);
    HamiltonianTwo(const Configuration &config, fs::path &path_cache,
                   std::shared_ptr<HamiltonianOne<Scalar>> hamiltonian_one1,
                   std::shared_ptr<HamiltonianOne<Scalar>> hamiltonian_one2);
    void calculate(const Configuration &conf_tot);

private:
    std::shared_ptr<HamiltonianOne<Scalar>> hamiltonian_one1; // TODO const HamiltonianOne
    std::shared_ptr<HamiltonianOne<Scalar>> hamiltonian_one2;
    double deltaE;
    int deltaN;
    int deltaL;
    int deltaJ;
    int deltaM;
    size_t nSteps_two;
    std::string species1, species2;
    double min_R, max_R;
    int multipoleexponent;
    bool samebasis;
    fs::path path_cache;
};

#ifdef USE_COMPLEX
extern template class HamiltonianTwo<std::complex<double>>;
#else
extern template class HamiltonianTwo<double>;
#endif

#endif // HAMILTONIAN_TWO_H
