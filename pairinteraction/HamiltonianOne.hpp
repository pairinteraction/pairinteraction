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

#ifndef HAMILTONIAN_ONE_H
#define HAMILTONIAN_ONE_H

#include "Basisnames.hpp"
#include "ConfParser.hpp"
#include "Hamiltonian.hpp"
#include "MatrixElements.hpp"
#include "SQLite.hpp"
#include "filesystem.hpp"

#include <cmath>
#include <iostream>
#include <memory>

template <typename Scalar>
class HamiltonianOne : public Hamiltonian<Scalar, BasisnamesOne> {
public:
    HamiltonianOne(const Configuration &config, fs::path &path_cache,
                   std::shared_ptr<BasisnamesOne> basis_one);
    const Configuration &getConf() const;

protected:
    void changeToSpherical(double val_x, double val_y, double val_z, double &val_p, double &val_m,
                           double &val_0);
    void changeToSpherical(double val_x, double val_y, double val_z, std::complex<double> &val_p,
                           std::complex<double> &val_m, std::complex<double> &val_0);
    void configure(const Configuration &config);
    void build();

private:
    Configuration basicconf;
    double deltaE;
    double min_E_x, min_E_y, min_E_z, max_E_x, max_E_y, max_E_z, min_B_x, min_B_y, min_B_z, max_B_x,
        max_B_y, max_B_z;
    size_t nSteps;
    bool diamagnetism;
    std::string species;
    fs::path path_cache;
};

extern template class HamiltonianOne<std::complex<double>>;
extern template class HamiltonianOne<double>;

#endif // HAMILTONIAN_ONE_H
