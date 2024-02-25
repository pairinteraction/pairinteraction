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

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "Hamiltonianmatrix.hpp"

#include <vector>

template <typename Scalar, class Basisnames>
class Hamiltonian {
public:
    Hamiltonian() = default;
    std::shared_ptr<Hamiltonianmatrix<Scalar>> get(size_t idx) { return matrix_diag[idx]; }
    std::shared_ptr<const Hamiltonianmatrix<Scalar>> get(size_t idx) const {
        return matrix_diag[idx];
    }
    std::shared_ptr<const Configuration> getParams(size_t idx) const { return params[idx]; }
    size_t size() const { return matrix_diag.size(); }
    std::shared_ptr<const Basisnames> names() const { return basis; }
    void removeUnnecessaryStates(std::vector<bool> &necessary) {
        basis->removeUnnecessaryStates(necessary);
        for (auto &p : matrix_diag) {
            p->removeUnnecessaryStates(necessary);
        }
    }

protected:
    std::vector<std::shared_ptr<Hamiltonianmatrix<Scalar>>> matrix_diag;
    std::vector<std::string> matrix_path;
    std::vector<std::shared_ptr<Configuration>> params;
    std::shared_ptr<Basisnames> basis;
};

#endif // HAMILTONIAN_H
