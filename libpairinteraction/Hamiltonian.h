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

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "Hamiltonianmatrix.h"

#include <vector>

template <class T>
class Hamiltonian {
public:
    Hamiltonian() {}
    std::shared_ptr<Hamiltonianmatrix> get(size_t idx) {
        return matrix_diag[idx];
    }
    std::shared_ptr<const Hamiltonianmatrix> get(size_t idx) const {
        return matrix_diag[idx];
    }
    std::shared_ptr<const Configuration> getParams(size_t idx) const {
        return params[idx];
    }
    size_t size() const {
        return matrix_diag.size();
    }
    std::shared_ptr<const T> names() const {
        return basis;
    }
    void removeUnnecessaryStates(std::vector<bool> &necessary) {
        basis->removeUnnecessaryStates(necessary);
        for (auto &p : matrix_diag) {
            p->removeUnnecessaryStates(necessary);
        }
    }

protected:
    std::vector<std::shared_ptr<Hamiltonianmatrix>> matrix_diag;
    std::vector<std::string> matrix_path;
    std::vector<std::shared_ptr<Configuration>> params;
    std::shared_ptr<T> basis;
};

#endif // HAMILTONIAN_H
