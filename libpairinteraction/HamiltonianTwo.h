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

#ifndef HAMILTONIAN_TWO_H
#define HAMILTONIAN_TWO_H

#include "dtypes.h"
#include "ConfParser.h"
#include "MatrixElements.h"
#include "SQLite.h"
#include "Basisnames.h"
#include "Hamiltonian.h"
#include "HamiltonianOne.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/hex.hpp>
#include <boost/math/special_functions/binomial.hpp>

class HamiltonianTwo : public Hamiltonian<BasisnamesTwo> {
public:
    HamiltonianTwo(const Configuration &config, boost::filesystem::path& path_cache, const std::shared_ptr<HamiltonianOne>& hamiltonian_one);
    HamiltonianTwo(const Configuration &config, boost::filesystem::path& path_cache, std::shared_ptr<HamiltonianOne> hamiltonian_one1, std::shared_ptr<HamiltonianOne> hamiltonian_one2);
    void calculate(const Configuration &conf_tot);

private:
    std::shared_ptr<HamiltonianOne> hamiltonian_one1; // TODO const HamiltonianOne
    std::shared_ptr<HamiltonianOne> hamiltonian_one2;
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
    boost::filesystem::path path_cache;
};


#endif // HAMILTONIAN_TWO_H
