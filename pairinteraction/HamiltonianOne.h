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

#ifndef HAMILTONIAN_ONE_H
#define HAMILTONIAN_ONE_H

#include "dtypes.h"
#include "ConfParser.h"
#include "MatrixElements.h"
#include "SQLite.h"
#include "Basisnames.h"
#include "Hamiltonian.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/hex.hpp>

class HamiltonianOne : public Hamiltonian<BasisnamesOne>{
public:
    HamiltonianOne(const Configuration &config, boost::filesystem::path& path_cache, std::shared_ptr<BasisnamesOne> basis_one);
    const Configuration& getConf() const;

protected:
    void changeToSpherical(double val_x, double val_y, double val_z, double& val_p, double& val_m, double& val_0);
    void changeToSpherical(double val_x, double val_y, double val_z, std::complex<double>& val_p, std::complex<double>& val_m, std::complex<double>& val_0);
    void configure(const Configuration &config);
    void build();

private:
    Configuration basicconf;
    double deltaE;
    double min_E_x,min_E_y,min_E_z,max_E_x,max_E_y,max_E_z,min_B_x,min_B_y,min_B_z,max_B_x,max_B_y,max_B_z;
    size_t nSteps;
    bool diamagnetism;
    std::string species;
    boost::filesystem::path path_cache;
};

#endif // HAMILTONIAN_ONE_H
