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

#include "ConfParser.h"
#include "HamiltonianOne.h"
#include "HamiltonianTwo.h"

#include <iostream>
#include <memory>
#include <boost/filesystem.hpp>

/*
///////////////////// TODOs /////////////////////

* parallelize construction of Hamiltonian
* construct only one half of symmetric matrices
* check why very small Bz fields (e.g. 1e-12) leads to a large basis -> numerical error?

*/

int compute(const char* config_name, const char* output_name) {
    std::cout << std::unitbuf;

    Eigen::setNbThreads(1); // TODO set it to setNbThreads(0) when Eigen's multithreading is needed

    boost::filesystem::path path_config = boost::filesystem::absolute(config_name);
    boost::filesystem::path path_cache  = boost::filesystem::absolute(output_name);

    // === Load configuration ===
    Configuration config;
    config.load_from_json(path_config.string());

    bool existAtom1 = config.count("species1") && config.count("n1") && config.count("l1") && config.count("j1") && config.count("m1");
    bool existAtom2 = config.count("species2") && config.count("n2") && config.count("l2") && config.count("j2") && config.count("m2");

    // === Solve the system ===
    bool combined = config["samebasis"].str() == "true";

    if (combined) {
        if (config["species1"].str() != config["species2"].str()) {
            std::cout << "species1 and species2 has to be the same in order to use the same basis set." << std::endl;
            return 1;
        }
        std::shared_ptr<HamiltonianOne> hamiltonian_one;
        if (existAtom1 && existAtom2) {
            std::cout << ">>TYP" << std::setw(7) << 3 << std::endl;
            auto basisnames_one = std::make_shared<BasisnamesOne>(BasisnamesOne::fromBoth(config));
            hamiltonian_one = std::make_shared<HamiltonianOne>(config, path_cache, basisnames_one);
        }
        std::shared_ptr<HamiltonianTwo> hamiltonian_two;
        if (existAtom1 && existAtom2 && config.count("minR")) {
            std::cout << ">>TYP" << std::setw(7) << 2 << std::endl;
            hamiltonian_two = std::make_shared<HamiltonianTwo>(config, path_cache, hamiltonian_one);
        }
    } else {
        std::shared_ptr<HamiltonianOne> hamiltonian_one1;
        if (existAtom1) {
            std::cout << ">>TYP" << std::setw(7) << 0 << std::endl;
            auto basisnames_one1 = std::make_shared<BasisnamesOne>(BasisnamesOne::fromFirst(config));
            hamiltonian_one1 = std::make_shared<HamiltonianOne>(config, path_cache, basisnames_one1);
        }
        std::shared_ptr<HamiltonianOne> hamiltonian_one2;
        if (existAtom2) {
            std::cout << ">>TYP" << std::setw(7) << 1 << std::endl;
            auto basisnames_one2 = std::make_shared<BasisnamesOne>(BasisnamesOne::fromSecond(config));
            hamiltonian_one2 = std::make_shared<HamiltonianOne>(config, path_cache, basisnames_one2);
        }
        std::shared_ptr<HamiltonianTwo> hamiltonian_two;
        if (existAtom1 && existAtom2 && config.count("minR")) {
            std::cout << ">>TYP" << std::setw(7) << 2 << std::endl;
            hamiltonian_two = std::make_shared<HamiltonianTwo>(config, path_cache, hamiltonian_one1, hamiltonian_one2);
        }
    }

    // === Communicate that everything has finished ===
    std::cout << ">>END" << std::endl;

    return 0;
}
