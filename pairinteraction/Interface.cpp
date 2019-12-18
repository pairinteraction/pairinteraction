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

#include "Interface.h"
#include "ConfParser.h"
#include "HamiltonianOne.h"
#include "HamiltonianTwo.h"
#include "filesystem.h"

#include <fmt/format.h>

#include <iostream>
#include <locale>
#include <memory>

/*
///////////////////// TODOs /////////////////////

* parallelize construction of Hamiltonian
* construct only one half of symmetric matrices
* check why very small Bz fields (e.g. 1e-12) leads to a large basis -> numerical error?

*/

int compute(const std::string &config_name, const std::string &output_name) {
    std::cout << std::unitbuf;

    Eigen::setNbThreads(1); // TODO set it to setNbThreads(0) when Eigen's multithreading is needed

    fs::path path_config = fs::absolute(config_name);
    fs::path path_cache = fs::absolute(output_name);

    // === Load configuration ===
    Configuration config;
    config.load_from_json(path_config.string());

    bool existAtom1 = (config.count("species1") != 0u) && (config.count("n1") != 0u) &&
        (config.count("l1") != 0u) && (config.count("j1") != 0u) && (config.count("m1") != 0u);
    bool existAtom2 = (config.count("species2") != 0u) && (config.count("n2") != 0u) &&
        (config.count("l2") != 0u) && (config.count("j2") != 0u) && (config.count("m2") != 0u);

    // === Solve the system ===
    bool combined = config["samebasis"].str() == "true";

    if (combined) {
        if (config["species1"].str() != config["species2"].str()) {
            std::cout
                << "species1 and species2 has to be the same in order to use the same basis set."
                << std::endl;
            return 1;
        }
        std::shared_ptr<HamiltonianOne> hamiltonian_one;
        if (existAtom1 && existAtom2) {
            std::cout << fmt::format(">>TYP{:7d}", 3) << std::endl;
            auto basisnames_one = std::make_shared<BasisnamesOne>(BasisnamesOne::fromBoth(config));
            hamiltonian_one = std::make_shared<HamiltonianOne>(config, path_cache, basisnames_one);
        }
        std::shared_ptr<HamiltonianTwo> hamiltonian_two;
        if (existAtom1 && existAtom2 && (config.count("minR") != 0u)) {
            std::cout << fmt::format(">>TYP{:7d}", 2) << std::endl;
            hamiltonian_two = std::make_shared<HamiltonianTwo>(config, path_cache, hamiltonian_one);
        }
    } else {
        std::shared_ptr<HamiltonianOne> hamiltonian_one1;
        if (existAtom1) {
            std::cout << fmt::format(">>TYP{:7d}", 0) << std::endl;
            auto basisnames_one1 =
                std::make_shared<BasisnamesOne>(BasisnamesOne::fromFirst(config));
            hamiltonian_one1 =
                std::make_shared<HamiltonianOne>(config, path_cache, basisnames_one1);
        }
        std::shared_ptr<HamiltonianOne> hamiltonian_one2;
        if (existAtom2) {
            std::cout << fmt::format(">>TYP{:7d}", 1) << std::endl;
            auto basisnames_one2 =
                std::make_shared<BasisnamesOne>(BasisnamesOne::fromSecond(config));
            hamiltonian_one2 =
                std::make_shared<HamiltonianOne>(config, path_cache, basisnames_one2);
        }
        std::shared_ptr<HamiltonianTwo> hamiltonian_two;
        if (existAtom1 && existAtom2 && (config.count("minR") != 0u)) {
            std::cout << fmt::format(">>TYP{:7d}", 2) << std::endl;
            hamiltonian_two = std::make_shared<HamiltonianTwo>(config, path_cache, hamiltonian_one1,
                                                               hamiltonian_one2);
        }
    }

    // === Communicate that everything has finished ===
    std::cout << ">>END" << std::endl;

    return 0;
}
