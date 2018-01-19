/*
 * Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
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

#include <boost/filesystem.hpp>
#include <stdexcept>
#include <iostream>

#include "libpairinteraction/dtypes.h"
#include "libpairinteraction/State.h"
#include "libpairinteraction/SystemOne.h"
#include "libpairinteraction/SystemTwo.h"

int main() {
    // Create cache directory
    boost::filesystem::path path_cache =
        boost::filesystem::temp_directory_path()/boost::filesystem::unique_path();
    if (boost::filesystem::create_directory(path_cache)) {
        std::cout << "Cache directory "
                    << boost::filesystem::absolute(path_cache).string()
                    << " created." << std::endl;
    } else {
        throw std::runtime_error("Could not create cache directory.");
    }

    // Setup states
    StateOne state_one("Rb", 61, 2, 1.5, 1.5);

    // Build one-atom system
    SystemOne system_one(state_one.element, path_cache.string());
    system_one.restrictEnergy(state_one.getEnergy() - 20, state_one.getEnergy() + 20);
    system_one.restrictN(state_one.n - 1, state_one.n + 1);
    system_one.restrictL(state_one.l - 1, state_one.l + 1);
    system_one.setEfield({{0, 0, 0.1}});
    system_one.setBfield({{0, 0, 1}});
    
    // Print Hamiltonian
    std::cout << system_one.getHamiltonianmatrix() << std::endl;
    
    // Diagonalize one-atom system
    system_one.diagonalize();

    // Print Hamiltonian
    std::cout << system_one.getHamiltonianmatrix() << std::endl;
    
    // Remove cache directory
    if (boost::filesystem::remove_all(path_cache)) {
        std::cout << "Cache directory "
                    << boost::filesystem::absolute(path_cache).string()
                    << " removed." << std::endl;
    } else {
        throw std::runtime_error("Could not remove cache directory.");
    }
}
