/*
 * Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
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

#include <boost/filesystem.hpp>
#include <iostream>
#include <pairinteraction/State>
#include <pairinteraction/System>

int main() {
    // Create cache directory if not already existing
    boost::filesystem::path path_cache = boost::filesystem::current_path() / "cache";
    if (boost::filesystem::create_directory(path_cache)) {
        std::cout << "Cache directory " << path_cache.string() << " created."
                  << "\n";
    }

    // Set up cache
    MatrixElementCache cache(path_cache.string());

    // Set up states
    StateOne state_one("Rb", 61, 2, 1.5, 1.5);

    // Build one-atom system
    SystemOne system_one(state_one.getSpecies(), cache);
    system_one.restrictEnergy(state_one.getEnergy() - 20, state_one.getEnergy() + 20);
    system_one.restrictN(state_one.getN() - 1, state_one.getN() + 1);
    system_one.restrictL(state_one.getL() - 1, state_one.getL() + 1);
    system_one.setEfield({{0, 0, 0.1}});
    system_one.setBfield({{0, 0, 1}});

    // Print Hamiltonian
    std::cout << system_one.getHamiltonian() << "\n";

    // Diagonalize one-atom system
    system_one.diagonalize();

    // Print Hamiltonian
    std::cout << system_one.getHamiltonian() << "\n";
}
