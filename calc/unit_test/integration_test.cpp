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

#include "dtypes.h"
#include "State.h"
#include "SystemOne.h"
#include "SystemTwo.h"

#define BOOST_TEST_MODULE Integration test
#include <boost/test/unit_test.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <stdexcept>
#include <fstream>
#include <iostream>

BOOST_AUTO_TEST_CASE( integration_test ) {

    // Create cache directory
    boost::filesystem::path path_cache = boost::filesystem::temp_directory_path() / boost::filesystem::unique_path();
    if (boost::filesystem::create_directory(path_cache)) {
        std::cout << "Cache directory " << boost::filesystem::absolute(path_cache).string() << " created." << std::endl;
    } else {
        throw std::runtime_error( "Could not create cache directory." );
    }

    // Setup states
    StateOne state_one(L"Rb", 61, 2, 1.5, 1.5);
    StateTwo state_two(state_one, state_one);

    // Build one-atom system
    SystemOne system_one(state_one.element, path_cache.wstring());
    system_one.restrictEnergy(state_one.getEnergy()-40, state_one.getEnergy()+40);
    system_one.restrictN(state_one.n-1, state_one.n+1);
    system_one.restrictL(state_one.l-1, state_one.l+1);

    BOOST_CHECK_EQUAL(system_one.getNumVectors(), 64);
    BOOST_CHECK_EQUAL(system_one.getNumStates(), 64);

    // Build two-atom system
    SystemTwo system_two(system_one, system_one, path_cache.wstring());
    system_two.restrictEnergy(state_two.getEnergy()-2, state_two.getEnergy()+2);
    system_two.setConservedParityUnderPermutation(ODD);
    system_two.setDistance(6);
    system_two.setAngle(0.9);

    BOOST_CHECK_EQUAL(system_two.getNumVectors(), 239);
    BOOST_CHECK_EQUAL(system_two.getNumStates(), 468);

    // Diagonalize two-atom system
    eigen_sparse_t hamiltonian = system_two.getHamiltonianmatrix();
    eigen_sparse_t basis = system_two.getCoefficients();

    // Prune results (without pruning, max_diff_hamiltonian might be infinity due to division by zero)
    hamiltonian.prune(1e-12,1);
    basis.prune(1e-12,1);

    /*std::ofstream ofs("../calc/unit_test/integration_test_referencedata.txt");
    boost::archive::text_oarchive oa(ofs);
    oa << hamiltonian << basis;
    ofs.close();*/

    // Load reference data
    eigen_sparse_t hamiltonian_reference;
    eigen_sparse_t basis_reference;
    std::ifstream ifs("./calc/unit_test/integration_test_referencedata.txt");
    boost::archive::text_iarchive ia(ifs);
    ia >> hamiltonian_reference >> basis_reference;

    // Compare current results to the reference data (the results have to be compared before diagonalization as the order of the eigenvectors is not fixed)
    eigen_sparse_double_t diff = (hamiltonian - hamiltonian_reference).cwiseQuotient(hamiltonian.cwiseMin(hamiltonian_reference)).cwiseAbs();
    double max_diff_hamiltonian = *std::max_element(diff.valuePtr(),diff.valuePtr()+diff.nonZeros());
    std::cout << "Maximum deviation from reference Hamiltonian: " << max_diff_hamiltonian*100 << " %" << std::endl;
    BOOST_CHECK_SMALL(max_diff_hamiltonian, 1e-6);

    diff = (basis - basis_reference).cwiseQuotient(basis.cwiseMin(basis_reference)).cwiseAbs();
    double max_diff_basis = *std::max_element(diff.valuePtr(),diff.valuePtr()+diff.nonZeros());
    std::cout << "Maximum deviation from reference basis: " << max_diff_basis*100 << " %" << std::endl;
    BOOST_CHECK_SMALL(max_diff_basis, 1e-6);

    // Diagonalize two-atom system
    system_two.diagonalize();

    // Delete cache directory
    boost::filesystem::remove_all(path_cache);

    // TODO check Hamiltonian for system_one in case of electric and magnetic fields
    // TODO call more methods to increase code covering
    // TODO cause exceptions and check whether they are handled correctly
}
