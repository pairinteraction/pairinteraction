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

#include "State.h"
#include "SystemOne.h"
#include "SystemTwo.h"
#include "dtypes.h"

#define BOOST_TEST_MODULE Integration test
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>

BOOST_AUTO_TEST_CASE(integration_test)
{

    constexpr bool dump_new_reference_data = false;

    ////////////////////////////////////////////////////////////////////
    /// Preparations ///////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // Create cache directory
    boost::filesystem::path path_cache =
        boost::filesystem::temp_directory_path() /
        boost::filesystem::unique_path();
    if (boost::filesystem::create_directory(path_cache)) {
        std::cout << "Cache directory "
                  << boost::filesystem::absolute(path_cache).string()
                  << " created." << std::endl;
    } else {
        throw std::runtime_error("Could not create cache directory.");
    }

    // Load reference data
    eigen_sparse_t hamiltonian_one_reference, hamiltonian_two_reference;
    eigen_sparse_t basis_one_reference, basis_two_reference;

    if (!dump_new_reference_data) {
        std::ifstream ifs(
            "./libpairinteraction/unit_test/integration_test_referencedata.txt");
        boost::archive::text_iarchive ia(ifs);
        ia >> hamiltonian_one_reference >> basis_one_reference >>
            hamiltonian_two_reference >> basis_two_reference;
        ifs.close();
    }

    // Setup states
    StateOne state_one("Rb", 61, 2, 1.5, 1.5);
    StateTwo state_two(state_one, state_one);

    ////////////////////////////////////////////////////////////////////
    /// Test system consisting of one atom /////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // Build one-atom system
    SystemOne system_one(state_one.element, path_cache.string());
    system_one.restrictEnergy(state_one.getEnergy() - 40,
                              state_one.getEnergy() + 40);
    system_one.restrictN(state_one.n - 1, state_one.n + 1);
    system_one.restrictL(state_one.l - 1, state_one.l + 1);
    system_one.setEfield({0, 0, 0.1});
    system_one.setBfield({0, 0, 1});

    // Check for correct dimensions
    BOOST_CHECK_EQUAL(system_one.getNumVectors(), 64);
    BOOST_CHECK_EQUAL(system_one.getNumStates(), 64);

    // Compare current results to the reference data (the results have to be
    // compared before diagonalization as the order of the eigenvectors is not
    // fixed)
    eigen_sparse_t hamiltonian_one = system_one.getHamiltonianmatrix();
    eigen_sparse_t basis_one = system_one.getCoefficients();
    hamiltonian_one.prune(1e-16, 1); // without pruning, max_diff_hamiltonian
                                     // might be infinity due to division by
                                     // zero
    basis_one.prune(1e-16, 1);       // without pruning, max_diff_basis might be
                                     // infinity due to division by zero

    eigen_sparse_double_t diff;
    double max_diff_hamiltonian, max_diff_basis;

    if (!dump_new_reference_data) {
        diff = (hamiltonian_one - hamiltonian_one_reference)
                   .cwiseQuotient(
                       hamiltonian_one.cwiseMin(hamiltonian_one_reference))
                   .cwiseAbs();
        max_diff_hamiltonian = *std::max_element(
            diff.valuePtr(), diff.valuePtr() + diff.nonZeros());
        std::cout << "One-atom system, relative maximum deviation from "
                     "reference Hamiltonian: "
                  << max_diff_hamiltonian << std::endl;
        BOOST_CHECK_SMALL(max_diff_hamiltonian, 1e-6);

        diff = (basis_one - basis_one_reference)
                   .cwiseQuotient(basis_one.cwiseMin(basis_one_reference))
                   .cwiseAbs();
        max_diff_basis = *std::max_element(diff.valuePtr(),
                                           diff.valuePtr() + diff.nonZeros());
        std::cout << "One-atom system, relative maximum deviation from "
                     "reference basis: "
                  << max_diff_basis << std::endl;
        BOOST_CHECK_SMALL(max_diff_basis, 1e-9);
    }

    // Diagonalize one-atom system
    system_one.diagonalize();

    ////////////////////////////////////////////////////////////////////
    /// Test system consisting of two atoms ////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // Build one-atom system (for this test, system_one has to be diagonal by
    // itself because diagonalization can lead to different order of
    // eigenvectors)
    system_one = SystemOne(state_one.element, path_cache.string());
    system_one.restrictEnergy(state_one.getEnergy() - 40,
                              state_one.getEnergy() + 40);
    system_one.restrictN(state_one.n - 1, state_one.n + 1);
    system_one.restrictL(state_one.l - 1, state_one.l + 1);

    // Build two-atom system
    SystemTwo system_two(system_one, system_one, path_cache.string());
    system_two.restrictEnergy(state_two.getEnergy() - 2,
                              state_two.getEnergy() + 2);
    system_two.setConservedParityUnderPermutation(ODD);
    system_two.setDistance(6);
    system_two.setAngle(0.9);

    // Check for correct dimensions
    BOOST_CHECK_EQUAL(system_two.getNumVectors(), 239);
    BOOST_CHECK_EQUAL(system_two.getNumStates(), 468);

    // Compare current results to the reference data (the results have to be
    // compared before diagonalization as the order of the eigenvectors is not
    // fixed)
    eigen_sparse_t hamiltonian_two = system_two.getHamiltonianmatrix();
    eigen_sparse_t basis_two = system_two.getCoefficients();
    hamiltonian_two.prune(1e-16, 1); // without pruning, max_diff_hamiltonian
                                     // might be infinity due to division by
                                     // zero
    basis_two.prune(1e-16, 1);       // without pruning, max_diff_basis might be
                                     // infinity due to division by zero

    if (!dump_new_reference_data) {
        diff = (hamiltonian_two - hamiltonian_two_reference)
                   .cwiseQuotient(
                       hamiltonian_two.cwiseMin(hamiltonian_two_reference))
                   .cwiseAbs();
        max_diff_hamiltonian = *std::max_element(
            diff.valuePtr(), diff.valuePtr() + diff.nonZeros());
        std::cout << "Two-atom system, relative maximum deviation from "
                     "reference Hamiltonian: "
                  << max_diff_hamiltonian << std::endl;
        BOOST_CHECK_SMALL(max_diff_hamiltonian, 1e-6);

        diff = (basis_two - basis_two_reference)
                   .cwiseQuotient(basis_two.cwiseMin(basis_two_reference))
                   .cwiseAbs();
        max_diff_basis = *std::max_element(diff.valuePtr(),
                                           diff.valuePtr() + diff.nonZeros());
        std::cout << "Two-atom system, relative maximum deviation from "
                     "reference basis: "
                  << max_diff_basis << std::endl;
        BOOST_CHECK_SMALL(max_diff_basis, 1e-9);
    }

    // Diagonalize two-atom system
    system_two.diagonalize();

    ////////////////////////////////////////////////////////////////////
    /// Clean up ///////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    if (dump_new_reference_data) {
        std::ofstream ofs(
            "../libpairinteraction/unit_test/integration_test_referencedata.txt");
        boost::archive::text_oarchive oa(ofs);
        oa << hamiltonian_one << basis_one << hamiltonian_two << basis_two;
        ofs.close();

        // We deliberately crash the test in case we are dumping new
        // data because we want to prevent accidentally dumping new
        // data when running in Travis CI.
        BOOST_CHECK_MESSAGE(!dump_new_reference_data,
                            "No tests were executed. Only dumping data!");
    }

    // Delete cache directory
    boost::filesystem::remove_all(path_cache);

    // TODO call more methods to increase code covering
    // TODO cause exceptions and check whether they are handled correctly
}
