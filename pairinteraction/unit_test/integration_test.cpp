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

#include "MatrixElementCache.h"
#include "State.h"
#include "SystemOne.h"
#include "SystemTwo.h"
#include "dtypes.h"
#include "filesystem.h"

#define BOOST_TEST_MODULE Integration test
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

struct F {
    F() : path_cache(fs::create_temp_directory()) {
        std::cout << "Cache directory " << fs::absolute(path_cache).string() << " created."
                  << std::endl;
    }
    ~F() {
        // Delete cache directory
        fs::remove_all(path_cache);
    }
    fs::path path_cache;
};

BOOST_FIXTURE_TEST_CASE(integration_test, F) // NOLINT
{

    constexpr bool dump_new_reference_data = false;
    constexpr double tolerance = 1e-6;

    ////////////////////////////////////////////////////////////////////
    /// Preparations ///////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // Set up cache
    MatrixElementCache cache(path_cache.string());

    // Load reference data
    eigen_sparse_t hamiltonian_one_reference, hamiltonian_two_reference;
    eigen_sparse_t basis_one_reference, basis_two_reference;

    if (!dump_new_reference_data) {
        std::ifstream ifs("./pairinteraction/unit_test/integration_test_referencedata.txt");
        boost::archive::text_iarchive ia(ifs);
        ia >> hamiltonian_one_reference >> basis_one_reference >> hamiltonian_two_reference >>
            basis_two_reference;
        ifs.close();
    }

    // Setup states
    StateOne state_one("Rb", 61, 2, 1.5, 1.5);
    StateTwo state_two(state_one, state_one);

    ////////////////////////////////////////////////////////////////////
    /// Test system consisting of one atom /////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // Build one-atom system
    SystemOne system_one(state_one.getSpecies(), cache);
    system_one.restrictEnergy(state_one.getEnergy() - 40, state_one.getEnergy() + 40);
    system_one.restrictN(state_one.getN() - 1, state_one.getN() + 1);
    system_one.restrictL(state_one.getL() - 1, state_one.getL() + 1);
    system_one.setEfield({{0, 0, 0.1}});
    system_one.setBfield({{0, 0, 1}});
    system_one.enableDiamagnetism(false);

    // Check for correct dimensions
    BOOST_CHECK_EQUAL(system_one.getNumBasisvectors(), 64);
    BOOST_CHECK_EQUAL(system_one.getNumStates(), 64);

    // Compare current results to the reference data (the results have to be
    // compared before diagonalization as the order of the eigenvectors is not
    // fixed)
    eigen_sparse_t hamiltonian_one = system_one.getHamiltonian();
    eigen_sparse_t basis_one = system_one.getBasisvectors();

    eigen_sparse_double_t diff;
    double max_diff_hamiltonian, max_diff_basis;

    if (!dump_new_reference_data) {
        diff = (hamiltonian_one - hamiltonian_one_reference)
                   .pruned(tolerance, 1) // without pruning, max_diff_hamiltonian
                                         // might be infinity due to division by
                                         // zero
                   .cwiseQuotient(hamiltonian_one.cwiseMin(hamiltonian_one_reference))
                   .cwiseAbs();
        max_diff_hamiltonian =
            *std::max_element(diff.valuePtr(), diff.valuePtr() + diff.nonZeros());
        std::cout << "One-atom system, relative maximum deviation from "
                     "reference Hamiltonian: "
                  << max_diff_hamiltonian << std::endl;
        BOOST_CHECK_SMALL(max_diff_hamiltonian, tolerance);

        diff = (basis_one - basis_one_reference)
                   .pruned(tolerance, 1) // without pruning, max_diff_hamiltonian
                                         // might be infinity due to division by
                                         // zero
                   .cwiseQuotient(basis_one.cwiseMin(basis_one_reference))
                   .cwiseAbs();
        max_diff_basis = *std::max_element(diff.valuePtr(), diff.valuePtr() + diff.nonZeros());
        std::cout << "One-atom system, relative maximum deviation from "
                     "reference basis: "
                  << max_diff_basis << std::endl;
        BOOST_CHECK_SMALL(max_diff_basis, tolerance);
    }

    // Diagonalize one-atom system
    system_one.diagonalize();

    ////////////////////////////////////////////////////////////////////
    /// Test system consisting of two atoms ////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // Build one-atom system (for this test, system_one has to be diagonal by
    // itself because diagonalization can lead to different order of
    // eigenvectors)
    // system_one = SystemOne(state_one.species, cache); // TODO  object of type 'SystemOne' cannot
    // be assigned because its copy assignment operator is implicitly deleted
    SystemOne system_one_new(state_one.getSpecies(), cache);
    system_one_new.restrictEnergy(state_one.getEnergy() - 40, state_one.getEnergy() + 40);
    system_one_new.restrictN(state_one.getN() - 1, state_one.getN() + 1);
    system_one_new.restrictL(state_one.getL() - 1, state_one.getL() + 1);

    // Build two-atom system
    SystemTwo system_two(system_one_new, system_one_new, cache);
    system_two.restrictEnergy(state_two.getEnergy() - 2, state_two.getEnergy() + 2);
    system_two.setConservedParityUnderPermutation(ODD);
    system_two.setDistance(6);
    system_two.setAngle(0.9);

    // Check for correct dimensions
    BOOST_CHECK_EQUAL(system_two.getNumBasisvectors(), 239);
    BOOST_CHECK_EQUAL(system_two.getNumStates(), 468);

    // Compare current results to the reference data (the results have to be
    // compared before diagonalization as the order of the eigenvectors is not
    // fixed)
    eigen_sparse_t hamiltonian_two = system_two.getHamiltonian();
    eigen_sparse_t basis_two = system_two.getBasisvectors();

    if (!dump_new_reference_data) {
        diff = (hamiltonian_two - hamiltonian_two_reference)
                   .pruned(tolerance, 1) // without pruning, max_diff_hamiltonian
                                         // might be infinity due to division by
                                         // zero
                   .cwiseQuotient(hamiltonian_two.cwiseMin(hamiltonian_two_reference))
                   .cwiseAbs();
        max_diff_hamiltonian =
            *std::max_element(diff.valuePtr(), diff.valuePtr() + diff.nonZeros());
        std::cout << "Two-atom system, relative maximum deviation from "
                     "reference Hamiltonian: "
                  << max_diff_hamiltonian << std::endl;
        BOOST_CHECK_SMALL(max_diff_hamiltonian, tolerance);

        diff = (basis_two - basis_two_reference)
                   .pruned(tolerance, 1) // without pruning, max_diff_hamiltonian
                                         // might be infinity due to division by
                                         // zero
                   .cwiseQuotient(basis_two.cwiseMin(basis_two_reference))
                   .cwiseAbs();
        max_diff_basis = *std::max_element(diff.valuePtr(), diff.valuePtr() + diff.nonZeros());
        std::cout << "Two-atom system, relative maximum deviation from "
                     "reference basis: "
                  << max_diff_basis << std::endl;
        BOOST_CHECK_SMALL(max_diff_basis, tolerance);
    }

    // Diagonalize two-atom system
    system_two.diagonalize();

    ////////////////////////////////////////////////////////////////////
    /// Clean up ///////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    if (dump_new_reference_data) {
        std::ofstream ofs("../pairinteraction/unit_test/integration_test_referencedata.txt");
        boost::archive::text_oarchive oa(ofs);
        oa << hamiltonian_one << basis_one << hamiltonian_two << basis_two;
        ofs.close();

        // ATTENTION
        // After generating integration_test_referencedata.txt, we possibly have to manually modify
        // the Boost serialization library version - the beginning of the file should read "22
        // serialization::archive 12". Otherwise, unsupported version exceptions are thrown with
        // older versions of Boost even though they are compatible.

        // We deliberately crash the test in case we are dumping new
        // data because we want to prevent accidentally dumping new
        // data when running in Travis CI.
        BOOST_CHECK_MESSAGE(!dump_new_reference_data, "No tests were executed. Only dumping data!");
    }

    // TODO call more methods to increase code covering
    // TODO cause exceptions and check whether they are handled correctly
}
