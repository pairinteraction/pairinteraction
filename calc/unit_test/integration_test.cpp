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
    system_one.setBfield({0,0,1});
    system_one.setEfield({0,0,0.1});

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
    system_two.diagonalize();
    eigen_vector_double_t eigenvalues = system_two.getDiagonal();
    eigen_sparse_t eigenvectors = system_two.getCoefficients();

    /*std::ofstream ofs("../calc/unit_test/integration_test_referencedata.txt");
    boost::archive::text_oarchive oa(ofs);
    oa << eigenvalues << eigenvectors;
    ofs.close();*/

    // Load reference data
    eigen_vector_double_t eigenvalues_reference;
    eigen_sparse_t eigenvectors_reference;
    std::ifstream ifs("./calc/unit_test/integration_test_referencedata.txt");
    boost::archive::text_iarchive ia(ifs);
    ia >> eigenvalues_reference >> eigenvectors_reference;

    // Compare current results to the reference data
    eigen_sparse_t diff = (eigenvectors - eigenvectors_reference).cwiseAbs();
    double max_diff_eigenvectors = *std::max_element(diff.valuePtr(),diff.valuePtr()+diff.nonZeros());
    std::cout << "Maximum deviation from reference eigenvectors: " << max_diff_eigenvectors << std::endl;
    BOOST_CHECK_CLOSE(max_diff_eigenvectors, 0, 1e-12);
    double max_diff_eigenvalues = (eigenvalues - eigenvalues_reference).cwiseAbs().maxCoeff();
    std::cout << "Maximum deviation from reference eigenvalues: " << max_diff_eigenvalues << std::endl;
    BOOST_CHECK_CLOSE(max_diff_eigenvalues, 0, 1e-12);

    // Delete cache directory
    boost::filesystem::remove_all(path_cache);

    // TODO call more methods to increase code covering
    // TODO cause exceptions and check whether they are handled correctly
}
