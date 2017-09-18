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

#include <iostream>

BOOST_AUTO_TEST_CASE( integration_test ) {
    StateOne state_one(L"Rb", 61, 2, 1.5, 1.5);
    StateTwo state_two(state_one, state_one);

    // TODO create database directory

    SystemOne system_one(state_one.element);
    system_one.restrictEnergy(state_one.getEnergy()-40, state_one.getEnergy()+40);
    system_one.restrictN(state_one.n-2, state_one.n+2);
    system_one.restrictL(state_one.l-2, state_one.l+2);
    system_one.setBfield({0,0,1});
    system_one.setEfield({0,0,0.1});

    SystemTwo system_two(system_one, system_one);
    system_two.restrictEnergy(state_two.getEnergy()-2, state_two.getEnergy()+2);
    system_two.setConservedParityUnderPermutation(ODD);
    system_two.setDistance(6);
    system_two.setAngle(0.9);

    system_two.diagonalize();

    // TODO delete database directory

    // TODO use BOOST_CHECK_EQUAL to verify the result
}
