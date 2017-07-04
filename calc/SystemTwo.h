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

#ifndef SYSTEMTWO_H
#define SYSTEMTWO_H

#include "State.h"
#include "SystemBase.h"
#include "SystemOne.h"

#include <codecvt>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/unordered_map.hpp>


class SystemTwo : public SystemBase<StateTwo> {
public:
    SystemTwo(const SystemOne &b1, const SystemOne &b2, std::wstring cachedir);
    SystemTwo(const SystemOne &b1, const SystemOne &b2, std::wstring cachedir, bool memory_saving);
    SystemTwo(const SystemOne &b1, const SystemOne &b2);
    SystemTwo(const SystemOne &b1, const SystemOne &b2, bool memory_saving);

    const std::array<std::wstring, 2>& getElement();
    std::vector<StateOne> getStatesFirst();
    std::vector<StateOne> getStatesSecond();
    void setDistance(double d);
    void setAngle(double a);

    void setConservedParityUnderPermutation(parity_t parity);

protected:
    void initializeBasis() override;
    void initializeInteraction() override;
    void addInteraction() override;
    void transformInteraction(const eigen_sparse_t &transformator) override;
    void deleteInteraction() override;

private:
    void addCoefficient(const size_t &row_1, const size_t &row_2, const size_t &col_new, const scalar_t &value_new, std::vector<size_t> &stateidentifier2row, std::vector<eigen_triplet_t> &coefficients_triplets, std::vector<double> &sqnorm_list);

    std::array<std::wstring, 2> element;
    SystemOne system1; // is needed in the initializeBasis method and afterwards deleted
    SystemOne system2; // is needed in the initializeBasis method and afterwards deleted

    std::unordered_map<int, eigen_sparse_t> interaction;

    double distance;
    double angle;
    int kappa_max;

    parity_t sym_permutation;

    ////////////////////////////////////////////////////////////////////
    /// Method for serialization ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        (void)version;

        ar & boost::serialization::base_object<SystemBase<StateTwo>>(*this);
        ar & element & system1 & system2;
        ar & distance & angle & kappa_max & sym_permutation;
        ar & interaction;
    }
};

#endif
