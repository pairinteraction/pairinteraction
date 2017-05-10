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

#include "Basis.h"


// Implementation of BasisOne

BasisOne::BasisOne(std::string const& element)
    : element(element)
{}

void BasisOne::initialize()
{
    // TODO check whether specified basis is finite

    size_t idx = 0;
    std::vector<eigen_triplet_t> coefficients_triplets; // TODO reserve states, energies, basis_triplets
    std::vector<int> range_adapted_n, range_adapted_l;
    std::vector<double> range_adapted_j, range_adapted_m;

    if (range_n.empty()) {
        range_adapted_n = std::vector<int>({}); // TODO if empty, calculate the range via the energies
    } else {
        range_adapted_n = range_n;
    }
    for (auto n : range_adapted_n) {

        if (range_l.empty()) {
            range_adapted_l.resize(n);
            std::iota(range_adapted_l.begin(), range_adapted_l.end(), 0);
        } else {
            range_adapted_l = range_l;
        }
        for (auto l : range_adapted_l) {
            if (l > n-1) continue;

            if (range_j.empty()) {
                range_adapted_j = (l == 0) ? std::vector<double>({l+0.5f}) : std::vector<double>({std::fabs(l-0.5f), l+0.5f});
            } else {
                range_adapted_j = range_j;
            }
            for (auto j : range_adapted_j) {
                if (std::fabs(j-l) != 0.5) continue;

                real_t energy = StateOne(element,n,l,j,0.5).getEnergy();
                if ((energy < energy_min && energy_min != std::numeric_limits<real_t>::lowest()) || (energy > energy_max  && energy_max != std::numeric_limits<real_t>::max())) continue; // TODO take into account numerical errors

                if (range_m.empty()) {
                    range_adapted_m.resize(2*j+1);
                    std::iota(range_adapted_m.begin(), range_adapted_m.end(), -j);
                } else {
                    range_adapted_m = range_m;
                }
                for (auto m : range_adapted_m) {
                    if (std::fabs(m) > j) continue;

                    states.push_back(StateOne(element,n,l,j,m));
                    energies.push_back(energy);
                    coefficients_triplets.push_back(eigen_triplet_t(idx,idx,1)); // TODO take into account symmetries

                    ++idx;
                }
            }
        }
    }

    coefficients.resize(idx,idx);
    coefficients.setFromTriplets(coefficients_triplets.begin(), coefficients_triplets.end());
    coefficients_triplets.clear();
}


// Implementation of BasisOne

BasisTwo::BasisTwo(const BasisOne &b1, const BasisOne &b2)
    : basis1(b1), basis2(b2)
{}

BasisOne BasisTwo::getFirstBasis() const { return basis1; }
void BasisTwo::setFirstBasis(const BasisOne &b) { basis1 = b; }

BasisOne BasisTwo::getSecondBasis() const { return basis2; }
void BasisTwo::setSecondBasis(const BasisOne &b) { basis2 = b; }

void BasisTwo::initialize()
{
    // Restrict one atom states to the allowed quantum numbers
    basis1.restrictN(range_n);
    basis1.restrictL(range_l);
    basis1.restrictJ(range_j);
    basis1.restrictM(range_m);
    basis2.restrictN(range_n);
    basis2.restrictL(range_l);
    basis2.restrictJ(range_j);
    basis2.restrictM(range_m);

    // Combine the one atom states
    // TODO
}

bool BasisTwo::checkNewBasisOne()
{
    return true; // TODO
}
