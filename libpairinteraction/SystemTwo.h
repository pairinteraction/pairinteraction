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

#include <boost/math/special_functions/binomial.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <cmath>
#include <set>
#include <type_traits>
#include <Eigen/Sparse>

class SystemTwo : public SystemBase<StateTwo> {
public:
    SystemTwo(const SystemOne &b1, const SystemOne &b2, MatrixElementCache &cache);
    SystemTwo(const SystemOne &b1, const SystemOne &b2, MatrixElementCache &cache,
              bool memory_saving);

    const std::array<std::string, 2> &getSpecies();
    std::vector<StateOne> getStatesFirst();
    std::vector<StateOne> getStatesSecond();
    void setDistance(double d);
    void setGTbool(bool GTboolean);
    void setSurfaceDistance(double d);
    void setAngle(double a);
    void setOrder(double o);

    void setConservedParityUnderPermutation(parity_t parity);
    void setConservedParityUnderInversion(parity_t parity);
    void setConservedParityUnderReflection(parity_t parity);
    void setConservedMomentaUnderRotation(const std::set<int> &momenta);

protected:
    void initializeBasis() override;
    void initializeInteraction() override;
    void addInteraction() override;
    void transformInteraction(const eigen_sparse_t &transformator) override;
    void deleteInteraction() override;
    eigen_sparse_t rotateStates(const std::vector<size_t> &states_indices, double alpha,
                                double beta, double gamma) override;
    eigen_sparse_t buildStaterotator(double alpha, double beta, double gamma) override;
    void incorporate(SystemBase<StateTwo> &system) override;
    void onStatesChange() override;

private:
    std::array<std::string, 2> species;
    SystemOne system1; // is needed in the initializeBasis method and afterwards deleted
    SystemOne system2; // is needed in the initializeBasis method and afterwards deleted

    std::unordered_map<int, eigen_sparse_t> interaction_angulardipole;
    std::unordered_map<int, eigen_sparse_t> interaction_multipole;

    double minimal_le_roy_radius;
    double distance;
    bool GTbool;
    double surface_distance;
    double angle;
    unsigned int ordermax;    
    
    std::unordered_map<std::string, eigen_sparse_t > dipoleMatrices; //Dipole matrices for GreenTensor method
    

    parity_t sym_permutation;
    parity_t sym_inversion;
    parity_t sym_reflection;
    std::set<int> sym_rotation;

    std::array<double, 4> angle_terms;

    ////////////////////////////////////////////////////////////////////
    /// Utility methods ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    void checkDistance(const double &distance);

    void addBasisvectors(const StateTwo &state, const size_t &col_new, const scalar_t &value_new,
                         std::vector<eigen_triplet_t> &basisvectors_triplets,
                         std::vector<double> &sqnorm_list);

    template <typename T>
    void addTriplet(std::vector<Eigen::Triplet<T> > &triplets, size_t r_idx, size_t c_idx,
                    T val) {
            triplets.emplace_back(r_idx, c_idx, val);
    }

    template <class T>
    void addRotated(const StateTwo &state, const size_t &idx,
                    std::vector<Eigen::Triplet<T>> &triplets, WignerD &wigner, const double &alpha,
                    const double &beta, const double &gamma) {
        // Check whether the angles are compatible to the used data type
        double tolerance = 1e-16;
        if (std::is_same<T, double>::value &&
            std::abs(std::remainder(alpha, 2 * M_PI)) > tolerance) {
            throw std::runtime_error(
                "If the Euler angle alpha is not a multiple of 2 pi, the Wigner D matrix element "
                "is complex and cannot be converted to double.");
        }
        if (std::is_same<T, double>::value &&
            std::abs(std::remainder(gamma, 2 * M_PI)) > tolerance) {
            throw std::runtime_error(
                "If the Euler angle gamma is not a multiple of 2 pi, the Wigner D matrix element "
                "is complex and cannot be converted to double.");
        }

        // Add rotated triplet entries
        std::vector<T> val2_vector;
        val2_vector.reserve(2 * state.getSecondState().getJ() + 1);

        for (float m2 = -state.getSecondState().getJ(); m2 <= state.getSecondState().getJ(); ++m2) {
            val2_vector.push_back(
                utils::convert<T>(wigner(state.getSecondState().getJ(), m2,
                                         state.getSecondState().getM(), -gamma, -beta, -alpha)));
        }

        for (float m1 = -state.getFirstState().getJ(); m1 <= state.getFirstState().getJ(); ++m1) {
            auto val1 =
                utils::convert<T>(wigner(state.getFirstState().getJ(), m1,
                                         state.getFirstState().getM(), -gamma, -beta, -alpha));

            for (float m2 = -state.getSecondState().getJ(); m2 <= state.getSecondState().getJ();
                 ++m2) {
                StateTwo newstate(state.getSpecies(), state.getN(), state.getL(), state.getJ(),
                                  {{m1, m2}});
                auto state_iter = states.get<1>().find(newstate);

                if (state_iter != states.get<1>().end()) {
                    auto val = val1 * val2_vector[m2 + state.getSecondState().getJ()];
                    triplets.push_back(Eigen::Triplet<T>(state_iter->idx, idx, val));
                } else {
                    std::cerr << "Warning: Incomplete rotation because the basis is lacking some "
                                 "Zeeman levels."
                              << std::endl;
                }
            }
        }
    }

    bool isRefelectionAndRotationCompatible();

    ////////////////////////////////////////////////////////////////////
    /// Method for serialization ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int /*version*/) {
        ar &boost::serialization::base_object<SystemBase<StateTwo>>(*this);
        ar &species &system1 &system2;
        ar &distance &angle &surface_distance &ordermax &sym_permutation &sym_inversion &sym_reflection &sym_rotation;
        ar &angle_terms; //TODO include GreenTensor, dipoleMatrices
        ar &interaction_angulardipole &interaction_multipole;
    }
};

#endif
