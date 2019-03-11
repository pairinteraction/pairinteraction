/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke, Johannes Block. All rights reserved.
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

#ifndef SYSTEMTWO_H
#define SYSTEMTWO_H

#include "State.h"
#include "SystemBase.h"
#include "SystemOne.h"

#include <Eigen/Sparse>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <cmath>
#include <set>
#include <type_traits>

class SystemTwo : public SystemBase<StateTwo> {
public:
    SystemTwo(const SystemOne &b1, const SystemOne &b2, MatrixElementCache &cache);
    SystemTwo(const SystemOne &b1, const SystemOne &b2, MatrixElementCache &cache,
              bool memory_saving);

    const std::array<std::string, 2> &getSpecies();
    std::vector<StateOne> getStatesFirst();
    std::vector<StateOne> getStatesSecond();
    void enableGreenTensor(bool GTboolean);
    void setSurfaceDistance(double d);
    void setAngle(double a);
    void setDistance(double d);
    void setDistanceVector(std::array<double, 3> d);
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
    std::unordered_map<int, eigen_sparse_t> interaction_greentensor_dd;
    std::unordered_map<int, eigen_sparse_t> interaction_greentensor_dq;
    std::unordered_map<int, eigen_sparse_t> interaction_greentensor_qd;

    double minimal_le_roy_radius;
    double distance, distance_x, distance_y, distance_z;
    bool GTbool;
    double surface_distance;
    unsigned int ordermax;

    parity_t sym_permutation;
    parity_t sym_inversion;
    parity_t sym_reflection;
    std::set<int> sym_rotation;

    std::unordered_map<int, double> angle_terms;
    std::unordered_map<int, double> greentensor_terms_dd;
    std::unordered_map<int, double> greentensor_terms_dq;
    std::unordered_map<int, double> greentensor_terms_qd;

    ////////////////////////////////////////////////////////////////////
    /// Utility methods ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    void checkDistance(const double &distance);

    void addBasisvectors(const StateTwo &state, const size_t &col_new, const scalar_t &value_new,
                         std::vector<eigen_triplet_t> &basisvectors_triplets,
                         std::vector<double> &sqnorm_list);

    template <typename T>
    void addTriplet(std::vector<Eigen::Triplet<T>> &triplets, size_t r_idx, size_t c_idx, T val) {
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
        ar &distance &distance_x &distance_y &distance_z &surface_distance &ordermax;
        ar &sym_permutation &sym_inversion &sym_reflection &sym_rotation;
        ar &angle_terms &greentensor_terms_dd &greentensor_terms_dq &greentensor_terms_qd;
        ar &interaction_angulardipole &interaction_multipole &interaction_greentensor_dd
            &interaction_greentensor_dq &interaction_greentensor_qd;
    }
};

#endif
