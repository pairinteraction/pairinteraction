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

#include "State.hpp"
#include "SystemBase.hpp"
#include "SystemOne.hpp"

#include "EigenCompat.hpp"
#include <Eigen/SparseCore>
#include <boost/math/special_functions/binomial.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>

#include <cmath>
#include <set>
#include <type_traits>

template <typename Scalar_>
class SystemTwo : public SystemBase<Scalar_, StateTwo> {
public:
    using Scalar = Scalar_;
    SystemTwo(const SystemOne<Scalar_> &b1, const SystemOne<Scalar_> &b2,
              MatrixElementCache &cache);
    SystemTwo(const SystemOne<Scalar_> &b1, const SystemOne<Scalar_> &b2, MatrixElementCache &cache,
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
    void setOneAtomBasisvectors(const std::vector<std::array<size_t, 2>> &indices);

protected:
    void initializeBasis() override;
    void initializeInteraction() override;
    void addInteraction() override;
    void transformInteraction(const Eigen::SparseMatrix<Scalar_> &transformator) override;
    void deleteInteraction() override;
    Eigen::SparseMatrix<Scalar_> rotateStates(const std::vector<size_t> &states_indices,
                                              double alpha, double beta, double gamma) override;
    Eigen::SparseMatrix<Scalar_> buildStaterotator(double alpha, double beta,
                                                   double gamma) override;
    void incorporate(SystemBase<Scalar_, StateTwo> &system) override;
    void onStatesChange() override;

private:
    std::array<std::string, 2> species;
    SystemOne<Scalar_> system1; // is needed in the initializeBasis method and afterwards deleted
    SystemOne<Scalar_> system2; // is needed in the initializeBasis method and afterwards deleted

    std::unordered_map<int, Eigen::SparseMatrix<Scalar_>> interaction_angulardipole;
    std::unordered_map<int, Eigen::SparseMatrix<Scalar_>> interaction_multipole;
    std::unordered_map<int, Eigen::SparseMatrix<Scalar_>> interaction_greentensor_dd;
    std::unordered_map<int, Eigen::SparseMatrix<Scalar_>> interaction_greentensor_dq;
    std::unordered_map<int, Eigen::SparseMatrix<Scalar_>> interaction_greentensor_qd;

    double minimal_le_roy_radius;
    double distance;
    double distance_x; // NOLINT
    double distance_y; // NOLINT
    double distance_z;
    bool GTbool; // NOLINT
    double surface_distance;
    unsigned int ordermax; // NOLINT

    parity_t sym_permutation; // NOLINT
    parity_t sym_inversion;   // NOLINT
    parity_t sym_reflection;  // NOLINT
    std::set<int> sym_rotation;

    std::unordered_map<int, double> angle_terms;
    std::unordered_map<int, double> greentensor_terms_dd;
    std::unordered_map<int, double> greentensor_terms_dq;
    std::unordered_map<int, double> greentensor_terms_qd;
    std::vector<std::array<size_t, 2>> one_atom_basisvectors_indices;

    ////////////////////////////////////////////////////////////////////
    /// Utility methods ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    void checkDistance(const double &distance);

    void addBasisvectors(const StateTwo &state, const size_t &col_new, const Scalar_ &value_new,
                         std::vector<Eigen::Triplet<Scalar_>> &basisvectors_triplets,
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
                auto state_iter = this->states.template get<1>().find(newstate);

                if (state_iter != this->states.template get<1>().end()) {
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

    friend class cereal::access;
    SystemTwo();

    template <class Archive>
    void serialize(Archive &ar, unsigned int /* version */) {
        ar &cereal::make_nvp("base_class", cereal::base_class<SystemBase<Scalar_, StateTwo>>(this));
        ar &CEREAL_NVP(species) & CEREAL_NVP(system1) & CEREAL_NVP(system2);
        ar &CEREAL_NVP(distance) & CEREAL_NVP(distance_x) & CEREAL_NVP(distance_y) &
            CEREAL_NVP(distance_z) & CEREAL_NVP(surface_distance) & CEREAL_NVP(ordermax);
        ar &CEREAL_NVP(sym_permutation) & CEREAL_NVP(sym_inversion) & CEREAL_NVP(sym_reflection) &
            CEREAL_NVP(sym_rotation);
        ar &CEREAL_NVP(angle_terms) & CEREAL_NVP(greentensor_terms_dd) &
            CEREAL_NVP(greentensor_terms_dq) & CEREAL_NVP(greentensor_terms_qd);
        ar &CEREAL_NVP(interaction_angulardipole) & CEREAL_NVP(interaction_multipole) &
            CEREAL_NVP(interaction_greentensor_dd) & CEREAL_NVP(interaction_greentensor_dq) &
            CEREAL_NVP(interaction_greentensor_qd);
    }
};

extern template class SystemTwo<std::complex<double>>;
extern template class SystemTwo<double>;

#ifndef SWIG
CEREAL_REGISTER_TYPE(SystemTwo<std::complex<double>>) // NOLINT
CEREAL_REGISTER_TYPE(SystemTwo<double>)               // NOLINT
#endif

#endif
