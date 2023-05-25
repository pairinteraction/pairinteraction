/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
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

#ifndef SYSTEMONE_H
#define SYSTEMONE_H

#include "State.hpp"
#include "Symmetry.hpp"
#include "SystemBase.hpp"
#include "WignerD.hpp"
#include "utils.hpp"

#include <array>
// clang-format off
#if __has_include (<boost/serialization/version.hpp>)
#    include <boost/serialization/version.hpp>
#endif
#if __has_include (<boost/serialization/library_version_type.hpp>)
#    include <boost/serialization/library_version_type.hpp>
#endif
// clang-format on
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <cmath>
#include <set>
#include <unordered_map>

template <typename Scalar_>
class SystemOne : public SystemBase<Scalar_, StateOne> {
public:
    using Scalar = Scalar_;
    SystemOne(std::string species, MatrixElementCache &cache);
    SystemOne(std::string species, MatrixElementCache &cache, bool memory_saving);

    const std::string &getSpecies() const;
    void setEfield(std::array<double, 3> field);
    void setBfield(std::array<double, 3> field);
    void setEfield(std::array<double, 3> field, std::array<double, 3> to_z_axis,
                   std::array<double, 3> to_y_axis);
    void setBfield(std::array<double, 3> field, std::array<double, 3> to_z_axis,
                   std::array<double, 3> to_y_axis);
    void setEfield(std::array<double, 3> field, double alpha, double beta, double gamma);
    void setBfield(std::array<double, 3> field, double alpha, double beta, double gamma);
    void enableDiamagnetism(bool enable);
    void setIonCharge(int c);
    void setRydIonOrder(unsigned int o);
    void setRydIonDistance(double d);
    void setConservedParityUnderReflection(parity_t parity);
    void setConservedMomentaUnderRotation(const std::set<float> &momenta);

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
    void incorporate(SystemBase<Scalar_, StateOne> &system) override;

private:
    std::array<double, 3> efield, bfield;
    std::unordered_map<int, Scalar_> efield_spherical, bfield_spherical;
    bool diamagnetism;
    std::unordered_map<std::array<int, 2>, Scalar_, utils::hash<std::array<int, 2>>>
        diamagnetism_terms;
    int charge;
    unsigned int ordermax;
    double distance;
    std::string species;

    std::unordered_map<int, Eigen::SparseMatrix<Scalar_>> interaction_efield;
    std::unordered_map<int, Eigen::SparseMatrix<Scalar_>> interaction_bfield;
    std::unordered_map<std::array<int, 2>, Eigen::SparseMatrix<Scalar_>,
                       utils::hash<std::array<int, 2>>>
        interaction_diamagnetism;
    std::unordered_map<int, Eigen::SparseMatrix<Scalar_>> interaction_multipole;
    parity_t sym_reflection;
    std::set<float> sym_rotation;

    ////////////////////////////////////////////////////////////////////
    /// Utility methods ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    void addSymmetrizedBasisvectors(const StateOne &state, size_t &idx, const double &energy,
                                    std::vector<Eigen::Triplet<Scalar_>> &basisvectors_triplets,
                                    std::vector<Eigen::Triplet<Scalar_>> &hamiltonian_triplets,
                                    parity_t &sym_reflection_local);

    void addBasisvectors(const StateOne &state, const size_t &idx, const Scalar_ &value,
                         std::vector<Eigen::Triplet<Scalar_>> &basisvectors_triplets);

    void changeToSphericalbasis(std::array<double, 3> field,
                                std::unordered_map<int, double> &field_spherical);
    void changeToSphericalbasis(std::array<double, 3> field,
                                std::unordered_map<int, std::complex<double>> &field_spherical);
    void addTriplet(std::vector<Eigen::Triplet<Scalar_>> &triplets, size_t r_idx, size_t c_idx,
                    Scalar_ val);
    void rotateVector(std::array<double, 3> &field, std::array<double, 3> &to_z_axis,
                      std::array<double, 3> &to_y_axis);
    void rotateVector(std::array<double, 3> &field, double alpha, double beta, double gamma);

    template <class T>
    void addRotated(const StateOne &state, const size_t &idx,
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
        for (float m = -state.getJ(); m <= state.getJ(); ++m) {
            StateOne newstate(state.getSpecies(), state.getN(), state.getL(), state.getJ(), m);
            auto state_iter = this->states.template get<1>().find(newstate);

            if (state_iter != this->states.template get<1>().end()) {
                // We calculate the matrix element <m|d_mM|state.getM()>|m>. Note that we must
                // invert the angles since they are given for a passive rotation, i.e. rotation of
                // the coordinate system, and the WignerD matrix is for an active rotation of spins.
                auto val =
                    utils::convert<T>(wigner(state.getJ(), m, state.getM(), -gamma, -beta, -alpha));
                triplets.push_back(Eigen::Triplet<T>(state_iter->idx, idx, val));
            } else {
                std::cerr << "Warning: Incomplete rotation because the basis is lacking some "
                             "Zeeman levels."
                          << std::endl;
            }
        }
    }

    bool isRefelectionAndRotationCompatible();

    ////////////////////////////////////////////////////////////////////
    /// Method for serialization ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, unsigned int /*version*/) {
        ar &boost::serialization::base_object<SystemBase<Scalar_, StateOne>>(*this);
        ar &species;
        ar &efield &bfield &diamagnetism &sym_reflection &sym_rotation;
        ar &efield_spherical &bfield_spherical &diamagnetism_terms;
        ar &interaction_efield &interaction_bfield &interaction_diamagnetism;
    }
};

extern template class SystemOne<std::complex<double>>;
extern template class SystemOne<double>;

#endif
