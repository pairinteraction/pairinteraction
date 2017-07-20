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

#ifndef SYSTEMONE_H
#define SYSTEMONE_H

#include "State.h"
#include "SystemBase.h"

#include <cmath>
#include <type_traits>
#include <array>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/unordered_map.hpp>

#ifndef SWIG
namespace std {
    template <>
    struct hash<std::array<int, 2>>
    {
        size_t operator()(const std::array<int, 2>& a) const
        {
            return boost::hash_value(a);
        }
    };
}
#endif

class SystemOne : public SystemBase<StateOne> {
public:
    SystemOne(std::wstring const& element, std::wstring cachedir);
    SystemOne(std::wstring const& element, std::wstring cachedir, bool memory_saving);
    SystemOne(std::wstring const& element);
    SystemOne(std::wstring const& element, bool memory_saving);

    const std::wstring& getElement() const;
    void setEfield(std::array<double, 3> field);
    void setBfield(std::array<double, 3> field);
    void setEfield(std::array<double, 3> field, std::array<double, 3> to_z_axis, std::array<double, 3> to_y_axis);
    void setBfield(std::array<double, 3> field, std::array<double, 3> to_z_axis, std::array<double, 3> to_y_axis);
    void setEfield(std::array<double, 3> field, double alpha, double beta, double gamma);
    void setBfield(std::array<double, 3> field, double alpha, double beta, double gamma);
    void setDiamagnetism(bool enable);

protected:
    void initializeBasis() override;
    void initializeInteraction() override;
    void addInteraction() override;
    void transformInteraction(const eigen_sparse_t &transformator) override;
    void deleteInteraction() override;
    eigen_sparse_t rotateStates(const std::vector<size_t> &states_indices, double alpha, double beta, double gamma) override;
    eigen_sparse_t buildStaterotator(double alpha, double beta, double gamma) override;
    void incorporate(SystemBase<StateOne> &system) override;

private:
    std::array<double, 3> efield, bfield;
    std::unordered_map<int, scalar_t>  efield_spherical, bfield_spherical;
    bool diamagnetism;
    std::unordered_map<std::array<int, 2>, scalar_t> diamagnetism_terms;
    std::wstring element;

    std::unordered_map<int, eigen_sparse_t> interaction_efield;
    std::unordered_map<int, eigen_sparse_t> interaction_bfield;
    std::unordered_map<std::array<int, 2>, eigen_sparse_t> interaction_diamagnetism;

    void changeToSphericalbasis(std::array<double, 3> field, std::unordered_map<int, double>& field_spherical);
    void changeToSphericalbasis(std::array<double, 3> field, std::unordered_map<int, std::complex<double>>& field_spherical);
    void addTriplet(std::vector<eigen_triplet_t> &triplets, const size_t r_idx, const size_t c_idx, const scalar_t val);
    void rotateVector(std::array<double, 3> &field, std::array<double, 3> &to_z_axis, std::array<double, 3> &to_y_axis);
    void rotateVector(std::array<double, 3> &field, double alpha, double beta, double gamma);

    ////////////////////////////////////////////////////////////////////
    /// Utility methods ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    template<class T>
    void addRotated(const StateOne &state, const size_t &idx, std::vector<Eigen::Triplet<T>> &triplets, WignerD &wigner, const double &alpha, const double &beta, const double &gamma) {
        // Check whether the angles are compatible to the used data type
        double tolerance = 1e-16;
        if (std::is_same<T, double>::value && std::abs(std::remainder(alpha,2*M_PI)) > tolerance) throw std::runtime_error( "If the Euler angle alpha is not a multiple of 2 pi, the Wigner D matrix element is complex and cannot be converted to double." );
        if (std::is_same<T, double>::value && std::abs(std::remainder(gamma,2*M_PI)) > tolerance) throw std::runtime_error( "If the Euler angle gamma is not a multiple of 2 pi, the Wigner D matrix element is complex and cannot be converted to double." );

        // Add rotated triplet entries
        StateOne newstate = state;

        for (float m = -state.j; m <= state.j; ++m) {
            newstate.m = m;
            auto state_iter = states.get<1>().find(newstate);

            if (state_iter != states.get<1>().end()) {
                T val = convert<T>(wigner(state.j, state.m, m, alpha, beta, gamma));
                triplets.push_back(Eigen::Triplet<T>(state_iter->idx, idx, val));
            } else {
                std::cerr << "Warning: Incomplete rotation because the basis is lacking some Zeeman levels." << std::endl;
            }
        }
    }

    template<class T, class S>
    T convert(const S &val) {
        return val;
    }

    ////////////////////////////////////////////////////////////////////
    /// Method for serialization ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        (void)version;

        ar & boost::serialization::base_object<SystemBase<StateOne>>(*this);
        ar & element;
        ar & efield & bfield & diamagnetism;
        ar & efield_spherical & bfield_spherical & diamagnetism_terms;
        ar & interaction_efield & interaction_bfield & interaction_diamagnetism;
    }
};

#endif
