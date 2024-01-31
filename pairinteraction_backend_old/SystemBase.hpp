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

#ifndef SYSTEMBASE_H
#define SYSTEMBASE_H

#include "State.hpp"
#include "serialization_eigen.hpp"
#include "utils.hpp"

#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index_container.hpp>
#include <cereal/types/complex.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/vector.hpp>

#include <algorithm>
#include <complex>
#include <functional>
#include <set>
#include <vector>

class MatrixElementCache;

template <class T>
class enumerated_state {
public:
    enumerated_state(size_t idx, T state) : idx(idx), state(std::move(state)) {}
    enumerated_state()
        : state() { // TODO remove and use
                    // http://www.boost.org/doc/libs/1_46_1/libs/serialization/doc/serialization.html#constructors
                    // instead
    }
    size_t idx{0};
    T state;

private:
    ////////////////////////////////////////////////////////////////////
    /// Method for serialization ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    friend class cereal::access;

    template <class Archive>
    void serialize(Archive &ar, unsigned int /* version */) {
        ar &CEREAL_NVP(idx) & CEREAL_NVP(state);
    }
};

#ifndef SWIG

namespace utils {
template <class T>
struct hash<enumerated_state<T>> {
    size_t operator()(enumerated_state<T> const &s) const { return std::hash<T>{}(s.state); }
};
} // namespace utils

#endif

template <class T>
struct states_set {
    typedef typename boost::multi_index_container<
        enumerated_state<T>,
        boost::multi_index::indexed_by<
            boost::multi_index::random_access<>,
            boost::multi_index::hashed_unique<
                boost::multi_index::member<enumerated_state<T>, T, &enumerated_state<T>::state>,
                std::hash<T>>>>
        type;
};

namespace cereal {

template <class Archive, typename Value, typename IndexSpecifierList, typename Allocator>
void save(
    Archive &ar,
    boost::multi_index::multi_index_container<Value, IndexSpecifierList, Allocator> const &vector) {
    ar(make_size_tag(static_cast<size_type>(vector.size())));
    for (auto &&v : vector) {
        ar(v);
    }
}

template <class Archive, typename Value, typename IndexSpecifierList, typename Allocator>
void load(Archive &ar,
          boost::multi_index::multi_index_container<Value, IndexSpecifierList, Allocator> &vector) {
    size_type size;
    ar(make_size_tag(size));
    vector.reserve(static_cast<size_type>(size));
    for (size_type i = 0; i < size; ++i) {
        Value v;
        ar(v);
        vector.push_back(v);
    }
}

} // namespace cereal

template <class Scalar_, class State_>
class SystemBase {
public:
    using Scalar = Scalar_;
    using State = State_;

    virtual ~SystemBase() = default;

    void setMinimalNorm(const double &threshold);

    ////////////////////////////////////////////////////////////////////
    /// Methods to restrict the number of states inside the basis //////
    ////////////////////////////////////////////////////////////////////

    void restrictEnergy(double e_min, double e_max);

    void restrictN(int n_min, int n_max);

    void restrictN(std::set<int> n);

    void restrictL(int l_min, int l_max);

    void restrictL(std::set<int> l);

    void restrictJ(float j_min, float j_max);

    void restrictJ(std::set<float> j);

    void restrictM(float m_min, float m_max);

    void restrictM(std::set<float> m);

    ////////////////////////////////////////////////////////////////////
    /// Method for adding user-defined states //////////////////////////
    ////////////////////////////////////////////////////////////////////

    void addStates(const State &s);

    void addStates(const std::set<State> &s);

    // TODO make it possible to just use added states, i.e. use no restrictions on quantum numbers
    // and energies

    ////////////////////////////////////////////////////////////////////
    /// Methods to get overlaps ////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    Eigen::VectorX<double> getOverlap(const State &generalizedstate);

    Eigen::VectorX<double> getOverlap(const std::vector<State> &generalizedstates);

    Eigen::VectorX<double> getOverlap(const size_t &state_index);

    Eigen::VectorX<double> getOverlap(const std::vector<size_t> &states_indices);

    Eigen::VectorX<double> getOverlap(const State &generalizedstate,
                                      std::array<double, 3> to_z_axis,
                                      std::array<double, 3> to_y_axis);

    Eigen::VectorX<double> getOverlap(const std::vector<State> &generalizedstates,
                                      std::array<double, 3> to_z_axis,
                                      std::array<double, 3> to_y_axis);

    Eigen::VectorX<double> getOverlap(const size_t &state_index, std::array<double, 3> to_z_axis,
                                      std::array<double, 3> to_y_axis);

    Eigen::VectorX<double> getOverlap(const std::vector<size_t> &states_indices,
                                      std::array<double, 3> to_z_axis,
                                      std::array<double, 3> to_y_axis);

    Eigen::VectorX<double> getOverlap(const State &generalizedstate, double alpha, double beta,
                                      double gamma);

    Eigen::VectorX<double> getOverlap(const size_t &state_index, double alpha, double beta,
                                      double gamma);

    Eigen::VectorX<double> getOverlap(const std::vector<State> &generalizedstates, double alpha,
                                      double beta, double gamma);

    Eigen::VectorX<double> getOverlap(const std::vector<size_t> &states_indices, double alpha,
                                      double beta, double gamma);

    ////////////////////////////////////////////////////////////////////
    /// Methods to get properties of the system ////////////////////////
    ////////////////////////////////////////////////////////////////////

    std::vector<State> getStates();

    const typename states_set<State>::type &getStatesMultiIndex();

    Eigen::SparseMatrix<Scalar> &getBasisvectors();

    Eigen::SparseMatrix<Scalar> &getHamiltonian();

    size_t getNumBasisvectors();

    size_t getNumStates();

    std::vector<State> getMainStates();

    std::array<std::vector<size_t>, 2> getConnections(SystemBase<Scalar, State> &system_to,
                                                      double threshold);

    ////////////////////////////////////////////////////////////////////
    /// Methods to build, transform, and destroy the system ////////////
    ////////////////////////////////////////////////////////////////////

    void buildHamiltonian();

    void buildInteraction();

    void buildBasis();

    void diagonalize(double energy_lower_bound, double energy_upper_bound);

    void diagonalize(double energy_lower_bound, double energy_upper_bound, double threshold);

    void diagonalize();

    void diagonalize(double threshold);

    void canonicalize();

    void unitarize();

    void rotate(std::array<double, 3> to_z_axis, std::array<double, 3> to_y_axis);

    void rotate(double alpha, double beta, double gamma);

    void add(SystemBase<Scalar, State> &system);

    void constrainBasisvectors(std::vector<size_t> indices_of_wanted_basisvectors);

    void applySchriefferWolffTransformation(SystemBase<Scalar, State> &system0);

    ////////////////////////////////////////////////////////////////////
    /// Methods to manipulate entries of the Hamiltonian ///////////////
    ////////////////////////////////////////////////////////////////////

    size_t getStateIndex(const State &searched_state);

    std::vector<size_t> getStateIndex(const std::vector<State> &searched_states);

    size_t getBasisvectorIndex(const State &searched_state);

    std::vector<size_t> getBasisvectorIndex(const std::vector<State> &searched_states);

    void forgetStatemixing();

    Scalar getHamiltonianEntry(const State &state_row, const State &state_col);

    void setHamiltonianEntry(const State &state_row, const State &state_col, Scalar value);

    void addHamiltonianEntry(const State &state_row, const State &state_col, Scalar value);

protected:
    SystemBase();

    SystemBase(MatrixElementCache &cache);

    SystemBase(MatrixElementCache &cache, bool memory_saving);

    virtual void initializeBasis() = 0;
    virtual void initializeInteraction() = 0;

    virtual void transformInteraction(const Eigen::SparseMatrix<Scalar> &transformator) = 0;
    virtual void addInteraction() = 0;
    virtual void deleteInteraction() = 0;

    virtual Eigen::SparseMatrix<Scalar> rotateStates(const std::vector<size_t> &states_indices,
                                                     double alpha, double beta, double gamma) = 0;
    virtual Eigen::SparseMatrix<Scalar> buildStaterotator(double alpha, double beta,
                                                          double gamma) = 0;
    virtual void incorporate(SystemBase<Scalar, State> &system) = 0;

    virtual void onStatesChange(){};

    MatrixElementCache *m_cache{nullptr};

    double threshold_for_sqnorm{0.05};

    double energy_min, energy_max;
    std::set<int> range_n, range_l;
    std::set<float> range_j, range_m;
    std::set<State> states_to_add;

    bool memory_saving{false};
    bool is_interaction_already_contained{false};
    bool is_new_hamiltonian_required{false};

    typename states_set<State>::type states;
    Eigen::SparseMatrix<Scalar> basisvectors;
    Eigen::SparseMatrix<Scalar> hamiltonian;
    Eigen::SparseMatrix<Scalar> basisvectors_unperturbed_cache;
    Eigen::SparseMatrix<Scalar> hamiltonian_unperturbed_cache;

    ////////////////////////////////////////////////////////////////////
    /// Helper method that shoul be called by the derived classes //////
    ////////////////////////////////////////////////////////////////////

    void onParameterChange();

    void onSymmetryChange();

    ////////////////////////////////////////////////////////////////////
    /// Helper methods to check diagonality and unitarity of a matrix //
    ////////////////////////////////////////////////////////////////////

    bool checkIsDiagonal(const Eigen::SparseMatrix<Scalar> &mat);

    bool checkIsUnitary(const Eigen::SparseMatrix<Scalar> &mat);

    ////////////////////////////////////////////////////////////////////
    /// Helper methods to check the validity of states /////////////////
    ////////////////////////////////////////////////////////////////////

    template <class V>
    bool checkIsQuantumnumberValid(V q, std::set<V> range_q) {
        return range_q.empty() || range_q.find(q) != range_q.end();
    }

    bool checkIsQuantumstateValid(const State &state);

    bool checkIsQuantumstateValid(const StateOne &state, bool a);

    bool checkIsQuantumstateValid(const StateTwo &state, std::array<bool, 2> a);

    bool checkIsEnergyValid(double e);

    ////////////////////////////////////////////////////////////////////
    /// Helper methods to change the set of basis vectors //////////////
    ////////////////////////////////////////////////////////////////////

    void applyLeftsideTransformator(std::vector<Eigen::Triplet<Scalar>> &triplets_transformator);

    void applyLeftsideTransformator(Eigen::SparseMatrix<Scalar> &transformator);

    void applyRightsideTransformator(std::vector<Eigen::Triplet<Scalar>> &triplets_transformator);

    void applyRightsideTransformator(Eigen::SparseMatrix<Scalar> &transformator);

    void
    removeRestrictedStates(std::function<bool(const enumerated_state<State> &)> checkIsValidEntry);

    ////////////////////////////////////////////////////////////////////
    /// Utility methods ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    template <class V>
    void range(std::set<V> &rset, V rmin, V rmax) {
        rset.clear();
        for (V r = rmin; r <= rmax; ++r) {
            rset.insert(r);
        }
    }

    Eigen::Matrix<double, 3, 3> buildRotator(std::array<double, 3> to_z_axis,
                                             std::array<double, 3> to_y_axis);

    Eigen::Matrix<double, 3, 3> buildRotator(const double &alpha, const double &beta,
                                             const double &gamma);

    Eigen::Matrix<double, 3, 1> getEulerAngles(const std::array<double, 3> &to_z_axis,
                                               const std::array<double, 3> &to_y_axis);

private:
    void forgetRestrictions();

    ////////////////////////////////////////////////////////////////////
    /// Method to update the system ////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    void updateEverything();

    ////////////////////////////////////////////////////////////////////
    /// Method for serialization ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    friend class cereal::access;

    template <class Archive>
    void serialize(Archive &ar, unsigned int /*version*/) {
        ar &cereal::make_nvp("cache", *m_cache) & CEREAL_NVP(threshold_for_sqnorm);
        ar &CEREAL_NVP(energy_min) & CEREAL_NVP(energy_max) & CEREAL_NVP(range_n) &
            CEREAL_NVP(range_l) & CEREAL_NVP(range_j) & CEREAL_NVP(range_m) &
            CEREAL_NVP(states_to_add);
        ar &CEREAL_NVP(memory_saving) & CEREAL_NVP(is_interaction_already_contained) &
            CEREAL_NVP(is_new_hamiltonian_required);
        ar &CEREAL_NVP(states) & CEREAL_NVP(basisvectors) & CEREAL_NVP(hamiltonian);
        ar &CEREAL_NVP(basisvectors_unperturbed_cache) & CEREAL_NVP(hamiltonian_unperturbed_cache);
    }
};

extern template class SystemBase<std::complex<double>, StateOne>;
extern template class SystemBase<std::complex<double>, StateTwo>;
extern template class SystemBase<double, StateOne>;
extern template class SystemBase<double, StateTwo>;

#endif
