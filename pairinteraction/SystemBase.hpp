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

#include "MatrixElementCache.hpp"
#include "State.hpp"
#include "WignerD.hpp"
#include "dtypes.hpp"
#include "serialization_eigen.hpp"
#include "serialization_path.hpp"
#include "utils.hpp"
#include <unsupported/Eigen/MatrixFunctions>

#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index_container.hpp>
// clang-format off
#if __has_include (<boost/serialization/version.hpp>)
#    include <boost/serialization/version.hpp>
#endif
#if __has_include (<boost/serialization/library_version_type.hpp>)
#    include <boost/serialization/library_version_type.hpp>
#endif
// clang-format on
#include <boost/serialization/complex.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>
#include <complex>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef WITH_INTEL_MKL
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include <mkl_lapacke.h>
#include <mkl_solvers_ee.h>
#elif defined EIGEN_USE_LAPACKE
#include <lapacke.h>
#endif // WITH_INTEL_MKL

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

    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        (void)version;

        ar &idx &state;
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

template <class Scalar, class State>
class SystemBase {
public:
    virtual ~SystemBase() = default;

    void setMinimalNorm(const double &threshold) { threshold_for_sqnorm = threshold; }

    ////////////////////////////////////////////////////////////////////
    /// Methods to restrict the number of states inside the basis //////
    ////////////////////////////////////////////////////////////////////

    void restrictEnergy(double e_min, double e_max) { // restricts diagonal entries
        energy_min = e_min;
        energy_max = e_max;
    }

    void restrictN(int n_min, int n_max) { this->range(range_n, n_min, n_max); }

    void restrictN(std::set<int> n) { range_n = n; }

    void restrictL(int l_min, int l_max) { this->range(range_l, l_min, l_max); }

    void restrictL(std::set<int> l) { range_l = l; }

    void restrictJ(float j_min, float j_max) { this->range(range_j, j_min, j_max); }

    void restrictJ(std::set<float> j) { range_j = j; }

    void restrictM(float m_min, float m_max) { this->range(range_m, m_min, m_max); }

    void restrictM(std::set<float> m) { range_m = m; }

    ////////////////////////////////////////////////////////////////////
    /// Method for adding user-defined states //////////////////////////
    ////////////////////////////////////////////////////////////////////

    void addStates(const State &s) { states_to_add.insert(s); }

    void addStates(const std::set<State> &s) { states_to_add.insert(s.begin(), s.end()); }

    // TODO make it possible to just use added states, i.e. use no restrictions on quantum numbers
    // and energies

    ////////////////////////////////////////////////////////////////////
    /// Methods to get overlaps ////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    eigen_vector_double_t getOverlap(const State &generalizedstate) {
        return this->getOverlap(generalizedstate, 0, 0, 0);
    }

    eigen_vector_double_t getOverlap(const std::vector<State> &generalizedstates) {
        return this->getOverlap(generalizedstates, 0, 0, 0);
    }

    eigen_vector_double_t getOverlap(const size_t &state_index) {
        return this->getOverlap(state_index, 0, 0, 0);
    }

    eigen_vector_double_t getOverlap(const std::vector<size_t> &states_indices) {
        return this->getOverlap(states_indices, 0, 0, 0);
    }

    eigen_vector_double_t getOverlap(const State &generalizedstate, std::array<double, 3> to_z_axis,
                                     std::array<double, 3> to_y_axis) {
        auto euler_zyz = this->getEulerAngles(to_z_axis, to_y_axis);
        return this->getOverlap(generalizedstate, euler_zyz[0], euler_zyz[1], euler_zyz[2]);
    }

    eigen_vector_double_t getOverlap(const std::vector<State> &generalizedstates,
                                     std::array<double, 3> to_z_axis,
                                     std::array<double, 3> to_y_axis) {
        auto euler_zyz = this->getEulerAngles(to_z_axis, to_y_axis);
        return this->getOverlap(generalizedstates, euler_zyz[0], euler_zyz[1], euler_zyz[2]);
    }

    eigen_vector_double_t getOverlap(const size_t &state_index, std::array<double, 3> to_z_axis,
                                     std::array<double, 3> to_y_axis) {
        auto euler_zyz = this->getEulerAngles(to_z_axis, to_y_axis);
        return this->getOverlap(state_index, euler_zyz[0], euler_zyz[1], euler_zyz[2]);
    }

    eigen_vector_double_t getOverlap(const std::vector<size_t> &states_indices,
                                     std::array<double, 3> to_z_axis,
                                     std::array<double, 3> to_y_axis) {
        auto euler_zyz = this->getEulerAngles(to_z_axis, to_y_axis);
        return this->getOverlap(states_indices, euler_zyz[0], euler_zyz[1], euler_zyz[2]);
    }

    eigen_vector_double_t getOverlap(const State &generalizedstate, double alpha, double beta,
                                     double gamma) {
        std::vector<State> generalizedstates({generalizedstate});
        return this->getOverlap(generalizedstates, alpha, beta, gamma);
    }

    eigen_vector_double_t getOverlap(const size_t &state_index, double alpha, double beta,
                                     double gamma) {
        std::vector<size_t> states_indices({state_index});
        return this->getOverlap(states_indices, alpha, beta, gamma);
    }

    eigen_vector_double_t getOverlap(const std::vector<State> &generalizedstates, double alpha,
                                     double beta, double gamma) {
        // Build basis
        this->buildBasis();

        // Determine indices of the specified states
        std::vector<size_t> states_indices;
        states_indices.reserve(generalizedstates.size());

        for (const auto &searched_state : generalizedstates) {
            if (utils::is_true(searched_state.isGeneralized())) {
                for (const auto &state : states) {
                    if (state.state ^ searched_state) { // Check whether state.state is contained in
                                                        // searched_state
                        states_indices.push_back(state.idx);
                    }
                }
            } else {
                auto state_iter = states.template get<1>().find(searched_state);
                if (state_iter != states.template get<1>().end()) {
                    states_indices.push_back(state_iter->idx);
                }
            }
        }

        // Get the overlap
        return this->getOverlap(states_indices, alpha, beta, gamma);
    }

    eigen_vector_double_t getOverlap(const std::vector<size_t> &states_indices, double alpha,
                                     double beta, double gamma) {
        // Build basis
        this->buildBasis();

        // Build state vectors
        Eigen::SparseMatrix<Scalar> overlap_states;
        if (alpha == 0 && beta == 0 && gamma == 0) {
            std::vector<Eigen::Triplet<Scalar>> overlap_states_triplets;
            overlap_states_triplets.reserve(states_indices.size());

            size_t current = 0;
            for (auto const &idx : states_indices) {
                overlap_states_triplets.emplace_back(idx, current++, 1);
            }

            overlap_states.resize(states.size(), states_indices.size());
            overlap_states.setFromTriplets(overlap_states_triplets.begin(),
                                           overlap_states_triplets.end());
        } else {
            overlap_states = this->rotateStates(states_indices, alpha, beta, gamma);
        }

        // Calculate overlap
        Eigen::SparseMatrix<Scalar> product = basisvectors.adjoint() * overlap_states;
        eigen_vector_double_t overlap = eigen_vector_double_t::Zero(product.rows());
        for (int k = 0; k < product.outerSize(); ++k) {
            for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple(product, k); triple;
                 ++triple) {
                overlap[triple.row()] += std::pow(std::abs(triple.value()), 2);
            }
        }

        return overlap;
    }

    ////////////////////////////////////////////////////////////////////
    /// Methods to get properties of the system ////////////////////////
    ////////////////////////////////////////////////////////////////////

    std::vector<State> getStates() {
        this->buildBasis();
        std::vector<State> states_converted;
        states_converted.reserve(states.size());
        for (const auto &s : states) {
            states_converted.push_back(std::move(s.state));
        }
        return states_converted;
    }

    const typename states_set<State>::type &
    getStatesMultiIndex() { // TODO @hmenke typemap for "const typename states_set<T>::type &"
        return states;
    }

    Eigen::SparseMatrix<Scalar> &getBasisvectors() {
        this->buildBasis();
        return basisvectors;
    }

    Eigen::SparseMatrix<Scalar> &getHamiltonian() {
        this->buildHamiltonian();
        return hamiltonian;
    }

    size_t getNumBasisvectors() {
        // Build basis
        this->buildBasis();

        // Check variables for consistency
        if ((basisvectors.outerSize() != basisvectors.cols()) ||
            (basisvectors.outerSize() != hamiltonian.rows()) ||
            (basisvectors.outerSize() != hamiltonian.cols())) {
            throw std::runtime_error("Inconsistent variables at " + std::string(__FILE__) + ":" +
                                     std::to_string(__LINE__) + ".");
        }

        return basisvectors.cols();
    }

    size_t getNumStates() {
        // Build basis
        this->buildBasis();

        // Check variables for consistency
        if ((basisvectors.innerSize() != basisvectors.rows()) ||
            (static_cast<size_t>(basisvectors.innerSize()) != states.size())) {
            throw std::runtime_error("Inconsistent variables at " + std::string(__FILE__) + ":" +
                                     std::to_string(__LINE__) + ".");
        }

        return basisvectors.rows();
    }

    std::vector<State> getMainStates() {
        // Build basis
        this->buildBasis();

        // Get states with the main contribution
        std::vector<State> states_with_maxval;
        states_with_maxval.reserve(basisvectors.cols());

        for (int k = 0; k < basisvectors.outerSize(); ++k) { // col == idx_vector
            double maxval = -1;
            size_t row_with_maxval;

            for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple(basisvectors, k);
                 triple; ++triple) {
                if (std::abs(triple.value()) > maxval) {
                    row_with_maxval = triple.row();
                    maxval = std::abs(triple.value());
                }
            }

            states_with_maxval.push_back(states[row_with_maxval].state);
        }

        return states_with_maxval;
    }

    std::array<std::vector<size_t>, 2> getConnections(SystemBase<Scalar, State> &system_to,
                                                      double threshold) {
        // Build basis
        this->buildBasis();
        system_to.buildBasis();

        double threshold_sqrt = std::sqrt(threshold);

        // Calculate transformator between the set of states
        std::vector<Eigen::Triplet<Scalar>> triplets_transformator;
        triplets_transformator.reserve(std::min(states.size(), system_to.states.size()));

        for (const auto &s : system_to.states) {
            auto state_iter = states.template get<1>().find(s.state);
            if (state_iter != states.template get<1>().end()) {
                size_t idx_from = state_iter->idx;
                triplets_transformator.emplace_back(idx_from, s.idx, 1);
            }
        }

        Eigen::SparseMatrix<Scalar> transformator(states.size(), system_to.states.size());
        transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

        // Calculate overlap
        eigen_sparse_double_t overlap_sqrt =
            (basisvectors.adjoint() * transformator * system_to.basisvectors).cwiseAbs();

        // Determine the indices of the maximal values within the rows
        std::vector<size_t> rows_with_maxval;
        rows_with_maxval.reserve(overlap_sqrt.cols());

        for (int k = 0; k < overlap_sqrt.outerSize(); ++k) {
            double maxval = -1;
            size_t row_with_maxval;

            for (eigen_iterator_double_t triple(overlap_sqrt, k); triple; ++triple) {
                if (triple.value() > maxval) {
                    row_with_maxval = triple.row();
                    maxval = triple.value();
                }
            }

            rows_with_maxval.push_back(row_with_maxval);
        }

        // Determine the maximal values within the columns and construct connections
        eigen_sparse_double_t overlap_sqrt_transposed = overlap_sqrt.transpose();
        std::array<std::vector<size_t>, 2> connections;
        connections[0].reserve(std::max(overlap_sqrt.rows(), overlap_sqrt.cols()));
        connections[1].reserve(std::max(overlap_sqrt.rows(), overlap_sqrt.cols()));

        for (int k = 0; k < overlap_sqrt_transposed.outerSize(); ++k) { // cols
            double maxval = threshold_sqrt;
            size_t idx_from;
            size_t idx_to;

            for (eigen_iterator_double_t triple(overlap_sqrt_transposed, k); triple; ++triple) {
                if (triple.value() > maxval) {
                    idx_from = triple.col();
                    idx_to = triple.row();
                    maxval = triple.value();
                }
            }

            if (maxval > threshold_sqrt && rows_with_maxval[idx_to] == idx_from) {
                connections[0].push_back(idx_from);
                connections[1].push_back(idx_to);
            }
        }

        connections[0].shrink_to_fit();
        connections[1].shrink_to_fit();

        return connections;
    }

    ////////////////////////////////////////////////////////////////////
    /// Methods to build, transform, and destroy the system ////////////
    ////////////////////////////////////////////////////////////////////

    void buildHamiltonian() {
        // Build basis, also constructs the Hamiltonian matrix without interaction
        this->buildBasis();

        // Initialize Hamiltonian matrix with interaction if required
        if (is_new_hamiltonian_required) {
            if (is_interaction_already_contained) {

                // Check variables for consistency
                if (memory_saving || basisvectors_unperturbed_cache.size() == 0 ||
                    hamiltonian_unperturbed_cache.size() == 0) {
                    throw std::runtime_error("Inconsistent variables at " + std::string(__FILE__) +
                                             ":" + std::to_string(__LINE__) + ".");
                }

                // Reset the Hamiltonian if it already contains the interaction
                basisvectors = basisvectors_unperturbed_cache;
                hamiltonian = hamiltonian_unperturbed_cache;
            } else if (!memory_saving) {

                // Store the Hamiltonian without interaction
                basisvectors_unperturbed_cache = basisvectors;
                hamiltonian_unperturbed_cache = hamiltonian;
            }

            // Build interaction
            this->initializeInteraction();

            // Add interaction to the Hamiltonian
            this->addInteraction();

            if (memory_saving) {
                // Delete the variables used to add the interaction to the Hamiltonian
                this->deleteInteraction();
            }

            is_interaction_already_contained = true;
            is_new_hamiltonian_required = false;
        }
    }

    void buildInteraction() {
        // Build basis
        this->buildBasis();

        // Initialize interaction
        this->initializeInteraction(); // this method checks by itself whether a new initialization
                                       // is required
    }

    void buildBasis() {
        // Check variables for consistency
        if (((hamiltonian.size() == 0) != states.empty()) ||
            ((hamiltonian.size() == 0) != (basisvectors.size() == 0))) {
            throw std::runtime_error("Inconsistent variables at " + std::string(__FILE__) + ":" +
                                     std::to_string(__LINE__) + ".");
        }

        // In case of no new basis restrictions and already initialized basis, there is nothing to
        // do
        if (!states.empty() && states_to_add.empty() && range_n.empty() && range_l.empty() &&
            range_j.empty() && range_m.empty() &&
            energy_min == std::numeric_limits<double>::lowest() &&
            energy_max == std::numeric_limits<double>::max()) { // TODO check for new threshold, too
            return;
        }

        // Check whether the basis does not exist
        if (states.empty()) {

            // Initialize the basis
            this->initializeBasis();

            // Forget basis restrictions as they were applied now
            this->forgetRestrictions();

        } else {
            this->updateEverything();
        }

        // Check whether the basis is empty
        if (basisvectors.rows() == 0) {
            throw std::runtime_error("The basis contains no states.");
        }
        if (basisvectors.cols() == 0) {
            throw std::runtime_error("The basis contains no vectors.");
        }
    }

    void diagonalize(double energy_lower_bound, double energy_upper_bound) {
        this->diagonalize(energy_lower_bound, energy_upper_bound, 0);
    }

    void diagonalize(double energy_lower_bound, double energy_upper_bound, double threshold) {
#ifdef WITH_INTEL_MKL
        this->buildHamiltonian();

        // Check if already diagonal
        if (checkIsDiagonal(hamiltonian)) {
            return;
        }

        // Estimate number of found eigenvalues
        std::vector<double> diagonal_max(hamiltonian.outerSize(), 0);
        std::vector<double> diagonal_val;
        std::vector<int> diagonal_idx;
        diagonal_val.reserve(hamiltonian.outerSize());
        diagonal_idx.reserve(hamiltonian.outerSize());
        for (int k = 0; k < hamiltonian.outerSize(); ++k) {
            double val = 0;
            for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple(hamiltonian, k); triple;
                 ++triple) {
                if (triple.row() == triple.col()) {
                    val = std::real(triple.value());
                } else if (triple.row() != triple.col()) {
                    diagonal_max[k] = std::max(diagonal_max[k], std::abs(triple.value()));
                }
            }
            diagonal_idx.push_back(k);
            diagonal_val.push_back(val);
        }

        MKL_INT m0 =
            std::count_if(diagonal_val.begin(), diagonal_val.end(),
                          [&energy_lower_bound, &energy_upper_bound](const double &val) {
                              return (val < energy_upper_bound) && (val > energy_lower_bound);
                          }) +
            std::count_if(diagonal_idx.begin(), diagonal_idx.end(),
                          [&energy_lower_bound, &energy_upper_bound, &diagonal_val,
                           &diagonal_max](const int &i) {
                              return ((diagonal_val[i] >= energy_upper_bound) ||
                                      (diagonal_val[i] <= energy_lower_bound)) &&
                                  (diagonal_val[i] < energy_upper_bound + diagonal_max[i]) &&
                                  (diagonal_val[i] > energy_lower_bound - diagonal_max[i]);
                          });

        // Inplace conversion of the Hamiltonian to CSR with one-based indexing
        hamiltonian = hamiltonian.transpose();
        hamiltonian.makeCompressed();
        for (int i = 0; i < hamiltonian.rows() + 1; ++i) {
            hamiltonian.outerIndexPtr()[i] += 1;
        }
        for (int i = 0; i < hamiltonian.nonZeros(); ++i) {
            hamiltonian.innerIndexPtr()[i] += 1;
        }

        // Set default parameters for the diagonalization as described at
        // https://software.intel.com/en-us/mkl-developer-reference-c-extended-eigensolver-input-parameters.

        std::vector<MKL_INT> fpm(128);
        feastinit(&fpm[0]);
        fpm[0] = 1;  // enables terminal output
        fpm[1] = 6;  // number of contour points
        fpm[26] = 0; // disables matrix checker
        fpm[3] = 5;  // maximum number of refinement loops allowed
        if (threshold != 0) {
            // Adapt the error trace stopping criteria (10-fpm[2])
            fpm[2] = std::min(std::round(-std::log10(threshold)), 12.);
        }

        // Do the diagonalization
        {
            MKL_INT n = hamiltonian.rows(); // size of the matrix
            MKL_INT m;                      // will contain the number of eigenvalues
            std::vector<Scalar> x(m0 * n);  // the first m columns will contain the eigenvectors
            {
                std::vector<double> e(m0); // will contain the first m eigenvalues
                {
                    char uplo = 'F';             // full matrix is stored
                    MKL_INT info;                // will contain return codes
                    double epsout;               // will contain relative error
                    MKL_INT loop;                // will contain number of used refinement
                    std::vector<double> res(m0); // will contain the residual errors

                    this->feast_csrev(&uplo, &n, hamiltonian.valuePtr(),
                                      hamiltonian.outerIndexPtr(), hamiltonian.innerIndexPtr(),
                                      &fpm[0], &epsout, &loop, &energy_lower_bound,
                                      &energy_upper_bound, &m0, &e[0], &x[0], &m, &res[0], &info);
                    if (info != 0) {
                        throw std::runtime_error("Diagonalization with FEAST failed.");
                    }
                }

                // Build the new hamiltonian
                hamiltonian.resize(m, m);
                hamiltonian.setZero();
                hamiltonian.reserve(m);
                for (int idx = 0; idx < m; ++idx) {
                    hamiltonian.insert(idx, idx) = e[idx];
                }
                hamiltonian.makeCompressed();
            }

            // Transform the basis vectors
            Eigen::SparseMatrix<Scalar> evecs =
                Eigen::Map<Eigen::MatrixX<Scalar>>(&x[0], n, m).sparseView();
            if (threshold == 0) {
                basisvectors = basisvectors * evecs;
            } else {
                basisvectors = (basisvectors * evecs).pruned(threshold, 1);
            }
        }
#else  // WITH_INTEL_MKL
        (void)energy_lower_bound;
        (void)energy_upper_bound;
        (void)threshold;
        throw std::runtime_error(
            "The method does not work because the program was compiled without MKL support.");
#endif // WITH_INTEL_MKL
    }

    void diagonalize() { this->diagonalize(0); }

    void diagonalize(double threshold) {
        this->buildHamiltonian();

        // Check if already diagonal
        if (checkIsDiagonal(hamiltonian)) {
            return;
        }

#if defined EIGEN_USE_LAPACKE || WITH_INTEL_MKL

        // Diagonalize hamiltonian
        char jobz = 'V';                          // eigenvalues and eigenvectors are computed
        char uplo = 'U';                          // full matrix is stored, upper is used
        int n = hamiltonian.cols();               // size of the matrix
        Eigen::MatrixX<Scalar> mat = hamiltonian; // matrix
        int lda = mat.outerStride();              // leading dimension
        eigen_vector_double_t evals(n);           // eigenvalues
        int info = LAPACKE_evd(LAPACK_COL_MAJOR, jobz, uplo, n, mat.data(), lda, evals.data());
        if (info != 0) {
            throw std::runtime_error("Diagonalization with LAPACKE failed.");
        }

        // Get eigenvectors
        Eigen::SparseMatrix<Scalar> evecs = mat.sparseView();

#else // EIGEN_USE_LAPACKE || WITH_INTEL_MKL

        // Diagonalize hamiltonian
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixX<Scalar>> eigensolver(hamiltonian);

        // Get eigenvalues and eigenvectors
        eigen_vector_double_t evals = eigensolver.eigenvalues();
        Eigen::SparseMatrix<Scalar> evecs = eigensolver.eigenvectors().sparseView();

#endif // EIGEN_USE_LAPACKE || WITH_INTEL_MKL

        // Build the new hamiltonian
        hamiltonian.setZero();
        hamiltonian.reserve(evals.size());
        for (int idx = 0; idx < evals.size(); ++idx) {
            hamiltonian.insert(idx, idx) = evals.coeffRef(idx);
        }
        hamiltonian.makeCompressed();

        // Transform the basis vectors
        if (threshold == 0) {
            basisvectors = basisvectors * evecs;
        } else {
            basisvectors = (basisvectors * evecs).pruned(threshold, 1);
        }

        // TODO call transformInteraction (see applyRightsideTransformator), perhaps not?
    }

    void canonicalize() {
        this->buildHamiltonian();

        // Transform the hamiltonian
        hamiltonian = basisvectors * hamiltonian * basisvectors.adjoint();

        // Transform the basis vectors
        basisvectors = basisvectors * basisvectors.adjoint();

        // TODO call transformInteraction (see applyRightsideTransformator), perhaps yes?
    }

    void unitarize() {
        this->buildHamiltonian();

        // Determine hash of the list of states
        size_t hashvalue_states =
            utils::hash_range(states.template get<0>().begin(), states.template get<0>().end());

        // Loop over basisvectors
        size_t num_basisvectors = basisvectors.outerSize();

        states.clear();
        size_t idx_new = 0;
        for (int k = 0; k < basisvectors.outerSize(); ++k) { // col == idx_vector
            size_t hashvalue = hashvalue_states;
            for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple(basisvectors, k);
                 triple; ++triple) {
                utils::hash_combine(hashvalue, triple.row());
                utils::hash_combine(hashvalue, triple.value());
            }

            // Rewrite basisvectors as states
            std::stringstream ss;
            ss << std::hex << hashvalue;
            states.push_back(
                enumerated_state<State>(idx_new++, this->createStateFromLabel<State>(ss.str())));
        }
        states.shrink_to_fit();

        if (num_basisvectors != states.size()) {
            throw std::runtime_error("A hash collision occurred.");
        }

        // Adapt basisvectors
        basisvectors.resize(states.size(), states.size());
        basisvectors.setZero();
        basisvectors.reserve(states.size());
        for (int idx = 0; idx < states.size(); ++idx) {
            basisvectors.insert(idx, idx) = 1;
        }
        basisvectors.makeCompressed();

        // Delete caches as they are no longer meaningful
        basisvectors_unperturbed_cache.resize(0, 0);
        hamiltonian_unperturbed_cache.resize(0, 0);
    }

    void rotate(std::array<double, 3> to_z_axis, std::array<double, 3> to_y_axis) {
        auto euler_zyz = this->getEulerAngles(to_z_axis, to_y_axis);
        this->rotate(euler_zyz[0], euler_zyz[1], euler_zyz[2]);
    }

    void rotate(double alpha, double beta, double gamma) { // TODO applyRightsideTransformator ?
        // Build Hamiltonian and basis
        this->buildHamiltonian();

        // Get the rotator for the basis states
        Eigen::SparseMatrix<Scalar> transformator = this->buildStaterotator(alpha, beta, gamma);

        // Rotate basis
        this->transformInteraction(basisvectors.adjoint());

        basisvectors = transformator * basisvectors;
        if (basisvectors_unperturbed_cache.size() != 0) {
            basisvectors_unperturbed_cache = transformator * basisvectors_unperturbed_cache;
        }

        this->transformInteraction(basisvectors);
    }

    void add(SystemBase<Scalar, State> &system) {
        // --- Build bases ---
        this->buildBasis();
        system.buildBasis();

        size_t size_this = basisvectors.cols();
        size_t size_other = system.basisvectors.cols();

        // --- Combine system specific variables ---
        this->incorporate(system);

        // --- Combine universal variables ---
        if (memory_saving != system.memory_saving) {
            throw std::runtime_error(
                "The value of the variable 'memory_saving' must be the same for both systems.");
        }
        if (is_interaction_already_contained != system.is_interaction_already_contained) {
            throw std::runtime_error("The value of the variable 'is_interaction_already_contained' "
                                     "must be the same for both systems.");
        }
        if (is_new_hamiltonian_required != system.is_new_hamiltonian_required) {
            throw std::runtime_error(
                "The value of the variable 'is_new_hamiltonian_required' must be the same "
                "for both systems.");
        }

        // --- Combine states and build transformators ---
        Eigen::SparseMatrix<Scalar> transformator;
        {
            std::vector<Eigen::Triplet<Scalar>> transformator_triplets;
            transformator_triplets.reserve(system.states.size());

            for (const auto &entry : system.states) {
                auto state_iter = states.template get<1>().find(entry.state);

                size_t newidx = states.size();
                if (state_iter == states.template get<1>().end()) {
                    states.push_back(enumerated_state<State>(newidx, entry.state));
                } else {
                    newidx = state_iter->idx;
                }

                transformator_triplets.emplace_back(newidx, entry.idx, 1);
            }

            transformator.resize(states.size(), system.basisvectors.rows());
            transformator.setFromTriplets(transformator_triplets.begin(),
                                          transformator_triplets.end());
        }

        // --- Combine basisvectors and basisvectors_unperturbed_cache ---
        basisvectors.conservativeResize(states.size(), size_this + size_other);
        basisvectors.rightCols(size_other) = transformator * system.basisvectors;

        if ((basisvectors_unperturbed_cache.size() != 0) !=
            (system.basisvectors_unperturbed_cache.size() != 0)) {
            throw std::runtime_error("Inconsistent variables at " + std::string(__FILE__) + ":" +
                                     std::to_string(__LINE__) + ".");
        }

        if (basisvectors_unperturbed_cache.size() != 0) {
            basisvectors_unperturbed_cache.conservativeResize(
                states.size(),
                basisvectors_unperturbed_cache.cols() +
                    system.basisvectors_unperturbed_cache.cols());
            basisvectors_unperturbed_cache.rightCols(system.basisvectors_unperturbed_cache.cols()) =
                transformator * system.basisvectors_unperturbed_cache;
        }

        // --- Check if basis vectors are orthogonal ---
        Eigen::SparseMatrix<Scalar> tmp =
            (basisvectors.leftCols(size_this).adjoint() * basisvectors.rightCols(size_other))
                .pruned(1e-12, 1);
        if (tmp.nonZeros() != 0) {
            throw std::runtime_error(
                "Two systems cannot be combined if their basis vectors are not orthogonal.");
        }

        // --- Combine hamiltonian and hamiltonian_unperturbed_cache ---
        Eigen::SparseMatrix<Scalar> shifter(hamiltonian.rows() + system.hamiltonian.rows(),
                                            system.hamiltonian.rows());
        {
            std::vector<Eigen::Triplet<Scalar>> shifter_triplets;
            shifter_triplets.reserve(system.hamiltonian.rows());

            for (size_t idx = 0; idx < system.hamiltonian.rows(); ++idx) {
                shifter_triplets.emplace_back(hamiltonian.rows() + idx, idx, 1);
            }

            shifter.setFromTriplets(shifter_triplets.begin(), shifter_triplets.end());
        }

        hamiltonian.conservativeResize(hamiltonian.rows() + system.hamiltonian.rows(),
                                       hamiltonian.cols() + system.hamiltonian.cols());
        hamiltonian.rightCols(system.hamiltonian.cols()) = shifter * system.hamiltonian;

        if ((hamiltonian_unperturbed_cache.size() != 0) !=
            (system.hamiltonian_unperturbed_cache.size() != 0)) {
            throw std::runtime_error("Inconsistent variables at " + std::string(__FILE__) + ":" +
                                     std::to_string(__LINE__) + ".");
        }

        if (hamiltonian_unperturbed_cache.size() != 0) {
            hamiltonian_unperturbed_cache.conservativeResize(
                hamiltonian_unperturbed_cache.rows() + system.hamiltonian_unperturbed_cache.rows(),
                hamiltonian_unperturbed_cache.cols() + system.hamiltonian_unperturbed_cache.cols());
            hamiltonian_unperturbed_cache.rightCols(system.hamiltonian_unperturbed_cache.cols()) =
                shifter * system.hamiltonian_unperturbed_cache;
        }
    }

    void constrainBasisvectors(std::vector<size_t> indices_of_wanted_basisvectors) {
        this->buildHamiltonian();

        // Check if indices are unique
        {
            std::set<size_t> tmp(indices_of_wanted_basisvectors.begin(),
                                 indices_of_wanted_basisvectors.end());
            if (tmp.size() < indices_of_wanted_basisvectors.size()) {
                throw std::runtime_error("Indices are occuring multiple times.");
            }
        }

        // Build transformator and remove vectors
        std::vector<Eigen::Triplet<Scalar>> triplets_transformator;
        triplets_transformator.reserve(indices_of_wanted_basisvectors.size());

        size_t idx_new = 0;
        for (const auto &idx : indices_of_wanted_basisvectors) {
            if (idx >= basisvectors.cols()) {
                throw std::runtime_error("A basis vector with index " + std::to_string(idx) +
                                         " could not be found.");
            }
            triplets_transformator.emplace_back(idx, idx_new++, 1);
        }

        this->applyRightsideTransformator(triplets_transformator);
    }

    void applySchriefferWolffTransformation(SystemBase<Scalar, State> &system0) {
        this->diagonalize();
        system0.buildHamiltonian();

        // Check that system, on which applySchriefferWolffTransformation() is called, is unitary
        if (!this->checkIsUnitary(basisvectors)) {
            throw std::runtime_error("The system, on which applySchriefferWolffTransformation() is "
                                     "called, is not unitary. Call unitarize() on the system.");
        }

        // Check that system0, i.e. the unperturbed system, is unitary
        if (!this->checkIsUnitary(system0.basisvectors)) {
            throw std::runtime_error(
                "The unperturbed system passed to applySchriefferWolffTransformation() is not "
                "unitary. Call unitarize() on the system.");
        }

        // Check that the unperturbed system is diagonal
        if (!this->checkIsDiagonal(system0.hamiltonian)) {
            throw std::runtime_error("The unperturbed system passed to "
                                     "applySchriefferWolffTransformation() is not diagonal.");
        }

        // --- Express the basis vectors of system0 as a linearcombination of all states ---

        // Combine states and build transformators
        Eigen::SparseMatrix<Scalar> transformator;
        {
            std::vector<Eigen::Triplet<Scalar>> transformator_triplets;
            transformator_triplets.reserve(system0.states.size());

            for (const auto &entry : system0.states) {
                auto state_iter = states.template get<1>().find(entry.state);

                // Check that all the states of system0 occur (since we already checked that the
                // number of states is the same, this ensures that all the states are the same)
                if (state_iter == states.template get<1>().end()) {
                    throw std::runtime_error(
                        "The unperturbed system contains states that are not occuring.");
                }
                size_t newidx = state_iter->idx;

                transformator_triplets.emplace_back(newidx, entry.idx, 1);
            }

            transformator.resize(states.size(), system0.states.size());
            transformator.setFromTriplets(transformator_triplets.begin(),
                                          transformator_triplets.end());
        }

        // Get basisvectors of system0
        Eigen::SparseMatrix<Scalar> low_energy_basis0 = transformator * system0.basisvectors;

        // --- Find basisvectors which corresponds to the basis vectors of system0 ---

        // Get overlaps between basis vectors
        eigen_vector_double_t overlap = (basisvectors.adjoint() * low_energy_basis0).cwiseAbs2() *
            eigen_vector_double_t::Ones(low_energy_basis0.cols());

        // Get indices of the low_energy_basis0.cols() largest entries and build transformator
        {
            std::vector<int> indices(basisvectors.cols());
            std::iota(indices.begin(), indices.end(), 0);
            std::nth_element(indices.begin(), indices.begin() + low_energy_basis0.cols(),
                             indices.end(),
                             [&overlap](int a, int b) { return overlap[a] > overlap[b]; });

            std::vector<Eigen::Triplet<Scalar>> transformator_triplets;
            transformator_triplets.reserve(low_energy_basis0.cols());

            for (int i = 0; i < low_energy_basis0.cols(); ++i) {
                transformator_triplets.emplace_back(indices[i], i, 1);
            }

            transformator.resize(basisvectors.cols(), low_energy_basis0.cols());
            transformator.setFromTriplets(transformator_triplets.begin(),
                                          transformator_triplets.end());
        }

        // Get corresponding basis vectors
        Eigen::SparseMatrix<Scalar> low_energy_basis = basisvectors * transformator;

        // --- Perform the Schrieffer Wolff transformation ---

        // Projector on the selected basis vectors (typically low-energy subspace)
        Eigen::SparseMatrix<Scalar> projector0 = low_energy_basis0 * low_energy_basis0.adjoint();
        Eigen::SparseMatrix<Scalar> projector = low_energy_basis * low_energy_basis.adjoint();

        // Reflection operator which change signes of the selectes basis vectors but act
        // trivially on others
        Eigen::MatrixX<Scalar> reflection0 =
            Eigen::MatrixX<Scalar>::Identity(states.size(), states.size()) - 2 * projector0;
        Eigen::MatrixX<Scalar> reflection =
            Eigen::MatrixX<Scalar>::Identity(states.size(), states.size()) - 2 * projector;

        // Direct rotator from low_energy_basis to low_energy_basis0
        Eigen::SparseMatrix<Scalar> rotator = (reflection0 * reflection).sqrt().sparseView();

        if (std::isnan(std::abs(rotator.coeffRef(0, 0)))) {
            throw std::runtime_error(
                "Error in calculating the matrix square root. Try to use the picomplex module or "
                "increase accuracy by choosing smaller thresholds.");
        }

        // Calculate the effective Hamiltonian
        transformator = basisvectors.adjoint() * rotator.adjoint() * projector0 * low_energy_basis0;

        this->applyRightsideTransformator(transformator);
    }

    ////////////////////////////////////////////////////////////////////
    /// Methods to manipulate entries of the Hamiltonian ///////////////
    ////////////////////////////////////////////////////////////////////

    size_t getStateIndex(const State &searched_state) {
        this->buildBasis();

        if (utils::is_true(searched_state.isGeneralized())) {
            throw std::runtime_error("The method must not be called on a generalized state.");
        }

        auto state_iter = states.template get<1>().find(searched_state);
        if (state_iter == states.template get<1>().end()) {
            throw std::runtime_error("The method is called on a non-existing state.");
        }

        return state_iter->idx;
    }

    std::vector<size_t> getStateIndex(const std::vector<State> &searched_states) {
        this->buildBasis();

        std::vector<size_t> indices;
        indices.reserve(searched_states.size());

        for (const auto &state : searched_states) {
            if (utils::is_true(state.isGeneralized())) {
                throw std::runtime_error("The method must not be called on a generalized state.");
            }

            auto state_iter = states.template get<1>().find(state);
            if (state_iter == states.template get<1>().end()) {
                throw std::runtime_error("The method is called on a non-existing state.");
            }
            indices.push_back(state_iter->idx);
        }

        return indices;
    }

    size_t getBasisvectorIndex(const State &searched_state) {
        this->buildBasis();

        size_t stateidx = this->getStateIndex(searched_state);

        double maxval = -1;
        size_t col_with_maxval = -1;
        for (int k = 0; k < basisvectors.outerSize(); ++k) { // col == idx_vector
            for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple(basisvectors, k);
                 triple; ++triple) {
                if (triple.row() == stateidx) {
                    if (std::abs(triple.value()) > maxval) {
                        col_with_maxval = triple.col();
                        maxval = std::abs(triple.value());
                    }
                    break;
                }
            }
        }

        return col_with_maxval;
    }

    std::vector<size_t> getBasisvectorIndex(const std::vector<State> &searched_states) {
        this->buildBasis();

        // Ensure that the states within searched_states are unique
        {
            std::set<State> tmp(searched_states.begin(), searched_states.end());
            if (tmp.size() < searched_states.size()) {
                throw std::runtime_error("States are occuring multiple times.");
            }
        }

        // Construct canonical basis vectors
        Eigen::SparseMatrix<Scalar> canonicalbasis;
        {
            std::vector<size_t> state_indices = this->getStateIndex(searched_states);
            std::vector<Eigen::Triplet<Scalar>> canonicalbasis_triplets;
            canonicalbasis_triplets.reserve(searched_states.size());
            for (size_t i = 0; i < state_indices.size(); ++i) {
                canonicalbasis_triplets.emplace_back(state_indices[i], i, 1);
            }
            canonicalbasis.resize(states.size(), searched_states.size());
            canonicalbasis.setFromTriplets(canonicalbasis_triplets.begin(),
                                           canonicalbasis_triplets.end());
        }

        // Get overlaps between basis vectors
        eigen_sparse_double_t overlap = (canonicalbasis.adjoint() * basisvectors).cwiseAbs2();
        eigen_vector_double_t overlap_total =
            eigen_vector_double_t::Ones(canonicalbasis.cols()).transpose() * overlap;

        // Get indices of the canonicalbasis.cols() largest entries
        std::vector<size_t> indices_available(basisvectors.cols());
        std::iota(indices_available.begin(), indices_available.end(), 0);
        std::nth_element(
            indices_available.begin(), indices_available.begin() + canonicalbasis.cols(),
            indices_available.end(),
            [&overlap_total](size_t a, size_t b) { return overlap_total[a] > overlap_total[b]; });
        indices_available.resize(canonicalbasis.cols());

        // Resort the indices to match the order of searched_states, TODO use Hungarian algorithm
        std::vector<size_t> indices(canonicalbasis.cols(), std::numeric_limits<size_t>::max());
        for (const auto &k : indices_available) {
            size_t row_with_maxval;
            double maxval = -1;
            for (eigen_iterator_double_t triple(overlap, k); triple; ++triple) {
                if (indices[triple.row()] == std::numeric_limits<size_t>::max() &&
                    triple.value() > maxval) {
                    row_with_maxval = triple.row();
                    maxval = triple.value();
                }
            }

            if (maxval == -1) {
                throw std::runtime_error(
                    "There is a state for which no basis vector could be found.");
            }

            indices[row_with_maxval] = k;
        }

        return indices;
    }

    void forgetStatemixing() {
        this->diagonalize();

        std::vector<Eigen::Triplet<Scalar>> basisvectors_triplets;
        basisvectors_triplets.reserve(basisvectors.cols());
        double threshold = std::sqrt(0.5);

        for (int k = 0; k < basisvectors.outerSize(); ++k) { // col == idx_vector
            for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple(basisvectors, k);
                 triple; ++triple) {
                if (std::abs(triple.value()) > threshold) {
                    basisvectors_triplets.emplace_back(triple.row(), triple.col(), 1);
                    break;
                }
            }
        }

        if (basisvectors_triplets.size() < basisvectors.cols()) {
            throw std::runtime_error(
                "The states are mixed too strongly for calling forgetStatemixing().");
        }

        basisvectors.setFromTriplets(basisvectors_triplets.begin(), basisvectors_triplets.end());
    }

    Scalar getHamiltonianEntry(const State &state_row, const State &state_col) {
        this->buildHamiltonian();

        size_t idx_row = this->getStateIndex(state_row);
        size_t idx_col = this->getStateIndex(state_col);

        return (basisvectors.row(idx_row) * hamiltonian).dot(basisvectors.row(idx_col));
    }

    void setHamiltonianEntry(const State &state_row, const State &state_col, Scalar value) {
        this->buildHamiltonian();

        size_t idx_row = this->getStateIndex(state_row);
        size_t idx_col = this->getStateIndex(state_col);

        value -= (basisvectors.row(idx_row) * hamiltonian).dot(basisvectors.row(idx_col));

        Eigen::SparseMatrix<Scalar> tmp(states.size(), states.size());
        tmp.reserve(2);
        tmp.insert(idx_row, idx_col) = value;
        if (idx_row != idx_col) {
            tmp.insert(idx_col, idx_row) = utils::conjugate(value);
        }
        tmp.makeCompressed();

        hamiltonian += basisvectors.adjoint() * tmp * basisvectors;
    }

    void addHamiltonianEntry(const State &state_row, const State &state_col, Scalar value) {
        this->buildHamiltonian();

        size_t idx_row = this->getStateIndex(state_row);
        size_t idx_col = this->getStateIndex(state_col);

        Eigen::SparseMatrix<Scalar> tmp(states.size(), states.size());
        tmp.reserve(2);
        tmp.insert(idx_row, idx_col) = value;
        if (idx_row != idx_col) {
            tmp.insert(idx_col, idx_row) = utils::conjugate(value);
        }
        tmp.makeCompressed();

        hamiltonian += basisvectors.adjoint() * tmp * basisvectors;
    }

protected:
    SystemBase(MatrixElementCache &cache)
        : cache(cache), threshold_for_sqnorm(0.05),
          energy_min(std::numeric_limits<double>::lowest()),
          energy_max(std::numeric_limits<double>::max()), memory_saving(false),
          is_interaction_already_contained(false), is_new_hamiltonian_required(false) {}

    SystemBase(MatrixElementCache &cache, bool memory_saving)
        : cache(cache), threshold_for_sqnorm(0.05),
          energy_min(std::numeric_limits<double>::lowest()),
          energy_max(std::numeric_limits<double>::max()), memory_saving(memory_saving),
          is_interaction_already_contained(false), is_new_hamiltonian_required(false) {}

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

    MatrixElementCache &cache;

    double threshold_for_sqnorm;

    double energy_min, energy_max;
    std::set<int> range_n, range_l;
    std::set<float> range_j, range_m;
    std::set<State> states_to_add;

    bool memory_saving;
    bool is_interaction_already_contained;
    bool is_new_hamiltonian_required;

    typename states_set<State>::type states;
    Eigen::SparseMatrix<Scalar> basisvectors;
    Eigen::SparseMatrix<Scalar> hamiltonian;
    Eigen::SparseMatrix<Scalar> basisvectors_unperturbed_cache;
    Eigen::SparseMatrix<Scalar> hamiltonian_unperturbed_cache;

    ////////////////////////////////////////////////////////////////////
    /// Helper method that shoul be called by the derived classes //////
    ////////////////////////////////////////////////////////////////////

    void onParameterChange() {
        // Check variables for consistency
        if ((basisvectors_unperturbed_cache.size() == 0) !=
            (hamiltonian_unperturbed_cache.size() == 0)) {
            throw std::runtime_error("Inconsistent variables at " + std::string(__FILE__) + ":" +
                                     std::to_string(__LINE__) + ".");
        }

        // Throw error if the Hamiltonian cannot be changed anymore
        if (is_interaction_already_contained && basisvectors_unperturbed_cache.size() == 0) {
            throw std::runtime_error(
                "If memory saving is activated or unitarize() has been called, one cannot change "
                "parameters after interaction was added to the Hamiltonian.");
        }

        is_new_hamiltonian_required = true;
    }

    void onSymmetryChange() {
        // Throw error if the Symmetry cannot be changed anymore
        if (!states.empty()) {
            throw std::runtime_error("One cannot change symmetries after the basis was built.");
        }
    }

    ////////////////////////////////////////////////////////////////////
    /// Helper methods to check diagonality and unitarity of a matrix //
    ////////////////////////////////////////////////////////////////////

    bool checkIsDiagonal(const Eigen::SparseMatrix<Scalar> &mat) {
        Eigen::SparseMatrix<Scalar> tmp = mat;
        tmp.prune(1e-12, 1);

        for (int k = 0; k < tmp.outerSize(); ++k) {
            for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple(tmp, k); triple;
                 ++triple) {
                if (triple.row() != triple.col()) {
                    return false;
                }
            }
        }
        return true;
    }

    bool checkIsUnitary(const Eigen::SparseMatrix<Scalar> &mat) {
        Eigen::SparseMatrix<Scalar> tmp = (mat.adjoint() * mat).pruned(1e-12, 1);

        if (tmp.nonZeros() != tmp.outerSize()) {
            return false;
        }

        for (int k = 0; k < tmp.outerSize(); ++k) {
            for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple(tmp, k); triple;
                 ++triple) {
                if (triple.row() != triple.col()) {
                    return false;
                }
                if (std::abs(std::real(triple.value()) - 1) > 1e-12 ||
                    std::abs(std::imag(triple.value())) > 1e-12) {
                    return false;
                }
            }
        }
        return true;
    }

    ////////////////////////////////////////////////////////////////////
    /// Helper methods to check the validity of states /////////////////
    ////////////////////////////////////////////////////////////////////

    template <class V>
    bool checkIsQuantumnumberValid(V q, std::set<V> range_q) {
        return range_q.empty() || range_q.find(q) != range_q.end();
    }

    bool checkIsQuantumstateValid(const State &state) {
        return checkIsQuantumstateValid(state, state.isArtificial());
    }

    bool checkIsQuantumstateValid(const State &state, bool a) {
        return a ||
            (checkIsQuantumnumberValid(state.getN(), range_n) &&
             checkIsQuantumnumberValid(state.getL(), range_l) &&
             checkIsQuantumnumberValid(state.getJ(), range_j) &&
             checkIsQuantumnumberValid(state.getM(), range_m));
    }

    bool checkIsQuantumstateValid(const State &state, std::array<bool, 2> a) {
        bool valid = true;
        for (int idx = 0; idx < 2; ++idx) {
            valid = valid &&
                (a[idx] ||
                 (checkIsQuantumnumberValid(state.getN(idx), range_n) &&
                  checkIsQuantumnumberValid(state.getL(idx), range_l) &&
                  checkIsQuantumnumberValid(state.getJ(idx), range_j) &&
                  checkIsQuantumnumberValid(state.getM(idx), range_m)));
        }
        return valid;
    }

    bool checkIsEnergyValid(double e) {
        return (e > energy_min || energy_min == std::numeric_limits<double_t>::lowest()) &&
            (e < energy_max ||
             energy_max ==
                 std::numeric_limits<double_t>::max()); // TODO take into account numerical errors
    }

    ////////////////////////////////////////////////////////////////////
    /// Helper methods to change the set of basis vectors //////////////
    ////////////////////////////////////////////////////////////////////

    void applyLeftsideTransformator(std::vector<Eigen::Triplet<Scalar>> &triplets_transformator) {
        Eigen::SparseMatrix<Scalar> transformator(triplets_transformator.size(),
                                                  basisvectors.rows());
        transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

        this->applyLeftsideTransformator(transformator);
    }

    void applyLeftsideTransformator(Eigen::SparseMatrix<Scalar> &transformator) {
        // Apply transformator in order to remove rows from the basisvector matrix (i.e. states)
        basisvectors = transformator * basisvectors;
        if (basisvectors_unperturbed_cache.size() != 0) {
            basisvectors_unperturbed_cache = transformator * basisvectors_unperturbed_cache;
        }
    }

    void applyRightsideTransformator(std::vector<Eigen::Triplet<Scalar>> &triplets_transformator) {
        Eigen::SparseMatrix<Scalar> transformator(basisvectors.cols(),
                                                  triplets_transformator.size());
        transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

        this->applyRightsideTransformator(transformator);
    }

    void applyRightsideTransformator(Eigen::SparseMatrix<Scalar> &transformator) {
        // Apply transformator in order to remove columns from the basisvector matrix (i.e. basis
        // vectors)
        basisvectors = basisvectors * transformator;
        if (basisvectors_unperturbed_cache.size() != 0) {
            basisvectors_unperturbed_cache = basisvectors_unperturbed_cache * transformator;
        }

        // Apply transformator in order to remove rows and columns from the matrices that help
        // constructing the total Hamiltonian
        this->transformInteraction(transformator);

        // Apply transformator in order to remove rows and columns from the total Hamiltonian
        hamiltonian = transformator.adjoint() * hamiltonian * transformator;
        if (hamiltonian_unperturbed_cache.size() != 0) {
            hamiltonian_unperturbed_cache =
                transformator.adjoint() * hamiltonian_unperturbed_cache * transformator;
        }
    }

    template <typename F>
    void removeRestrictedStates(F &&checkIsValidEntry) {
        // Validate the call signature of the passed in object
        static_assert(std::is_same<decltype(checkIsValidEntry(
                                       std::declval<enumerated_state<State> const &>())),
                                   bool>::value,
                      "Function signature mismatch!  Expected `bool "
                      "checkIsValidEntry(enumerated_state<State> const &)'");

        // Build transformator and remove states
        typename states_set<State>::type states_new;
        states_new.reserve(states.size());
        std::vector<Eigen::Triplet<Scalar>> triplets_transformator;
        triplets_transformator.reserve(states.size());

        size_t idx_new = 0;

        for (const auto &entry : states) {
            if (checkIsValidEntry(entry)) {
                states_new.push_back(enumerated_state<State>(idx_new, entry.state));
                triplets_transformator.emplace_back(idx_new, entry.idx, 1);
                ++idx_new;
            }
        }

        states_new.shrink_to_fit();
        states = states_new;

        this->applyLeftsideTransformator(triplets_transformator);
    }

    ////////////////////////////////////////////////////////////////////
    /// Utility methods ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    template <class S>
    S createStateFromLabel(const std::string &label) const;

    template <class V>
    void range(std::set<V> &rset, V rmin, V rmax) {
        rset.clear();
        for (V r = rmin; r <= rmax; ++r) {
            rset.insert(r);
        }
    }

    Eigen::Matrix<double, 3, 3> buildRotator(std::array<double, 3> to_z_axis,
                                             std::array<double, 3> to_y_axis) {
        auto to_z_axis_mapped = Eigen::Map<Eigen::Matrix<double, 3, 1>>(&to_z_axis[0]).normalized();
        auto to_y_axis_mapped = Eigen::Map<Eigen::Matrix<double, 3, 1>>(&to_y_axis[0]).normalized();

        double tolerance = 1e-16;
        if (std::abs(to_z_axis_mapped.dot(to_y_axis_mapped)) > tolerance) {
            throw std::runtime_error("The z-axis and the y-axis are not orhogonal.");
        }

        Eigen::Matrix<double, 3, 3> transformator;
        transformator << to_y_axis_mapped.cross(to_z_axis_mapped), to_y_axis_mapped,
            to_z_axis_mapped;

        return transformator;
    }

    Eigen::Matrix<double, 3, 3> buildRotator(const double &alpha, const double &beta,
                                             const double &gamma) {
        Eigen::Matrix<double, 3, 3> transformator;

        transformator = Eigen::AngleAxisd(alpha, Eigen::Matrix<double, 3, 1>::UnitZ()) *
            Eigen::AngleAxisd(beta, Eigen::Matrix<double, 3, 1>::UnitY()) *
            Eigen::AngleAxisd(gamma, Eigen::Matrix<double, 3, 1>::UnitZ());

        return transformator;
    }

    Eigen::Matrix<double, 3, 1> getEulerAngles(const std::array<double, 3> &to_z_axis,
                                               const std::array<double, 3> &to_y_axis) {
        Eigen::Matrix<double, 3, 3> rotator = this->buildRotator(to_z_axis, to_y_axis);
        Eigen::Matrix<double, 3, 1> euler_zyz = rotator.eulerAngles(2, 1, 2);
        return euler_zyz; // alpha, beta, gamma
    }

private:
    void forgetRestrictions() {
        energy_min = std::numeric_limits<double>::lowest();
        energy_max = std::numeric_limits<double>::max();
        range_n.clear();
        range_l.clear();
        range_j.clear();
        range_m.clear();
        states_to_add.clear();
    }

#ifdef WITH_INTEL_MKL
    void feast_csrev(const char *uplo, const MKL_INT *n, const double *a, const MKL_INT *ia,
                     const MKL_INT *ja, MKL_INT *fpm, double *epsout, MKL_INT *loop,
                     const double *emin, const double *emax, MKL_INT *m0, double *e, double *x,
                     MKL_INT *m, double *res, MKL_INT *info) {
        dfeast_scsrev(uplo, n, a, ia, ja, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info);
    }

    void feast_csrev(const char *uplo, const MKL_INT *n, const MKL_Complex16 *a, const MKL_INT *ia,
                     const MKL_INT *ja, MKL_INT *fpm, double *epsout, MKL_INT *loop,
                     const double *emin, const double *emax, MKL_INT *m0, double *e,
                     MKL_Complex16 *x, MKL_INT *m, double *res, MKL_INT *info) {
        zfeast_hcsrev(uplo, n, a, ia, ja, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info);
    }
#endif // WITH_INTEL_MKL

#if defined EIGEN_USE_LAPACKE || WITH_INTEL_MKL
    int LAPACKE_evd(const int matrix_layout, const char jobz, const char uplo, const lapack_int n,
                    double *a, const lapack_int lda, double *w) {
        return LAPACKE_dsyevd(matrix_layout, jobz, uplo, n, a, lda, w);
    }

    int LAPACKE_evd(const int matrix_layout, const char jobz, const char uplo, const lapack_int n,
                    lapack_complex_double *a, const lapack_int lda, double *w) {
        return LAPACKE_zheevd(matrix_layout, jobz, uplo, n, a, lda, w);
    }
#endif // EIGEN_USE_LAPACKE || WITH_INTEL_MKL

    ////////////////////////////////////////////////////////////////////
    /// Method to update the system ////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    void updateEverything() {

        if (!range_n.empty() || !range_l.empty() || !range_j.empty() || !range_m.empty()) {

            ////////////////////////////////////////////////////////////////////
            /// Remove restricted states ///////////////////////////////////////
            ////////////////////////////////////////////////////////////////////

            // Remove states if the quantum numbers are not allowed
            removeRestrictedStates([=](const enumerated_state<State> &entry) -> bool {
                return this->checkIsQuantumstateValid(entry.state);
            });

            // Comunicate that the list of states has changed
            this->onStatesChange();
        }

        if (!range_n.empty() || !range_l.empty() || !range_j.empty() || !range_m.empty() ||
            energy_min != std::numeric_limits<double>::lowest() ||
            energy_max != std::numeric_limits<double>::max()) {

            ////////////////////////////////////////////////////////////////////
            /// Remove vectors with restricted energies or too small norm //////
            ////////////////////////////////////////////////////////////////////

            // Build transformator and remove vectors (if their energy is not allowed or the squared
            // norm too small)
            std::vector<Eigen::Triplet<Scalar>> triplets_transformator;
            triplets_transformator.reserve(basisvectors.cols());

            size_t idx_new = 0;
            for (int idx = 0; idx < basisvectors.cols(); ++idx) { // idx = col = num basis vector

                if (checkIsEnergyValid(std::real(hamiltonian.coeff(idx, idx)))) {
                    double_t sqnorm = 0;

                    // Calculate the square norm of the columns of the basisvector matrix
                    for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple(basisvectors,
                                                                                    idx);
                         triple; ++triple) {
                        sqnorm += std::pow(std::abs(triple.value()), 2);
                    }
                    if (sqnorm > threshold_for_sqnorm) {
                        triplets_transformator.emplace_back(idx, idx_new++, 1);
                    }
                }
            }

            this->applyRightsideTransformator(triplets_transformator);

            ////////////////////////////////////////////////////////////////////
            /// Remove states that barely occur within the vectors /////////////
            ////////////////////////////////////////////////////////////////////

            // Calculate the square norm of the rows of the basisvector matrix
            std::vector<double> sqnorm_list(basisvectors.rows(), 0);
            for (int k = 0; k < this->basisvectors.cols(); ++k) {
                for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple(basisvectors, k);
                     triple; ++triple) {
                    sqnorm_list[triple.row()] += std::pow(std::abs(triple.value()), 2);
                }
            }

            // Remove states if the squared norm is to small
            removeRestrictedStates([=](const enumerated_state<State> &entry) -> bool {
                return sqnorm_list[entry.idx] > threshold_for_sqnorm;
            });

            // Comunicate that the list of states has changed
            this->onStatesChange();
        }

        if (!states_to_add.empty()) {
            throw std::runtime_error(
                "States cannot be added anymore once the basis vectors have been created.");
        }

        // Forget basis restrictions as they were applied now
        this->forgetRestrictions();
    }

    ////////////////////////////////////////////////////////////////////
    /// Method for serialization ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int /*version*/) {
        ar &cache &threshold_for_sqnorm;
        ar &energy_min &energy_max &range_n &range_l &range_j &range_m &states_to_add;
        ar &memory_saving &is_interaction_already_contained &is_new_hamiltonian_required;
        ar &states &basisvectors &hamiltonian;
        ar &basisvectors_unperturbed_cache &hamiltonian_unperturbed_cache;
    }
};

#endif
