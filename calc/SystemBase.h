#ifndef SYSTEMBASE_H
#define SYSTEMBASE_H

#include "pyutils.h"
#include "dtypes.h"
#include "State.h"

#include <vector>
#include <numeric>
#include <string>
#include <set>
#include <limits>
#include <unordered_set>
#include <memory>
#include <string>
#include <stdexcept>
#include <complex>
#include <functional>
#include <boost/filesystem.hpp>
#include <iterator>
#include <algorithm>
#include <string>

template<class T> class SystemBase {
public:
    SystemBase(std::vector<T> states, std::string cachedir) : // TODO
        cachedir(boost::filesystem::absolute(cachedir)), energy_min(std::numeric_limits<double>::lowest()), energy_max(std::numeric_limits<double>::max()), states(states), memory_saving(false), is_interaction_already_contained(false), is_new_hamiltonianmatrix_required(false) {
        throw std::runtime_error( "Method yet not implemented." );
    }

    virtual ~SystemBase() = default;

    // TODO setThresholdForSqnorm()

    void setArtificialstates(const std::vector<T> &s) { // TODO [dummystates]
        states_artifical = s;
    }

    ////////////////////////////////////////////////////////////////////
    /// Methods to restrict the number of states inside the basis //////
    ////////////////////////////////////////////////////////////////////

    void restrictEnergy(double e_min, double e_max) { // restricts diagonal entries
        energy_min = e_min;
        energy_max = e_max;
    }

    void restrictN(int n_min, int n_max) {
        n_min = std::max(n_min, 1);
        this->range(range_n, n_min, n_max);
    }

    void restrictN(std::set<int> n) {
        range_n = n;
    }

    void restrictL(int l_min, int l_max) {
        l_min = std::max(l_min, 0);
        this->range(range_l, l_min, l_max);
    }

    void restrictL(std::set<int> l) {
        range_l = l;
    }

    void restrictJ(float j_min, float j_max) {
        j_min = std::fmax(j_min, 0.5);
        this->range(range_j, j_min, j_max);
    }

    void restrictJ(std::set<float> j) {
        range_j = j;
    }

    void restrictM(float m_min, float m_max) {
        this->range(range_m, m_min, m_max);
    }

    void restrictM(std::set<float> m) {
        range_m = m;
    }

    ////////////////////////////////////////////////////////////////////
    /// Methods to get properties of the system ////////////////////////
    ////////////////////////////////////////////////////////////////////

    const std::vector<T>& getStates() {
        this->buildBasis();
        return states;
    }

    eigen_sparse_t& getCoefficients() {
        this->buildBasis();
        return coefficients;
    }

    numpy::array getDiagonal() {
        this->buildHamiltonian();
        eigen_vector_double_t diagonal = hamiltonianmatrix.diagonal().real();
        return numpy::copy(diagonal);
    }

    eigen_sparse_t& getHamiltonianmatrix() {
        this->buildHamiltonian();
        return hamiltonianmatrix;
    }

    size_t getNumVectors() {
        // Build basis
        this->buildBasis();

        // Check variables for consistency
        if ( (coefficients.outerSize() != coefficients.cols()) || (coefficients.outerSize() != hamiltonianmatrix.rows()) || (coefficients.outerSize() != hamiltonianmatrix.cols()) ) {
            throw std::runtime_error( "Inconsistent variables at " + std::string(__FILE__) + ":" + std::to_string(__LINE__) + ".");
        }

        return coefficients.cols();
    }

    size_t getNumStates() {
        // Build basis
        this->buildBasis();

        // Check variables for consistency
        if ( (coefficients.innerSize() != coefficients.rows()) || (static_cast<size_t>(coefficients.innerSize()) != states.size()) ) {
            throw std::runtime_error( "Inconsistent variables at " + std::string(__FILE__) + ":" + std::to_string(__LINE__) + ".");
        }

        return coefficients.rows();
    }

    ////////////////////////////////////////////////////////////////////
    /// Methods to build, transform, and destroy the system ////////////
    ////////////////////////////////////////////////////////////////////

    void buildHamiltonian() {
        // Build basis, also constructs the Hamiltonian matrix without interaction
        this->buildBasis();

        // Initialize Hamiltonian matrix with interaction if required
        if (is_new_hamiltonianmatrix_required) {
            if (is_interaction_already_contained) {

                // Check variables for consistency
                if (memory_saving || coefficients_unperturbed_cache.size() == 0 || hamiltonianmatrix_unperturbed_cache.size() == 0) {
                    throw std::runtime_error( "Inconsistent variables at " + std::string(__FILE__) + ":" + std::to_string(__LINE__) + ".");
                }

                // Reset the Hamiltonian if it already contains the interaction
                coefficients = coefficients_unperturbed_cache;
                hamiltonianmatrix = hamiltonianmatrix_unperturbed_cache;
            } else if (!memory_saving) {

                // Store the Hamiltonian without interaction
                coefficients_unperturbed_cache = coefficients;
                hamiltonianmatrix_unperturbed_cache = hamiltonianmatrix;
            }

            // Build interaction
            this->initializeInteraction();

            // Add interaction to the Hamiltonian
            this->addInteraction();

            if (memory_saving){
                // Delete the variables used to add the interaction to the Hamiltonian
                this->deleteInteraction();
            }

            is_interaction_already_contained = true;
            is_new_hamiltonianmatrix_required = false;
        }
    }

    void buildInteraction() {
        // Build basis
        this->buildBasis();

        // Initialize interaction
        this->initializeInteraction(); // this method checks by itself whether a new initialization is required
    }

    void buildBasis() {
        // Check variables for consistency
        if ( ((hamiltonianmatrix.size() == 0) != states.empty()) || ((hamiltonianmatrix.size() == 0) != (coefficients.size() == 0)) ) {
            throw std::runtime_error( "Inconsistent variables at " + std::string(__FILE__) + ":" + std::to_string(__LINE__) + ".");
        }

        // In case of no new basis restrictions and already initialized basis, there is nothing to do
        if (states.size() != 0 && range_n.empty() && range_l.empty() && range_j.empty() && range_m.empty() &&
                energy_min == std::numeric_limits<double>::lowest() && energy_max == std::numeric_limits<double>::max()) { // TODO check for new threshold, too
            return;
        }

        // Check whether the basis does not exist
        if (states.size() == 0) {

            // Initialize the basis
            this->initializeBasis();

            // Forget basis restrictions as they were applied now
            this->forgetRestrictions();

        } else {
            this->updateEverything();
        }

        // Add dummy states
        if (states_artifical.size() > 0) { // TODO [dummystates]
            size_t row = coefficients.rows();
            size_t col = coefficients.cols();
            hamiltonianmatrix.conservativeResize(hamiltonianmatrix.rows()+states_artifical.size(), hamiltonianmatrix.cols()+states_artifical.size());
            coefficients.conservativeResize(coefficients.rows()+states_artifical.size(), coefficients.cols()+states_artifical.size());

            states.reserve(states.size()+states_artifical.size());
            coefficients.reserve(coefficients.size()+states_artifical.size());
            for (const auto &state : states_artifical) {
                states.push_back(state);
                coefficients.insert(row++, col++) = 1;
            }
            coefficients.makeCompressed();
            states_artifical.clear();
        }

        // Check whether the basis is empty
        if (coefficients.rows() == 0) {
            throw std::runtime_error( "The basis is contains no states." );
        }
        if (coefficients.cols() == 0) {
            throw std::runtime_error( "The basis is contains no vectors." );
        }
    }

    void diagonalize() {
        this->buildHamiltonian();

        // Check if already diagonal
        if (checkIsDiagonal(hamiltonianmatrix)) return;

        // Diagonalize hamiltonianmatrix
        eigen_dense_t densemat(hamiltonianmatrix);
        Eigen::SelfAdjointEigenSolver<eigen_dense_t> eigensolver(densemat);
        // TODO improve storage management, e.g. delete densemat right away

        // Get eigenvalues and eigenvectors
        eigen_vector_double_t evals = eigensolver.eigenvalues();
        eigen_sparse_t evecs = eigensolver.eigenvectors().sparseView(1e-4,0.5); // TODO use threshold dependence squared

        // Build the new hamiltonianmatrix
        hamiltonianmatrix.setZero();
        hamiltonianmatrix.reserve(evals.size());
        for (int idx = 0; idx < evals.size(); ++idx) {
            hamiltonianmatrix.insert(idx, idx) = evals.coeffRef(idx);
        }
        hamiltonianmatrix.makeCompressed();

        // Transform the basis vectors
        coefficients = (coefficients * evecs).pruned(1e-4,0.5); // TODO use threshold dependence
    }

    void canonicalize() {
        this->buildHamiltonian();

        // Transform the hamiltonianmatrix
        hamiltonianmatrix = coefficients*hamiltonianmatrix*coefficients.adjoint();

        // Transform the basis vectors
        coefficients = coefficients * coefficients.adjoint();
    }

    // TODO rotate()

    ////////////////////////////////////////////////////////////////////
    /// Methods to manipulate individual entries of the Hamiltonian ////
    ////////////////////////////////////////////////////////////////////

    size_t getStateindex(const T &state) {
        this->buildBasis();

        return std::distance(states.begin(), std::find(states.begin(), states.end(), state));
    }

    scalar_t getHamiltonianentry(const T &state_row, const T &state_col) {
        this->buildHamiltonian();

        size_t idx_row = getStateindex(state_row);
        size_t idx_col = getStateindex(state_col);

        eigen_sparse_t tmp = coefficients*hamiltonianmatrix*coefficients.adjoint(); // TODO check whether canonicalization successful by calculating checkIsDiagonal((coefficients*coefficients.adjoint()).prune()) TODO [dummystates]

        return tmp.coeff(idx_row, idx_col);

    }

    void setHamiltonianentry(const T &state_row, const T &state_col, scalar_t value) {
        this->buildHamiltonian();

        size_t idx_row = getStateindex(state_row);
        size_t idx_col = getStateindex(state_col);

        eigen_sparse_t tmp = coefficients*hamiltonianmatrix*coefficients.adjoint(); // TODO check whether canonicalization successful by calculating checkIsDiagonal((coefficients*coefficients.adjoint()).prune()) TODO [dummystates]
        tmp.coeffRef(idx_row, idx_col) = value; // TODO check whether this also works if the element does not exist TODO [dummystates]
        if (idx_row != idx_col) tmp.coeffRef(idx_col, idx_row) = this->conjugate(value);

        hamiltonianmatrix = coefficients.adjoint()*tmp*coefficients;
    }

    void addHamiltonianentry(const T &state_row, const T &state_col, scalar_t value) {
        this->buildHamiltonian();

        size_t idx_row = getStateindex(state_row);
        size_t idx_col = getStateindex(state_col);

        eigen_sparse_t tmp(states.size(), states.size());
        tmp.reserve(2);
        tmp.insert(idx_row, idx_col) = value;
        if (idx_row != idx_col) tmp.insert(idx_col, idx_row) = this->conjugate(value);
        tmp.makeCompressed();

        hamiltonianmatrix += coefficients.adjoint()*tmp*coefficients;
    }

protected:
    SystemBase(std::string cachedir) : cachedir(boost::filesystem::absolute(cachedir)), energy_min(std::numeric_limits<double>::lowest()), energy_max(std::numeric_limits<double>::max()), memory_saving(false), is_interaction_already_contained(false), is_new_hamiltonianmatrix_required(false) {
    }

    SystemBase(std::string cachedir, bool memory_saving) : cachedir(boost::filesystem::absolute(cachedir)), energy_min(std::numeric_limits<double>::lowest()), energy_max(std::numeric_limits<double>::max()), memory_saving(memory_saving), is_interaction_already_contained(false), is_new_hamiltonianmatrix_required(false) {
    }

    SystemBase() : cachedir(""), energy_min(std::numeric_limits<double>::lowest()), energy_max(std::numeric_limits<double>::max()), memory_saving(false), is_interaction_already_contained(false), is_new_hamiltonianmatrix_required(false) {
    }

    SystemBase(bool memory_saving) : cachedir(""), energy_min(std::numeric_limits<double>::lowest()), energy_max(std::numeric_limits<double>::max()), memory_saving(memory_saving), is_interaction_already_contained(false), is_new_hamiltonianmatrix_required(false) {
    }

    virtual void initializeBasis() = 0;
    virtual void initializeInteraction() = 0;

    virtual void transformInteraction(const eigen_sparse_t &transformator) = 0;
    virtual void addInteraction() = 0;
    virtual void deleteInteraction() = 0;

    boost::filesystem::path cachedir;

    double energy_min, energy_max;
    std::set<int> range_n, range_l;
    std::set<float> range_j, range_m;

    bool memory_saving;
    bool is_interaction_already_contained;
    bool is_new_hamiltonianmatrix_required;

    std::vector<T> states;
    eigen_sparse_t coefficients;
    eigen_sparse_t hamiltonianmatrix;
    eigen_sparse_t coefficients_unperturbed_cache;
    eigen_sparse_t hamiltonianmatrix_unperturbed_cache;

    ////////////////////////////////////////////////////////////////////
    /// Helper method that shoul be called by the derived classes //////
    ////////////////////////////////////////////////////////////////////

    void onParameterChange () {
        // Check variables for consistency
        if ((coefficients_unperturbed_cache.size() == 0) != (hamiltonianmatrix_unperturbed_cache.size() == 0)) {
            throw std::runtime_error( "Inconsistent variables at " + std::string(__FILE__) + ":" + std::to_string(__LINE__) + ".");
        }

        // Throw error if the Hamiltonian cannot be changed anymore
        if (is_interaction_already_contained && coefficients_unperturbed_cache.size() == 0) {
            throw std::runtime_error( "If memory saving is activated, one cannot change parameters after interaction was added to the Hamiltonian." );
        }

        is_new_hamiltonianmatrix_required = true;
    }

    void onSymmetryChange () {
        // Throw error if the Symmetry cannot be changed anymore
        if (!states.empty()) {
            throw std::runtime_error( "One cannot change symmetries after the basis was built." );
        }
    }

    ////////////////////////////////////////////////////////////////////
    /// Helper method to check diagonality of a matrix /////////////////
    ////////////////////////////////////////////////////////////////////

    bool checkIsDiagonal (const eigen_sparse_t &mat) {
        for (int k=0; k<mat.outerSize(); ++k) {
            for (eigen_iterator_t triple(mat,k); triple; ++triple) {
                if (triple.row() != triple.col()) return false;
            }
        }
        return true;
    }

    ////////////////////////////////////////////////////////////////////
    /// Helper methods to check the validity of states /////////////////
    ////////////////////////////////////////////////////////////////////

    template<class V>
    bool checkIsQuantumnumberValid(V q, std::set<V> range_q) {
        return range_q.empty() || range_q.find(q) != range_q.end();
    }

    template<class V>
    bool checkIsQuantumnumberValid(std::array<V, 2> q, std::set<V> range_q) {
        return range_q.empty() || (range_q.find(q[0]) != range_q.end() && range_q.find(q[1]) != range_q.end());
    }

    bool checkIsQuantumstateValid(T state) {
        return checkIsQuantumnumberValid(state.n, range_n) && checkIsQuantumnumberValid(state.l, range_l) && checkIsQuantumnumberValid(state.j, range_j) && checkIsQuantumnumberValid(state.m, range_m);
    }

    bool checkIsEnergyValid(double e) {
        return (e > energy_min || energy_min == std::numeric_limits<double_t>::lowest()) && (e < energy_max  || energy_max == std::numeric_limits<double_t>::max()); // TODO take into account numerical errors
    }

    ////////////////////////////////////////////////////////////////////
    /// Helper methods to change the set of basis vectors //////////////
    ////////////////////////////////////////////////////////////////////

    void applyLeftsideTransformator(std::vector<eigen_triplet_t> triplets_transformator) {
        eigen_sparse_t transformator(triplets_transformator.size(),coefficients.rows());
        transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

        // Apply transformator in order to remove rows from the coefficient matrix (i.e. states)
        coefficients = transformator*coefficients;
        if (coefficients_unperturbed_cache.size() != 0) coefficients_unperturbed_cache = transformator*coefficients_unperturbed_cache;
    }

    void applyRightsideTransformator(std::vector<eigen_triplet_t> triplets_transformator) {
        eigen_sparse_t transformator(coefficients.cols(),triplets_transformator.size());
        transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

        // Apply transformator in order to remove columns from the coefficient matrix (i.e. basis vectors)
        coefficients = coefficients*transformator;
        if (coefficients_unperturbed_cache.size() != 0) coefficients_unperturbed_cache = coefficients_unperturbed_cache*transformator;

        // Apply transformator in order to remove rows and columns from the matrices that help constructing the total Hamiltonianmatrix
        this->transformInteraction(transformator);

        // Apply transformator in order to remove rows and columns from the total Hamiltonianmatrix
        hamiltonianmatrix = transformator.adjoint()*hamiltonianmatrix*transformator;
        if (hamiltonianmatrix_unperturbed_cache.size() != 0) hamiltonianmatrix_unperturbed_cache = transformator.adjoint()*hamiltonianmatrix_unperturbed_cache*transformator;
    }

    void removeRestrictedStates(std::function<bool(size_t)> checkIsValidIndex) {
        // Build transformator and remove states
        std::vector<T> states_new;
        states_new.reserve(states.size());
        std::vector<eigen_triplet_t> triplets_transformator;
        triplets_transformator.reserve(states.size());

        size_t idx_new = 0;
        for (size_t idx = 0; idx < states.size(); ++idx) {
            if (checkIsValidIndex(idx)) {
                states_new.push_back(states[idx]);
                triplets_transformator.push_back(eigen_triplet_t(idx_new++,idx,1));
            }
        }

        states_new.shrink_to_fit();
        states = states_new;
        this->applyLeftsideTransformator(triplets_transformator);
    }

    ////////////////////////////////////////////////////////////////////
    /// Utility methods ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    std::complex<double> conjugate(const std::complex<double> &val) {
        return std::conj(val);
    }

    double conjugate(const double &val) {
        return val;
    }

    double real(const double &val) {
        return val;
    }

    double real(const std::complex<double> &val) {
        return val.real();
    }

    template<class V>
    void range(std::set<V> &rset, V rmin, V rmax) {
        rset.clear();
        for (V r = rmin; r <= rmax; ++r) {
            rset.insert(r);
        }
    }

private:
    std::vector<T> states_artifical; // TODO think about whether one can use here real states, too TODO [dummystates]

    void forgetRestrictions() {
        energy_min = std::numeric_limits<double>::lowest();
        energy_max = std::numeric_limits<double>::max();
        range_n.clear();
        range_l.clear();
        range_j.clear();
        range_m.clear();
    }

    ////////////////////////////////////////////////////////////////////
    /// Method to update the system ////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    void updateEverything() {

        if (!range_n.empty() || !range_l.empty() || !range_j.empty() || !range_m.empty()) {

            ////////////////////////////////////////////////////////////////////
            /// Remove restricted states ///////////////////////////////////////
            ////////////////////////////////////////////////////////////////////

            // Remove states if the quantum numbers are not allowed
            removeRestrictedStates([=](size_t idx) -> bool { return this->checkIsQuantumstateValid(states[idx]); } );
        }

        if (!range_n.empty() || !range_l.empty() || !range_j.empty() || !range_m.empty() ||
                energy_min != std::numeric_limits<double>::lowest() ||
                energy_max != std::numeric_limits<double>::max()) { // TODO also check for a new value of threshold_for_sqnorm

            ////////////////////////////////////////////////////////////////////
            /// Remove vectors with restricted energies or too small norm //////
            ////////////////////////////////////////////////////////////////////

            // Build transformator and remove vectors (if their energy is not allowed or the squared norm to small)
            std::vector<eigen_triplet_t> triplets_transformator;
            triplets_transformator.reserve(coefficients.cols());

            size_t idx_new = 0;
            for (int idx=0; idx<coefficients.cols(); ++idx) { // idx = col = num basis vector

                if (checkIsEnergyValid(this->real(hamiltonianmatrix.coeff(idx, idx)))) {
                    double_t sqnorm = 0;

                    // Calculate the square norm of the columns of the coefficient matrix
                    for (eigen_iterator_t triple(coefficients,idx); triple; ++triple) {
                        sqnorm += std::pow(std::abs(triple.value()),2);
                    }
                    if (sqnorm > 0.05) {
                        triplets_transformator.push_back(eigen_triplet_t(idx,idx_new++,1));
                    }
                }
            }

            this->applyRightsideTransformator(triplets_transformator);

            ////////////////////////////////////////////////////////////////////
            /// Remove states that barely occur within the vectors /////////////
            ////////////////////////////////////////////////////////////////////

            // Calculate the square norm of the rows of the coefficient matrix
            std::vector<double> sqnorm_list(coefficients.rows(),0);
            for (int k=0; k<this->coefficients.cols(); ++k) {
                for (eigen_iterator_t triple(coefficients,k); triple; ++triple) {
                    sqnorm_list[triple.row()] += std::pow(std::abs(triple.value()),2);
                }
            }

            // Remove states if the squared norm is to small
            removeRestrictedStates([=](size_t idx) -> bool { return sqnorm_list[idx] > 0.05; } );
        }

        // Forget basis restrictions as they were applied now
        this->forgetRestrictions();
    }
};

#endif
