#ifndef SYSTEMBASE_H
#define SYSTEMBASE_H

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

template<class T> class SystemBase {
public:
    SystemBase(std::vector<T> states, std::string cachedir) : // TODO
        cachedir(boost::filesystem::absolute(cachedir)), energy_min(std::numeric_limits<double>::lowest()), energy_max(std::numeric_limits<double>::max()), states(states), hamiltonianhelpers_missing(true), hamiltonian_missing(true) {
        throw std::runtime_error( "Method yet not implemented." );
    }

    virtual ~SystemBase() = default;

    // TODO setThresholdForSqnorm()

    ////////////////////////////////////////////////////////////////////
    /// Methods to restrict the number of states inside the basis //////
    ////////////////////////////////////////////////////////////////////

    void restrictEnergy(double e_min, double e_max) { // restricts diagonal entries
        energy_min = e_min;
        energy_max = e_max;
    }

    void restrictN(int n_min, int n_max) {
        n_min = std::max(n_min, 1);
        range_n.clear();
        for (int n = n_min; n <= n_max; ++n) {
            range_n.insert(n);
        }
    }

    void restrictN(std::set<int> n) {
        range_n = n;
    }

    void restrictL(int l_min, int l_max) {
        l_min = std::max(l_min, 0);
        range_l.clear();
        for (int l = l_min; l <= l_max; ++l) {
            range_l.insert(l);
        }
    }

    void restrictL(std::set<int> l) {
        range_l = l;
    }

    void restrictJ(float j_min, float j_max) {
        j_min = std::fmax(j_min, 0.5);
        range_j.clear();
        for (float j = j_min; j <= j_max; ++j) {
            range_j.insert(j);
        }
    }

    void restrictJ(std::set<float> j) {
        range_j = j;
    }

    void restrictM(float m_min, float m_max) {
        range_m.clear();
        for (float m = m_min; m <= m_max; ++m) {
            range_m.insert(m);
        }
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

    const eigen_sparse_t& getCoefficients() {
        this->buildBasis();
        return coefficients;
    }

    const std::vector<double> getDiagonal() {
        this->buildHamiltonian();
        eigen_vector_double_t diagonal = hamiltonianmatrix.diagonal().real();
        return std::vector<double>(diagonal.data(), diagonal.data()+diagonal.size());
    }

    const eigen_sparse_t& getHamiltonianmatrix() {
        this->buildHamiltonian();
        return hamiltonianmatrix;
    }

    size_t getNumVectors() {
        // coefficients.outerSize() == coefficients.cols() == hamiltonianmatrix.rows() == hamiltonianmatrix.cols()
        // TODO check for new restrictions and not yet build coefficients, while keeping getNumVectors work in the buildBasis method --> use coefficients cols in the buildBasis method?
        return coefficients.cols();
    }

    size_t getNumStates() {
        // coefficients.innerSize() == coefficients.rows() == states.size()
        // TODO checkIsNewRestrictions
        // TODO check for new restrictions and not yet build coefficients, while keeping getNumStates work in the buildBasis method --> use coefficients rows in the buildBasis method?
        return coefficients.rows();
    }

    ////////////////////////////////////////////////////////////////////
    /// Methods to build, transform, and destroy the system ////////////
    ////////////////////////////////////////////////////////////////////

    void buildHamiltonian() {
        // After executing this method, the interaction strength is fixed // TODO throw error if a change is done nevertheless
        // This method overrides the variable "hamiltonianmatrix" with the total hamiltonian
        if (hamiltonian_missing) {
            this->buildHamiltonianhelpers();
            this->initializeHamiltonian();
            hamiltonian_missing = false;
        } else {
            this->updateEverything();
        }

        // Forget basis restrictions as they were applied now
        this->forgetRestrictions();
    }

    void buildHamiltonianhelpers() {
        // After executing this method, the interaction form is fixed // TODO throw error if a change is done nevertheless
        if (hamiltonianhelpers_missing) {
            this->buildBasis();
            this->initializeHamiltonianhelpers();
            hamiltonianhelpers_missing = false;
        } else {
            this->updateEverything();
        }
    }

    void buildBasis() {
        // Handle the case of no basis restrictions
        if (range_n.empty() && range_l.empty() && range_j.empty() && range_m.empty() &&
                energy_min == std::numeric_limits<double>::lowest() && energy_max == std::numeric_limits<double>::max()) {
            // TODO if basis not initialized, throw error infinite basis; else return
            return;
        }

        // After executing this method, the element fixed // TODO throw error if a change is done nevertheless
        if (hamiltonianmatrix.size() == 0 && states.empty() && coefficients.size() == 0) {
            this->initializeBasis();
            // Forget basis restrictions as they were applied now
            this->forgetRestrictions();
        } else if (hamiltonianmatrix.size() != 0 && !states.empty() && coefficients.size() != 0) {
            this->updateEverything();
        } else {
            throw std::runtime_error( "Method yet not implemented." ); // TODO implement a method that also works if !states.empty()
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

    // void addCoupling(T state_row, T state_col, scalar_t, strength)

    void clear() {
        states.clear();
        hamiltonianmatrix.resize(0,0);
        coefficients.resize(0,0);
        this->forgetRestrictions();
    }

protected:
    SystemBase(std::string cachedir) : cachedir(boost::filesystem::absolute(cachedir)), energy_min(std::numeric_limits<double>::lowest()), energy_max(std::numeric_limits<double>::max()), hamiltonianhelpers_missing(true), hamiltonian_missing(true) {
    }

    virtual void initializeBasis() = 0;
    virtual void initializeHamiltonianhelpers() = 0;
    virtual void initializeHamiltonian() = 0; // this method should also delete the hamiltonian helpers as afterwards, they won't be needed anymore

    virtual void transformHamiltonianhelpers(const eigen_sparse_t &transformator) = 0;

    boost::filesystem::path cachedir;

    double energy_min, energy_max;
    std::set<int> range_n, range_l;
    std::set<float> range_j, range_m;

    eigen_sparse_t hamiltonianmatrix;
    std::vector<T> states;
    eigen_sparse_t coefficients;

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
    }

    void applyRightsideTransformator(std::vector<eigen_triplet_t> triplets_transformator) {
        eigen_sparse_t transformator(coefficients.cols(),triplets_transformator.size());
        transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

        // Apply transformator in order to remove columns from the coefficient matrix (i.e. basis vectors)
        coefficients = coefficients*transformator;

        // Apply transformator in order to remove rows and columns from the current Hamiltonianmatrix
        hamiltonianmatrix = transformator.adjoint()*hamiltonianmatrix*transformator;

        // Apply transformator in order to remove rows and columns from the matrices that help constructing the total Hamiltonianmatrix
        if (!hamiltonianhelpers_missing && hamiltonian_missing) this->transformHamiltonianhelpers(transformator); // if the hamiltonian is not missing anymore, the helpers are not needed
    }

    void removeRestrictedStates(std::function<bool(size_t)> checkIsValidIndex) {
        // Build transformator and remove states
        std::vector<T> states_new;
        states_new.reserve(this->getNumStates());
        std::vector<eigen_triplet_t> triplets_transformator;
        triplets_transformator.reserve(this->getNumStates());

        size_t idx_new = 0;
        for (size_t idx = 0; idx < this->getNumStates(); ++idx) {
            if (checkIsValidIndex(idx)) {
                states_new.push_back(states[idx]);
                triplets_transformator.push_back(eigen_triplet_t(idx_new++,idx,1));
            }
        }

        states_new.shrink_to_fit();
        states = states_new;
        this->applyLeftsideTransformator(triplets_transformator);
    }

private:
    bool hamiltonianhelpers_missing;
    bool hamiltonian_missing;

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
            triplets_transformator.reserve(this->getNumVectors());

            size_t idx_new = 0;
            for (size_t idx=0; idx<this->getNumVectors(); ++idx) { // idx = col = num basis vector

                if (checkIsEnergyValid(std::complex<double>(hamiltonianmatrix.coeff(idx, idx)).real())) { // std::complex... makes this code work in case of a real hamiltonianmatrix, too
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
            std::vector<double> sqnorm_list(this->getNumStates(),0);
            for (size_t k=0; k<this->getNumVectors(); ++k) {
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
