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

template<class T> class SystemBase {
public:
    SystemBase(std::vector<T> states) : // TODO
        energy_min(std::numeric_limits<double>::lowest()), energy_max(std::numeric_limits<double>::max()), states(states) {
        throw std::runtime_error( "Method yet not implemented." );
    }

    virtual ~SystemBase() = default;

    bool checkIsDiagonal() {
        if (hamiltonianmatrix.size() == 0 && energies.empty()) throw std::runtime_error( "The system is empty." );
        if (hamiltonianmatrix.size() != 0 && !energies.empty()) throw std::runtime_error( "The system contains redundant data." );
        return !energies.empty();
    }

    // TODO setThresholdForSqnorm()

    ////////////////////////////////////////////////////////////////////
    /// Methods to restrict the number of states inside the basis //////
    ////////////////////////////////////////////////////////////////////

    void restrictEnergy(double e_min, double e_max) {
        energy_min = e_min;
        energy_max = e_max;
    }
    void restrictN(int n_min, int n_max) {
        n_min = std::max(n_min, 1);
        for (int n = n_min; n <= n_max; ++n) {
            range_n.insert(n);
        }
    }
    void restrictN(std::set<int> n) {
        range_n = n;
    }
    void restrictL(int l_min, int l_max) {
        l_min = std::max(l_min, 0);
        for (int l = l_min; l <= l_max; ++l) {
            range_l.insert(l);
        }
    }
    void restrictL(std::set<int> l) {
        range_l = l;
    }
    void restrictJ(float j_min, float j_max) {
        j_min = std::fmax(j_min, 0.5);
        for (float j = j_min; j <= j_max; ++j) {
            range_j.insert(j);
        }
    }
    void restrictJ(std::set<float> j) {
        range_j = j;
    }
    void restrictM(float m_min, float m_max) {
        for (float m = m_min; m <= m_max; ++m) {
            range_m.insert(m);
        }
    }
    void restrictM(std::set<float> m) {
        range_m = m;
    }
    const std::vector<T>& getStates() {
        this->build();
        return states;
    }
    const eigen_sparse_t& getCoefficients() {
        this->build();
        return coefficients;
    }
    const std::vector<double>& getEnergies() {
        if (!checkIsDiagonal()) throw std::runtime_error( "The method getEnergies() does not work in case of a not diagonalized system, use getHamiltonianmatrix() instead." );
        this->build();
        return energies;
    }
    const eigen_sparse_t& getHamiltonianmatrix() {
        if (checkIsDiagonal()) throw std::runtime_error( "The method getHamiltonianmatrix() does not work in case of a diagonalized system, use getEnergies() instead." );
        this->build();
        return hamiltonianmatrix;
    }

    ////////////////////////////////////////////////////////////////////
    /// Methods to build, transform, and destroy the system ////////////
    ////////////////////////////////////////////////////////////////////

    void build() {
        // Generate the basis from scratch or update the basis
        if (energies.empty() && hamiltonianmatrix.size() == 0
                && states.empty() && coefficients.size() == 0) {
            initialize();
        } else if (((checkIsDiagonal() && !energies.empty()) || (!checkIsDiagonal() && hamiltonianmatrix.size() != 0))
                   && !states.empty() && coefficients.size() != 0) {
            update();
        } else {
            throw std::runtime_error( "Method yet not implemented." ); // TODO implement a method that also works if !states.empty()
        }

        // TODO add method initialize, update rest of Hamiltonian

        // TODO also build the hamiltonianmatrix

        // Forget all previous restrictions as they were applied now
        this->forgetRestrictions();
    }

    void clear() {
        energies.clear();
        states.clear();
        coefficients.resize(0,0);
        this->forgetRestrictions();
    }

    // TODO diagonalize() should delete all the storage and make the initial Hamiltonianmatrix the original, one changeFields should throw an error
    // TODO canonicalize() should delete all the storage and make the initial Hamiltonianmatrix the original, one changeFields should throw an error
    // TODO rotate()

    ////////////////////////////////////////////////////////////////////
    /// Methods to get the size of the system //////////////////////////
    ////////////////////////////////////////////////////////////////////

    size_t getNumVectors() {
        // TODO check energies.size() == coefficients.outerSize() == coefficients.cols() == hamiltonianmatrix.rows() == hamiltonianmatrix.cols()
        return coefficients.cols();
    }

    size_t getNumStates() {
        // TODO check states.size() == coefficients.innerSize() == coefficients.rows()
        return coefficients.rows();
    }

protected:
    virtual void initialize() = 0;

    SystemBase() : energy_min(std::numeric_limits<double>::lowest()), energy_max(std::numeric_limits<double>::max()) {
    }

    void forgetRestrictions() {
        energy_min = std::numeric_limits<double>::lowest();
        energy_max = std::numeric_limits<double>::max();
        range_n.clear();
        range_l.clear();
        range_j.clear();
        range_m.clear();
    }

    double energy_min, energy_max;
    std::set<int> range_n, range_l;
    std::set<float> range_j, range_m;

    std::vector<double> energies;
    std::vector<T> states;
    eigen_sparse_t coefficients;
    eigen_sparse_t hamiltonianmatrix;

private:

    ////////////////////////////////////////////////////////////////////
    /// Method to update the system ////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    void update() {
        if (!range_n.empty() || !range_l.empty() || !range_j.empty() || !range_m.empty()) {

            ////////////////////////////////////////////////////////////////////
            /// Remove restricted states ///////////////////////////////////////
            ////////////////////////////////////////////////////////////////////

            // Build transformator and remove states (if the quantum numbers are not allowed)
            std::vector<T> states_new;
            states_new.reserve(this->getNumStates());
            std::vector<eigen_triplet_t> triplets_transformator;
            triplets_transformator.reserve(this->getNumStates());

            size_t idx_new = 0;
            for (size_t idx = 0; idx < this->getNumStates(); ++idx) {
                if (checkIsQuantumstateValid(states[idx])) {
                    states_new.push_back(states[idx]);
                    triplets_transformator.push_back(eigen_triplet_t(idx_new++,idx,1));
                }
            }

            states_new.shrink_to_fit();
            states = states_new;
            this->applyLeftsideTransformator(triplets_transformator);
        }

        if (!range_n.empty() || !range_l.empty() || !range_j.empty() || !range_m.empty() ||
                energy_min != std::numeric_limits<double>::lowest() ||
                energy_max != std::numeric_limits<double>::max()) { // TODO also check for a new value of threshold_for_sqnorm

            ////////////////////////////////////////////////////////////////////
            /// Remove vectors with restricted energies or too small norm //////
            ////////////////////////////////////////////////////////////////////

            // Build transformator and remove energies (if the qenergy is not allowed or the squared norm to small)
            std::vector<double> energies_new;
            if (checkIsDiagonal()) energies_new.reserve(this->getNumVectors());
            std::vector<eigen_triplet_t> triplets_transformator;
            triplets_transformator.reserve(this->getNumVectors());

            size_t idx_new = 0;
            for (size_t idx=0; idx<this->getNumVectors(); ++idx) { // idx = col = num basis vector

                if (!checkIsDiagonal() || checkIsEnergyValid(energies[idx])) {
                    double_t sqnorm = 0;

                    // Calculate the square norm of the columns of the coefficient matrix
                    for (eigen_iterator_t triple(coefficients,idx); triple; ++triple) {
                        sqnorm += std::pow(std::abs(triple.value()),2);
                    }
                    if (sqnorm > 0.05) {
                        triplets_transformator.push_back(eigen_triplet_t(idx,idx_new++,1));
                        if (checkIsDiagonal()) energies_new.push_back(energies[idx]);
                    }
                }
            }

            if (checkIsDiagonal()) {
                energies_new.shrink_to_fit();
                energies = energies_new;
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

            // Build transformator and remove states (if the squared norm is to small)
            std::vector<T> states_new;
            states_new.reserve(this->getNumStates());
            triplets_transformator.clear();
            triplets_transformator.reserve(this->getNumStates());

            idx_new = 0;
            for (size_t idx = 0; idx < this->getNumStates(); ++idx) {
                if (sqnorm_list[idx] > 0.05) {
                    states_new.push_back(states[idx]);
                    triplets_transformator.push_back(eigen_triplet_t(idx_new++,idx,1));
                }
            }

            states_new.shrink_to_fit();
            states = states_new;
            this->applyLeftsideTransformator(triplets_transformator);
        }
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
        return (e > energy_min || energy_min == std::numeric_limits<double_t>::lowest()) && (e < energy_max  || energy_max == std::numeric_limits<double_t>::max());
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
        if (!checkIsDiagonal()) hamiltonianmatrix = transformator.adjoint()*hamiltonianmatrix*transformator;

        // Apply transformator in order to remove rows and columns from the matrices that help constructing the total Hamiltonianmatrix
        // TODO
    }

};

#endif
