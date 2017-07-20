#ifndef SYSTEMBASE_H
#define SYSTEMBASE_H

#include "dtypes.h"
#include "State.h"
#include "serialization_eigen.h"
#include "serialization_path.h"
#include "WignerD.h"

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
#include <unordered_map>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>

template<class T>
class enumerated_state {
public:
    enumerated_state(size_t idx, T state) : idx(idx), state(state) {
    }
    enumerated_state() : idx(0), state() { // TODO remove and use http://www.boost.org/doc/libs/1_46_1/libs/serialization/doc/serialization.html#constructors instead
    }
    size_t idx;
    T state;

private:
    ////////////////////////////////////////////////////////////////////
    /// Method for serialization ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        (void)version;

        ar & idx & state;
    }
};

template<class T>
using states_set = boost::multi_index_container<
enumerated_state<T>,
boost::multi_index::indexed_by<
boost::multi_index::random_access<>,
boost::multi_index::hashed_unique<boost::multi_index::member<enumerated_state<T>,T,&enumerated_state<T>::state>, std::hash<T> >
>
>;

template<class T> class SystemBase {
public:
    SystemBase(std::vector<T> states, std::wstring cachedir) : // TODO
        cachedir(boost::filesystem::absolute(cachedir)), energy_min(std::numeric_limits<double>::lowest()), energy_max(std::numeric_limits<double>::max()), memory_saving(false), is_interaction_already_contained(false), is_new_hamiltonianmatrix_required(false) {
        (void)states;
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
    
    const std::vector<T> getStates() { // TODO @hmenke typemap for "const state_set<T>&"
        this->buildBasis();
        std::vector<T> states_converted;
        states_converted.reserve(states.size());
        for (const auto& s: states) {
            states_converted.push_back(s.state);
        }
        return states_converted;
    }
    
    eigen_sparse_t& getCoefficients() {
        this->buildBasis();
        return coefficients;
    }

    eigen_vector_double_t getDiagonal() {
        this->buildHamiltonian();
        return hamiltonianmatrix.diagonal().real();
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

    eigen_vector_double_t getOverlap(const T &generalizedstate) {
        // Build basis
        this->buildBasis();

        // Determine indices of the specified states
        std::vector<size_t> states_indices;
        for (const auto &state: states) {
            if (state.state ^ generalizedstate) { // Check whether state.state is contained in generalizedstate
                states_indices.push_back(state.idx);
            }
        }

        // Build state vectors
        std::vector<eigen_triplet_t> overlap_states_triplets;
        overlap_states_triplets.reserve(states_indices.size());

        size_t current = 0;
        for (auto const &idx: states_indices) {
            overlap_states_triplets.push_back(eigen_triplet_t(idx, current++, 1));
        }

        eigen_sparse_t overlap_states(states.size(), states_indices.size());
        overlap_states.setFromTriplets(overlap_states_triplets.begin(), overlap_states_triplets.end());
        overlap_states_triplets.clear();

        // Calculate overlap
        eigen_sparse_t product = coefficients.adjoint()*overlap_states;
        eigen_vector_double_t overlap = eigen_vector_double_t::Zero(product.rows());
        for (int k=0; k<product.outerSize(); ++k) {
            for (eigen_iterator_t triple(product, k); triple; ++triple) {
                overlap[triple.row()] += std::pow(std::abs(triple.value()),2);
            }
        }

        return overlap;
    }

    eigen_vector_double_t getOverlap(const T &state, std::array<double, 3> to_z_axis, std::array<double, 3> to_y_axis) {
        // Get Euler angles
        Eigen::Matrix<double,3,3> rotator = this->buildRotator(to_z_axis, to_y_axis);
        Eigen::Matrix<double,3,1> euler_zyz = rotator.eulerAngles(2, 1, 2);
        double alpha = euler_zyz[0];
        double beta = euler_zyz[1];
        double gamma = euler_zyz[2];

        return this->getOverlap(state, alpha, beta, gamma);
    }

    eigen_vector_double_t getOverlap(const T &generalizedstate, double alpha, double beta, double gamma) {
        // Build basis
        this->buildBasis();

        // Determine indices of the specified states
        std::vector<size_t> states_indices;
        for (const auto &state: states) {
            if (state.state ^ generalizedstate) { // Check whether state.state is contained in generalizedstate
                states_indices.push_back(state.idx);
            }
        }

        // Build rotated state vectors
        eigen_sparse_t overlap_states_rotated = this->rotateStates(states_indices, alpha, beta, gamma);

        // Calculate overlap
        eigen_sparse_t product = coefficients.adjoint()*overlap_states_rotated;
        eigen_vector_double_t overlap = eigen_vector_double_t::Zero(product.rows());
        for (int k=0; k<product.outerSize(); ++k) {
            for (eigen_iterator_t triple(product, k); triple; ++triple) {
                overlap[triple.row()] += std::pow(std::abs(triple.value()),2);
            }
        }

        return overlap;
    }

    std::vector<T> getMainStates() {
        // Build basis
        this->buildBasis();

        // Get states with the main contribution
        std::vector<T> states_with_maxval;
        states_with_maxval.reserve(coefficients.cols());

        for (int k=0; k<coefficients.outerSize(); ++k) { // col == idx_vector
            double maxval = -1;
            size_t row_with_maxval;

            for (eigen_iterator_t triple(coefficients, k); triple; ++triple) {
                if (std::abs(triple.value()) > maxval) {
                    row_with_maxval = triple.row();
                    maxval = std::abs(triple.value());
                }
            }

            states_with_maxval.push_back(states[row_with_maxval].state);
        }

        return states_with_maxval;
    }

    std::array<std::vector<size_t>,2> getConnections(SystemBase<T> &system_to, double threshold) {
        // Build basis
        this->buildBasis();

        double threshold_sqrt = std::sqrt(threshold);

        // Calculate transformator between the set of states
        std::vector<eigen_triplet_t> triplets_transformator;
        triplets_transformator.reserve(std::min(this->getNumStates(),system_to.getNumStates()));

        for (size_t idx_to = 0; idx_to < system_to.getNumStates(); ++idx_to) {
            const T &state = system_to.getStates()[idx_to];
            auto state_iter =  states.template get<1>().find(state);
            if (state_iter != states.template get<1>().end()) {
                size_t idx_from = state_iter->idx;
                triplets_transformator.push_back(eigen_triplet_t(idx_from,idx_to,1));
            }
        }

        eigen_sparse_t transformator(this->getNumStates(),system_to.getNumStates());
        transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

        // Calculate overlap
        eigen_sparse_double_t overlap_sqrt = (coefficients.adjoint()*transformator*system_to.getCoefficients()).cwiseAbs();

        // Determine the indices of the maximal values within the rows
        std::vector<size_t> rows_with_maxval;
        rows_with_maxval.reserve(overlap_sqrt.cols());

        for (int k=0; k<overlap_sqrt.outerSize(); ++k) {
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
        std::array<std::vector<size_t>,2> connections;
        connections[0].reserve(std::max(overlap_sqrt.rows(), overlap_sqrt.cols()));
        connections[1].reserve(std::max(overlap_sqrt.rows(), overlap_sqrt.cols()));

        for (int k=0; k<overlap_sqrt_transposed.outerSize(); ++k) { // cols
            double maxval = threshold_sqrt;
            size_t idx_from;
            size_t idx_to;

            for (eigen_iterator_double_t triple(overlap_sqrt_transposed,k); triple; ++triple) {
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
                states.push_back(enumerated_state<T>(row, state));
                coefficients.insert(row++, col++) = 1;
            }
            coefficients.makeCompressed();
            states_artifical.clear();
        }
        
        // Check whether the basis is empty
        if (coefficients.rows() == 0) {
            throw std::runtime_error( "The basis contains no states." );
        }
        if (coefficients.cols() == 0) {
            throw std::runtime_error( "The basis contains no vectors." );
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

    void rotate(std::array<double, 3> to_z_axis, std::array<double, 3> to_y_axis) {
        // Get Euler angles
        Eigen::Matrix<double,3,3> rotator = this->buildRotator(to_z_axis, to_y_axis);
        Eigen::Matrix<double,3,1> euler_zyz = rotator.eulerAngles(2, 1, 2);
        double alpha = euler_zyz[0];
        double beta = euler_zyz[1];
        double gamma = euler_zyz[2];

        this->rotate(alpha, beta, gamma);
    }

    void rotate(double alpha, double beta, double gamma) {
        // Build basis
        this->buildBasis();

        // Get the rotator for the basis states
        eigen_sparse_t transformator = this->buildStaterotator(alpha, beta, gamma);

        // Rotate basis
        this->transformInteraction(coefficients.adjoint());

        coefficients = transformator * coefficients;
        if (coefficients_unperturbed_cache.size() != 0) coefficients_unperturbed_cache = transformator*coefficients_unperturbed_cache;

        this->transformInteraction(coefficients);
    }

    void derotate(std::array<double, 3> to_z_axis, std::array<double, 3> to_y_axis) {
        // Get Euler angles
        Eigen::Matrix<double,3,3> rotator = this->buildRotator(to_z_axis, to_y_axis);
        Eigen::Matrix<double,3,1> euler_zyz = rotator.eulerAngles(2, 1, 2);
        double alpha = euler_zyz[0];
        double beta = euler_zyz[1];
        double gamma = euler_zyz[2];

        this->derotate(alpha, beta, gamma);
    }

    void derotate(double alpha, double beta, double gamma) {
        // Build basis
        this->buildBasis();

        // Get the rotator for the basis states
        eigen_sparse_t transformator = this->buildStaterotator(-alpha, -beta, -gamma);

        // Rotate basis
        this->transformInteraction(coefficients.adjoint());

        coefficients = transformator * coefficients;
        if (coefficients_unperturbed_cache.size() != 0) coefficients_unperturbed_cache = transformator*coefficients_unperturbed_cache;

        this->transformInteraction(coefficients);
    }

    void add(SystemBase<T> &system) {
        // --- Build Hamiltonians ---
        this->buildHamiltonian();
        system.buildHamiltonian();

        // --- Combine system specific variables ---
        this->incorporate(system);

        // --- Combine universal variables ---
        if (memory_saving != system.memory_saving) throw std::runtime_error("The value of the variable 'memory_saving' must be the same for both systems.");
        if (is_interaction_already_contained != system.is_interaction_already_contained) throw std::runtime_error("The value of the variable 'is_interaction_already_contained' must be the same for both systems.");
        if (is_new_hamiltonianmatrix_required != system.is_new_hamiltonianmatrix_required) throw std::runtime_error("The value of the variable 'is_new_hamiltonianmatrix_required' must be the same for both systems.");

        // --- Combine states and build transformators ---
        std::vector<eigen_triplet_t> transformator_triplets;
        transformator_triplets.reserve(system.states.size());

        for (const auto &entry : system.states) {
            auto state_iter =  states.template get<1>().find(entry.state);

            size_t newidx = states.size();
            if (state_iter == states.template get<1>().end()) {
                states.push_back(enumerated_state<T>(newidx, entry.state));
            } else {
                newidx = state_iter->idx;
            }

            transformator_triplets.push_back(eigen_triplet_t(newidx, entry.idx, 1));
        }

        eigen_sparse_t transformator(states.size(), system.coefficients.rows());
        transformator.setFromTriplets(transformator_triplets.begin(), transformator_triplets.end());
        transformator_triplets.clear();

        // --- Combine coefficients and coefficients_unperturbed_cache ---
        coefficients.conservativeResize(states.size(), coefficients.cols()+system.coefficients.cols());
        coefficients.rightCols(system.coefficients.cols()) = transformator*system.coefficients;

        if ((coefficients_unperturbed_cache.size() != 0) != (system.coefficients_unperturbed_cache.size() != 0)) throw std::runtime_error( "Inconsistent variables at " + std::string(__FILE__) + ":" + std::to_string(__LINE__) + ".");

        if (coefficients_unperturbed_cache.size() != 0) {
            coefficients_unperturbed_cache.conservativeResize(states.size(), coefficients_unperturbed_cache.cols()+system.coefficients_unperturbed_cache.cols());
            coefficients_unperturbed_cache.rightCols(system.coefficients_unperturbed_cache.cols()) = transformator*system.coefficients_unperturbed_cache;
        }

        // --- Combine hamiltonianmatrix and hamiltonianmatrix_unperturbed_cache ---
        std::vector<eigen_triplet_t> shifter_triplets;
        shifter_triplets.reserve(system.hamiltonianmatrix.rows());

        for (size_t idx = 0; idx < system.hamiltonianmatrix.rows(); ++idx) {
            shifter_triplets.push_back(eigen_triplet_t(hamiltonianmatrix.rows()+idx, idx, 1));
        }

        eigen_sparse_t shifter(hamiltonianmatrix.rows()+system.hamiltonianmatrix.rows(), system.hamiltonianmatrix.rows());
        shifter.setFromTriplets(shifter_triplets.begin(), shifter_triplets.end());
        shifter_triplets.clear();

        hamiltonianmatrix.conservativeResize(hamiltonianmatrix.rows()+system.hamiltonianmatrix.rows(), hamiltonianmatrix.cols()+system.hamiltonianmatrix.cols());
        hamiltonianmatrix.rightCols(system.hamiltonianmatrix.cols()) = shifter*system.hamiltonianmatrix;

        if ((hamiltonianmatrix_unperturbed_cache.size() != 0) != (system.hamiltonianmatrix_unperturbed_cache.size() != 0)) throw std::runtime_error( "Inconsistent variables at " + std::string(__FILE__) + ":" + std::to_string(__LINE__) + ".");

        if (hamiltonianmatrix_unperturbed_cache.size() != 0) {
            hamiltonianmatrix_unperturbed_cache.conservativeResize(hamiltonianmatrix_unperturbed_cache.rows()+system.hamiltonianmatrix_unperturbed_cache.rows(), hamiltonianmatrix_unperturbed_cache.cols()+system.hamiltonianmatrix_unperturbed_cache.cols());
            hamiltonianmatrix_unperturbed_cache.rightCols(system.hamiltonianmatrix_unperturbed_cache.cols()) = shifter*system.hamiltonianmatrix_unperturbed_cache;
        }
    }

    ////////////////////////////////////////////////////////////////////
    /// Methods to manipulate individual entries of the Hamiltonian ////
    ////////////////////////////////////////////////////////////////////
    
    size_t getStateindex(const T &state) {
        this->buildBasis();

        auto state_iter =  states.template get<1>().find(state);

        if (state_iter == states.template get<1>().end()) {
            throw std::runtime_error("The method getStateindex() is called on a non-existing state.");
        }

        return state_iter->idx;
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
    SystemBase(std::wstring cachedir) : cachedir(boost::filesystem::absolute(cachedir)), energy_min(std::numeric_limits<double>::lowest()), energy_max(std::numeric_limits<double>::max()), memory_saving(false), is_interaction_already_contained(false), is_new_hamiltonianmatrix_required(false) {
    }
    
    SystemBase(std::wstring cachedir, bool memory_saving) : cachedir(boost::filesystem::absolute(cachedir)), energy_min(std::numeric_limits<double>::lowest()), energy_max(std::numeric_limits<double>::max()), memory_saving(memory_saving), is_interaction_already_contained(false), is_new_hamiltonianmatrix_required(false) {
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

    virtual eigen_sparse_t rotateStates(const std::vector<size_t> &states_indices, double alpha, double beta, double gamma) = 0;
    virtual eigen_sparse_t buildStaterotator(double alpha, double beta, double gamma) = 0;
    virtual void incorporate(SystemBase<T> &system) = 0;

    boost::filesystem::path cachedir;
    
    double energy_min, energy_max;
    std::set<int> range_n, range_l;
    std::set<float> range_j, range_m;
    
    bool memory_saving;
    bool is_interaction_already_contained;
    bool is_new_hamiltonianmatrix_required;
    
    states_set<T> states;
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
    
    bool checkIsQuantumstateValid(const T &state) {
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
    
    void removeRestrictedStates(std::function<bool(const enumerated_state<T> &)> checkIsValidEntry) {
        // Build transformator and remove states
        states_set<T> states_new;
        states_new.reserve(states.size());
        std::vector<eigen_triplet_t> triplets_transformator;
        triplets_transformator.reserve(states.size());
        
        size_t idx_new = 0;

        for (const auto &entry: states) {
            if (checkIsValidEntry(entry)) {
                states_new.push_back(enumerated_state<T>(idx_new, entry.state));
                triplets_transformator.push_back(eigen_triplet_t(idx_new,entry.idx,1));
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

    Eigen::Matrix<double,3,3> buildRotator(std::array<double, 3> &to_z_axis, std::array<double, 3> &to_y_axis) {
        auto to_z_axis_mapped = Eigen::Map<Eigen::Matrix<double,3,1>>(&to_z_axis[0]).normalized();
        auto to_y_axis_mapped = Eigen::Map<Eigen::Matrix<double,3,1>>(&to_y_axis[0]).normalized();

        double tolerance = 1e-16;
        if (std::abs(to_z_axis_mapped.dot(to_y_axis_mapped)) > tolerance) throw std::runtime_error( "The z-axis and the y-axis are not orhogonal." );

        Eigen::Matrix<double,3,3> transformator;
        transformator << to_y_axis_mapped.cross(to_z_axis_mapped), to_y_axis_mapped, to_z_axis_mapped;

        return transformator;
    }

    Eigen::Matrix<double,3,3> buildRotator(const double &alpha, const double &beta, const double &gamma) {
        Eigen::Matrix<double,3,3> transformator;

        transformator = Eigen::AngleAxisd(alpha, Eigen::Matrix<double,3,1>::UnitZ())
                * Eigen::AngleAxisd(beta, Eigen::Matrix<double,3,1>::UnitY())
                * Eigen::AngleAxisd(gamma, Eigen::Matrix<double,3,1>::UnitZ());

        return transformator;
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
            removeRestrictedStates([=](const enumerated_state<T> &entry) -> bool { return this->checkIsQuantumstateValid(entry.state); } );
        }
        
        if (!range_n.empty() || !range_l.empty() || !range_j.empty() || !range_m.empty() ||
                energy_min != std::numeric_limits<double>::lowest() ||
                energy_max != std::numeric_limits<double>::max()) { // TODO also check for a new value of threshold_for_sqnorm
            
            ////////////////////////////////////////////////////////////////////
            /// Remove vectors with restricted energies or too small norm //////
            ////////////////////////////////////////////////////////////////////
            
            // Build transformator and remove vectors (if their energy is not allowed or the squared norm too small)
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
            removeRestrictedStates([=](const enumerated_state<T> &entry) -> bool { return sqnorm_list[entry.idx] > 0.05; } );
        }
        
        // Forget basis restrictions as they were applied now
        this->forgetRestrictions();
    }

    ////////////////////////////////////////////////////////////////////
    /// Method for serialization ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        (void)version;

        ar & cachedir;
        ar & energy_min & energy_max & range_n & range_l & range_j & range_m;
        ar & memory_saving & is_interaction_already_contained & is_new_hamiltonianmatrix_required;
        ar & states & coefficients & hamiltonianmatrix;
        ar & coefficients_unperturbed_cache & hamiltonianmatrix_unperturbed_cache;
    }
};

#endif
