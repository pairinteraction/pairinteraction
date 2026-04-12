// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/system/System.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/system/SystemPair.hpp"
#include "pairinteraction/utils/TaskControl.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/SparseCore>
#include <algorithm>
#include <complex>
#include <limits>
#include <memory>
#include <numeric>
#include <oneapi/tbb.h>
#include <optional>
#include <spdlog/spdlog.h>

namespace pairinteraction {
template <typename Derived>
System<Derived>::System(std::shared_ptr<const basis_t> basis)
    : basis(std::move(basis)),
      matrix(static_cast<Eigen::Index>(this->basis->get_number_of_states()),
             static_cast<Eigen::Index>(this->basis->get_number_of_states())) {}

template <typename Derived>
std::shared_ptr<const typename System<Derived>::basis_t> System<Derived>::get_basis() const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return basis;
}

template <typename Derived>
std::shared_ptr<const typename System<Derived>::basis_t> System<Derived>::get_eigenbasis() const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    if (!this->is_diagonal()) {
        throw std::runtime_error("The Hamiltonian has not been diagonalized yet.");
    }
    return basis;
}

template <typename Derived>
Eigen::VectorX<typename System<Derived>::real_t> System<Derived>::get_eigenenergies() const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    if (!this->is_diagonal()) {
        throw std::runtime_error("The Hamiltonian has not been diagonalized yet.");
    }
    return matrix.diagonal().real();
}

template <typename Derived>
const Eigen::SparseMatrix<typename System<Derived>::scalar_t, Eigen::RowMajor> &
System<Derived>::get_matrix() const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return matrix;
}

template <typename Derived>
const Transformation<typename System<Derived>::scalar_t> &
System<Derived>::get_transformation() const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return basis->get_transformation();
}

template <typename Derived>
Transformation<typename System<Derived>::scalar_t>
System<Derived>::get_rotator(real_t alpha, real_t beta, real_t gamma) const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return basis->get_rotator(alpha, beta, gamma);
}

template <typename Derived>
Sorting System<Derived>::get_sorter(const std::vector<TransformationType> &labels) const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }

    basis->perform_sorter_checks(labels);

    auto it = std::find(labels.begin(), labels.end(), TransformationType::SORT_BY_ENERGY);
    std::vector<TransformationType> before_energy(labels.begin(), it);
    bool contains_energy = (it != labels.end());
    std::vector<TransformationType> after_energy(contains_energy ? it + 1 : labels.end(),
                                                 labels.end());

    Sorting transformation;
    transformation.matrix.resize(matrix.rows());
    transformation.matrix.setIdentity();

    if (!before_energy.empty()) {
        basis->get_sorter_without_checks(before_energy, transformation);
    }

    if (contains_energy) {
        std::vector<real_t> energies_of_states;
        energies_of_states.reserve(matrix.rows());
        for (int i = 0; i < matrix.rows(); ++i) {
            energies_of_states.push_back(std::real(matrix.coeff(i, i)));
        }

        std::stable_sort(
            transformation.matrix.indices().data(),
            transformation.matrix.indices().data() + transformation.matrix.indices().size(),
            [&](int i, int j) { return energies_of_states[i] < energies_of_states[j]; });

        transformation.transformation_type.push_back(TransformationType::SORT_BY_ENERGY);
    }

    if (!after_energy.empty()) {
        basis->get_sorter_without_checks(after_energy, transformation);
    }

    if (labels != transformation.transformation_type) {
        throw std::invalid_argument("The states could not be sorted by all the requested labels.");
    }

    return transformation;
}

template <typename Derived>
std::vector<IndicesOfBlock>
System<Derived>::get_indices_of_blocks(const std::vector<TransformationType> &labels) const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }

    basis->perform_sorter_checks(labels);

    std::set<TransformationType> unique_labels(labels.begin(), labels.end());
    basis->perform_blocks_checks(unique_labels);

    auto it = unique_labels.find(TransformationType::SORT_BY_ENERGY);
    bool contains_energy = (it != unique_labels.end());
    if (contains_energy) {
        unique_labels.erase(it);
    }

    IndicesOfBlocksCreator blocks_creator({0, static_cast<size_t>(matrix.rows())});

    if (!unique_labels.empty()) {
        basis->get_indices_of_blocks_without_checks(unique_labels, blocks_creator);
    }

    if (contains_energy && matrix.rows() > 0) {
        scalar_t last_energy = std::real(matrix.coeff(0, 0));
        for (int i = 0; i < matrix.rows(); ++i) {
            if (std::real(matrix.coeff(i, i)) != last_energy) {
                blocks_creator.add(i);
                last_energy = std::real(matrix.coeff(i, i));
            }
        }
    }

    return blocks_creator.create();
}

template <typename Derived>
System<Derived> &System<Derived>::transform(const Transformation<scalar_t> &transformation) {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    if (matrix.cols() != 0) {
        matrix = transformation.matrix.adjoint() * matrix * transformation.matrix;
        basis = basis->transformed(transformation);
    }

    // A transformed system might have lost its block-diagonalizability if the
    // transformation was not a sorting
    blockdiagonalizing_labels.clear();

    return *this;
}

template <typename Derived>
System<Derived> &System<Derived>::transform(const Sorting &transformation) {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }

    if (matrix.cols() != 0) {
        matrix = matrix.twistedBy(transformation.matrix.inverse());
        basis = basis->transformed(transformation);
    }

    return *this;
}

template <typename Derived>
System<Derived> &System<Derived>::diagonalize(const DiagonalizerInterface<scalar_t> &diagonalizer,
                                              std::optional<real_t> min_eigenenergy,
                                              std::optional<real_t> max_eigenenergy, double rtol) {
    task_checkpoint("Preparing Hamiltonian...");

    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }

    if (this->is_diagonal()) {
        return *this;
    }

    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> eigenvectors;
    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> eigenenergies;

    // Sort the Hamiltonian according to the block structure
    if (!blockdiagonalizing_labels.empty()) {
        auto sorter = get_sorter(blockdiagonalizing_labels);
        matrix = matrix.twistedBy(sorter.matrix.inverse());
        basis = basis->transformed(sorter);
    }

    // Get the indices of the blocks
    auto blocks = get_indices_of_blocks(blockdiagonalizing_labels);

    assert((blockdiagonalizing_labels.empty() && blocks.size() == 1) ||
           !blockdiagonalizing_labels.empty());

    SPDLOG_DEBUG("Diagonalizing the Hamiltonian with {} blocks.", blocks.size());

    // Diagonalize the blocks in parallel
    std::vector<Eigen::VectorX<real_t>> eigenenergies_blocks(blocks.size());
    std::vector<Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>> eigenvectors_blocks(blocks.size());
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<size_t>(0, blocks.size()), [&](const auto &range) {
            for (size_t idx = range.begin(); idx != range.end(); ++idx) {
                task_checkpoint("Diagonalizing Hamiltonian blocks...");
                auto eigensys = min_eigenenergy.has_value() || max_eigenenergy.has_value()
                    ? diagonalizer.eigh(matrix.block(blocks[idx].start, blocks[idx].start,
                                                     blocks[idx].size(), blocks[idx].size()),
                                        min_eigenenergy, max_eigenenergy, rtol)
                    : diagonalizer.eigh(matrix.block(blocks[idx].start, blocks[idx].start,
                                                     blocks[idx].size(), blocks[idx].size()),
                                        rtol);
                eigenvectors_blocks[idx] = eigensys.eigenvectors;
                eigenenergies_blocks[idx] = eigensys.eigenvalues;
            }
        });

    // Get the number of non-zeros per row of the combined eigenvector matrix
    std::vector<Eigen::Index> non_zeros_per_inner_index;
    non_zeros_per_inner_index.reserve(matrix.rows());
    Eigen::Index num_rows = 0;
    Eigen::Index num_cols = 0;
    for (const auto &matrix : eigenvectors_blocks) {
        task_checkpoint("Collecting eigenvectors...");
        for (int i = 0; i < matrix.outerSize(); ++i) {
            non_zeros_per_inner_index.push_back(matrix.outerIndexPtr()[i + 1] -
                                                matrix.outerIndexPtr()[i]);
        }
        num_rows += matrix.rows();
        num_cols += matrix.cols();
    }

    assert(static_cast<size_t>(num_rows) == basis->get_number_of_kets());
    assert(static_cast<size_t>(num_cols) <= basis->get_number_of_states());

    eigenvectors.resize(num_rows, num_cols);
    eigenenergies.resize(num_cols, num_cols);

    if (num_cols > 0) {
        // Get the combined eigenvector matrix (in case of an restricted energy range, it is not
        // square)
        eigenvectors.reserve(non_zeros_per_inner_index);
        Eigen::Index offset_rows = 0;
        Eigen::Index offset_cols = 0;
        for (const auto &matrix : eigenvectors_blocks) {
            task_checkpoint("Combining eigenvectors...");
            for (Eigen::Index i = 0; i < matrix.outerSize(); ++i) {
                for (typename Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>::InnerIterator it(
                         matrix, i);
                     it; ++it) {
                    eigenvectors.insert(it.row() + offset_rows, it.col() + offset_cols) =
                        it.value();
                }
            }
            offset_rows += matrix.rows();
            offset_cols += matrix.cols();
        }
        eigenvectors.makeCompressed();

        assert(
            eigenvectors.nonZeros() ==
            std::accumulate(non_zeros_per_inner_index.begin(), non_zeros_per_inner_index.end(), 0));

        // Get the combined eigenenergy matrix
        eigenenergies.reserve(Eigen::VectorXi::Constant(num_cols, 1));
        Eigen::Index offset = 0;
        for (const auto &matrix : eigenenergies_blocks) {
            task_checkpoint("Combining eigenenergies...");
            for (int i = 0; i < matrix.size(); ++i) {
                eigenenergies.insert(i + offset, i + offset) = matrix(i);
            }
            offset += matrix.size();
        }
        eigenenergies.makeCompressed();

        // Fix phase ambiguity
        std::vector<scalar_t> map_col_to_max(num_cols, 0);
        for (int row = 0; row < eigenvectors.outerSize(); ++row) {
            task_checkpoint("Normalizing eigenvector phases...");
            for (typename Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>::InnerIterator it(
                     eigenvectors, row);
                 it; ++it) {
                if (std::abs(it.value()) > std::abs(map_col_to_max[it.col()])) {
                    map_col_to_max[it.col()] = it.value();
                }
            }
        }

        Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> phase_matrix;
        phase_matrix.resize(num_cols, num_cols);
        phase_matrix.reserve(Eigen::VectorXi::Constant(num_cols, 1));
        for (int i = 0; i < num_cols; ++i) {
            phase_matrix.insert(i, i) = std::abs(map_col_to_max[i]) / map_col_to_max[i];
        }
        phase_matrix.makeCompressed();

        task_checkpoint("Applying eigenvector phases...");
        eigenvectors = eigenvectors * phase_matrix;
    }

    // Store the diagonalized hamiltonian
    matrix = eigenenergies;
    basis = basis->transformed(eigenvectors);

    hamiltonian_is_diagonal = true;

    return *this;
}

template <typename Derived>
bool System<Derived>::is_diagonal() const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }

    if (!hamiltonian_is_diagonal) {
        real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();

        for (int row = 0; row < matrix.outerSize(); ++row) {
            for (typename Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>::InnerIterator it(matrix,
                                                                                           row);
                 it; ++it) {
                if (it.row() != it.col() && std::abs(it.value()) > numerical_precision) {
                    return false;
                }
            }
        }

        hamiltonian_is_diagonal = true;
    }

    return true;
}

// Explicit instantiation
template class System<SystemAtom<double>>;
template class System<SystemAtom<std::complex<double>>>;
template class System<SystemPair<double>>;
template class System<SystemPair<std::complex<double>>>;
} // namespace pairinteraction
