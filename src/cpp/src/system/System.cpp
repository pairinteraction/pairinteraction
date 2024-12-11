#include "pairinteraction/system/System.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisCombined.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"
#include "pairinteraction/operator/OperatorCombined.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/system/SystemCombined.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/SparseCore>
#include <complex>
#include <limits>
#include <memory>
#include <numeric>
#include <oneapi/tbb.h>
#include <spdlog/spdlog.h>

namespace pairinteraction {
template <typename Derived>
System<Derived>::System(std::shared_ptr<const basis_t> basis)
    : hamiltonian(std::make_unique<typename System<Derived>::operator_t>(std::move(basis))) {}

template <typename Derived>
System<Derived>::System(const System &other)
    : hamiltonian(std::make_unique<typename System<Derived>::operator_t>(*other.hamiltonian)),
      hamiltonian_requires_construction(other.hamiltonian_requires_construction),
      hamiltonian_is_diagonal(other.hamiltonian_is_diagonal),
      blockdiagonalizing_labels(other.blockdiagonalizing_labels) {}

template <typename Derived>
System<Derived>::System(System &&other) noexcept
    : hamiltonian(std::move(other.hamiltonian)),
      hamiltonian_requires_construction(other.hamiltonian_requires_construction),
      hamiltonian_is_diagonal(other.hamiltonian_is_diagonal),
      blockdiagonalizing_labels(std::move(other.blockdiagonalizing_labels)) {}

template <typename Derived>
System<Derived> &System<Derived>::operator=(const System &other) {
    if (this != &other) {
        hamiltonian = std::make_unique<operator_t>(*other.hamiltonian);
        hamiltonian_requires_construction = other.hamiltonian_requires_construction;
        hamiltonian_is_diagonal = other.hamiltonian_is_diagonal,
        blockdiagonalizing_labels = other.blockdiagonalizing_labels;
    }
    return *this;
}

template <typename Derived>
System<Derived> &System<Derived>::operator=(System &&other) noexcept {
    if (this != &other) {
        hamiltonian = std::move(other.hamiltonian);
        hamiltonian_requires_construction = other.hamiltonian_requires_construction;
        hamiltonian_is_diagonal = other.hamiltonian_is_diagonal,
        blockdiagonalizing_labels = std::move(other.blockdiagonalizing_labels);
    }
    return *this;
}

template <typename Derived>
System<Derived>::~System() = default;

template <typename Derived>
const Derived &System<Derived>::derived() const {
    return static_cast<const Derived &>(*this);
}

template <typename Derived>
std::shared_ptr<const typename System<Derived>::basis_t> System<Derived>::get_basis() const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return hamiltonian->get_basis();
}

template <typename Derived>
std::shared_ptr<const typename System<Derived>::basis_t> System<Derived>::get_eigenbasis() const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    if (!hamiltonian_is_diagonal) {
        throw std::runtime_error("The Hamiltonian has not been diagonalized yet.");
    }
    return hamiltonian->get_basis();
}

template <typename Derived>
Eigen::VectorX<typename System<Derived>::real_t> System<Derived>::get_eigenvalues() const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    if (!hamiltonian_is_diagonal) {
        throw std::runtime_error("The Hamiltonian has not been diagonalized yet.");
    }
    return hamiltonian->get_matrix().diagonal().real();
}

template <typename Derived>
const Eigen::SparseMatrix<typename System<Derived>::scalar_t, Eigen::RowMajor> &
System<Derived>::get_matrix() const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return hamiltonian->get_matrix();
}

template <typename Derived>
const Transformation<typename System<Derived>::scalar_t> &
System<Derived>::get_transformation() const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return hamiltonian->get_transformation();
}

template <typename Derived>
Transformation<typename System<Derived>::scalar_t>
System<Derived>::get_rotator(real_t alpha, real_t beta, real_t gamma) const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return hamiltonian->get_rotator(alpha, beta, gamma);
}

template <typename Derived>
Sorting System<Derived>::get_sorter(const std::vector<TransformationType> &labels) const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return hamiltonian->get_sorter(labels);
}

template <typename Derived>
std::vector<IndicesOfBlock>
System<Derived>::get_indices_of_blocks(const std::vector<TransformationType> &labels) const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return hamiltonian->get_indices_of_blocks(labels);
}

template <typename Derived>
System<Derived> &System<Derived>::transform(const Transformation<scalar_t> &transformation) {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    hamiltonian = std::make_unique<operator_t>(hamiltonian->transformed(transformation));

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
    hamiltonian = std::make_unique<operator_t>(hamiltonian->transformed(transformation));

    return *this;
}

template <typename Derived>
System<Derived> &System<Derived>::diagonalize(const DiagonalizerInterface<scalar_t> &diagonalizer,
                                              int precision,
                                              const Range<real_t> &eigenvalue_range) {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }

    if (hamiltonian_is_diagonal) {
        return *this;
    }

    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> eigenvectors;
    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> eigenvalues;

    // Sort the Hamiltonian according to the block structure
    if (!blockdiagonalizing_labels.empty()) {
        auto sorter = hamiltonian->get_sorter(blockdiagonalizing_labels);
        hamiltonian = std::make_unique<operator_t>(hamiltonian->transformed(sorter));
    }

    // Get the indices of the blocks
    auto blocks = hamiltonian->get_indices_of_blocks(blockdiagonalizing_labels);

    assert((blockdiagonalizing_labels.empty() && blocks.size() == 1) ||
           !blockdiagonalizing_labels.empty());

    SPDLOG_DEBUG("Diagonalizing the Hamiltonian with {} blocks.", blocks.size());

    // Diagonalize the blocks in parallel
    std::vector<Eigen::VectorX<real_t>> eigenvalues_blocks(blocks.size());
    std::vector<Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>> eigenvectors_blocks(blocks.size());
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<size_t>(0, blocks.size()), [&](const auto &range) {
            for (size_t idx = range.begin(); idx != range.end(); ++idx) {
                auto eigensys = eigenvalue_range.is_finite()
                    ? diagonalizer.eigh(
                          hamiltonian->get_matrix().block(blocks[idx].start, blocks[idx].start,
                                                          blocks[idx].size(), blocks[idx].size()),
                          eigenvalue_range.min(), eigenvalue_range.max(), precision)
                    : diagonalizer.eigh(
                          hamiltonian->get_matrix().block(blocks[idx].start, blocks[idx].start,
                                                          blocks[idx].size(), blocks[idx].size()),
                          precision);
                eigenvectors_blocks[idx] = eigensys.eigenvectors;
                eigenvalues_blocks[idx] = eigensys.eigenvalues;
            }
        });

    // Get the number of non-zeros per row of the combined eigenvector matrix
    std::vector<Eigen::Index> non_zeros_per_inner_index;
    non_zeros_per_inner_index.reserve(hamiltonian->get_matrix().rows());
    Eigen::Index num_rows = 0;
    Eigen::Index num_cols = 0;
    for (const auto &matrix : eigenvectors_blocks) {
        for (int i = 0; i < matrix.outerSize(); ++i) {
            non_zeros_per_inner_index.push_back(matrix.outerIndexPtr()[i + 1] -
                                                matrix.outerIndexPtr()[i]);
        }
        num_rows += matrix.rows();
        num_cols += matrix.cols();
    }

    // Get the combined eigenvector matrix (in case of an restricted energy range, it is not
    // square)
    eigenvectors.resize(num_rows, num_cols);
    eigenvectors.reserve(non_zeros_per_inner_index);
    Eigen::Index offset_rows = 0;
    Eigen::Index offset_cols = 0;
    for (const auto &matrix : eigenvectors_blocks) {
        for (Eigen::Index i = 0; i < matrix.outerSize(); ++i) {
            for (typename Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>::InnerIterator it(matrix,
                                                                                           i);
                 it; ++it) {
                eigenvectors.insert(it.row() + offset_rows, it.col() + offset_cols) = it.value();
            }
        }
        offset_rows += matrix.rows();
        offset_cols += matrix.cols();
    }
    eigenvectors.makeCompressed();

    assert(eigenvectors.nonZeros() ==
           std::accumulate(non_zeros_per_inner_index.begin(), non_zeros_per_inner_index.end(), 0));

    // Get the combined eigenvalue matrix
    eigenvalues.resize(num_cols, num_cols);
    eigenvalues.reserve(Eigen::VectorXi::Constant(num_cols, 1));
    Eigen::Index offset = 0;
    for (const auto &matrix : eigenvalues_blocks) {
        for (int i = 0; i < matrix.size(); ++i) {
            eigenvalues.insert(i + offset, i + offset) = matrix(i);
        }
        offset += matrix.size();
    }
    eigenvalues.makeCompressed();

    // Fix phase ambiguity
    std::vector<scalar_t> map_col_to_max(num_cols, 0);
    for (int row = 0; row < eigenvectors.outerSize(); ++row) {
        for (typename Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>::InnerIterator it(eigenvectors,
                                                                                       row);
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

    eigenvectors = eigenvectors * phase_matrix;

    // Store the diagonalized hamiltonian
    hamiltonian->get_matrix() = eigenvalues;
    hamiltonian->get_basis() = hamiltonian->get_basis()->transformed(eigenvectors);

    hamiltonian_is_diagonal = true;

    return *this;
}

template <typename Derived>
bool System<Derived>::is_diagonal() const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return hamiltonian_is_diagonal;
}

// Explicit instantiation
template class System<SystemAtom<float>>;
template class System<SystemAtom<double>>;
template class System<SystemAtom<std::complex<float>>>;
template class System<SystemAtom<std::complex<double>>>;
template class System<SystemCombined<float>>;
template class System<SystemCombined<double>>;
template class System<SystemCombined<std::complex<float>>>;
template class System<SystemCombined<std::complex<double>>>;
} // namespace pairinteraction
