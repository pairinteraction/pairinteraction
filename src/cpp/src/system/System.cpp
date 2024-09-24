#include "pairinteraction/system/System.hpp"

#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"

#include <Eigen/SparseCore>
#include <complex>
#include <limits>
#include <memory>
#include <oneapi/tbb.h>

namespace pairinteraction {
template <typename Derived>
System<Derived>::System(std::shared_ptr<const basis_t> basis)
    : hamiltonian(std::make_unique<typename System<Derived>::operator_t>(basis)) {}

template <typename Derived>
System<Derived>::System(const System &other)
    : hamiltonian(std::make_unique<typename System<Derived>::operator_t>(*other.hamiltonian)) {}

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
IndicesOfBlocks
System<Derived>::get_indices_of_blocks(const std::vector<TransformationType> &labels) const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return hamiltonian->get_indices_of_blocks(labels);
}

template <typename Derived>
Derived System<Derived>::transformed(const Transformation<scalar_t> &transformation) const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    Derived transformed = derived();
    transformed.hamiltonian =
        std::make_unique<operator_t>(hamiltonian->transformed(transformation));
    return transformed;
}

template <typename Derived>
Derived System<Derived>::transformed(const Sorting &transformation) const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    Derived transformed = derived();
    transformed.hamiltonian =
        std::make_unique<operator_t>(hamiltonian->transformed(transformation));
    return transformed;
}

template <typename Derived>
System<Derived> &System<Derived>::diagonalize(const DiagonalizerInterface<scalar_t> &diagonalizer,
                                              int precision,
                                              const Range<real_t> &eigenvalue_range) {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }

    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> eigenvectors;

    // Figure out whether the Hamiltonian can be block-diagonalized
    bool is_blockdiagonizable = !blockdiagonalizing_labels.empty();
    if (is_blockdiagonizable) {
        std::vector<TransformationType> labels{blockdiagonalizing_labels.begin(),
                                               blockdiagonalizing_labels.end()};
        for (const auto &label : hamiltonian->get_transformation().transformation_type) {
            if (!utils::is_sorting(label) && label != TransformationType::IDENTITY) {
                is_blockdiagonizable = false;
                break;
            }
        }
    }

    // Block-diagonalize the Hamiltonian if possible
    if (is_blockdiagonizable) {
        std::vector<TransformationType> labels{blockdiagonalizing_labels.begin(),
                                               blockdiagonalizing_labels.end()};

        // Sort the Hamiltonian according to the block structure
        auto sorter = hamiltonian->get_sorter(labels);
        hamiltonian = std::make_unique<operator_t>(hamiltonian->transformed(sorter));

        // Get the indices of the blocks
        auto blocks =
            hamiltonian->get_indices_of_blocks(labels)
                .get(); // TODO: get_indices_of_blocks should make use of a IndicesOfBlocksCreator

        // Split the Hamiltonian matrix into blocks
        std::vector<Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>> matrices;
        matrices.reserve(blocks.size());
        for (const auto &block : blocks) {
            matrices.push_back(hamiltonian->get_matrix().block(block.start, block.start,
                                                               block.size(), block.size()));
        }

        // Diagonalize the blocks in parallel
        oneapi::tbb::parallel_for(oneapi::tbb::blocked_range(matrices.begin(), matrices.end()),
                                  [&](const auto &range) {
                                      for (auto &matrix : range) {
                                          auto eigensys = eigenvalue_range.is_finite()
                                              ? diagonalizer.eigh(matrix, eigenvalue_range.min(),
                                                                  eigenvalue_range.max(), precision)
                                              : diagonalizer.eigh(matrix, precision);
                                          matrix = eigensys.eigenvectors;
                                      }
                                  });

        // Get the combined eigenvector matrix
        eigenvectors.resize(hamiltonian->get_matrix().rows(), hamiltonian->get_matrix().cols());

        // Get the number of non-zeros per row
        std::vector<Eigen::Index> non_zeros_per_inner_index;
        non_zeros_per_inner_index.reserve(hamiltonian->get_matrix().rows());
        for (const auto &matrix : matrices) {
            for (int i = 0; i < matrix.outerSize(); ++i) {
                non_zeros_per_inner_index.push_back(matrix.outerIndexPtr()[i + 1] -
                                                    matrix.outerIndexPtr()[i]);
            }
        }

        // Allocate memory for the combined eigenvector matrix
        eigenvectors.reserve(non_zeros_per_inner_index);

        // Fill the combined eigenvector matrix
        Eigen::Index offset = 0;
        for (const auto &matrix : matrices) {
            for (Eigen::Index i = 0; i < matrix.outerSize(); ++i) {
                for (typename Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>::InnerIterator it(
                         matrix, i);
                     it; ++it) {
                    eigenvectors.insert(it.row() + offset, it.col() + offset) = it.value();
                }
            }
            offset += matrix.rows();
        }

    } else {
        // Diagonalize the full Hamiltonian at once
        auto eigensys = eigenvalue_range.is_finite()
            ? diagonalizer.eigh(hamiltonian->get_matrix(), eigenvalue_range.min(),
                                eigenvalue_range.max(), precision)
            : diagonalizer.eigh(hamiltonian->get_matrix(), precision);

        // Get the eigenvector matrix
        eigenvectors = std::move(eigensys.eigenvectors);
    }

    // Store the diagonalized hamiltonian (possible future optimization: use
    // eigensys.eigenvalues directly instead of transforming the hamiltonian with the
    // eigenvectors, get rid of values smaller than the precision)
    hamiltonian = std::make_unique<operator_t>(hamiltonian->transformed(eigenvectors));

    return *this;
}

// Explicit instantiation
template class System<SystemAtom<float>>;
template class System<SystemAtom<double>>;
template class System<SystemAtom<std::complex<float>>>;
template class System<SystemAtom<std::complex<double>>>;
} // namespace pairinteraction
