// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/system/System.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/enums/SorterType.hpp"
#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/system/SystemPair.hpp"
#include "pairinteraction/utils/TaskControl.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/SparseCore>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fmt/format.h>
#include <limits>
#include <memory>
#include <numeric>
#include <oneapi/tbb.h>
#include <optional>
#include <spdlog/spdlog.h>
#include <stdexcept>
#include <vector>

namespace pairinteraction {
template <typename Derived>
System<Derived>::System(std::shared_ptr<const basis_t> basis)
    : basis(std::move(basis)),
      matrix(static_cast<Eigen::Index>(this->basis->get_number_of_states()),
             static_cast<Eigen::Index>(this->basis->get_number_of_states())) {
    // The unperturbed Hamiltonian is constructed from the energies of the kets
    for (const auto &ket : this->basis->get_kets()) {
        if (static_cast<real_t>(ket->get_energy()) == std::numeric_limits<real_t>::max()) {
            throw std::invalid_argument("The kets of the basis do not have a well-defined energy.");
        }
    }
}

template <typename Derived>
void System<Derived>::ensure_hamiltonian_constructed() const {
    if (!hamiltonian_requires_construction) {
        return;
    }
    construct_hamiltonian();
    hamiltonian_requires_construction = false;
}

template <typename Derived>
void System<Derived>::prune_couplings_between_blocks() const {
    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();

    // Get the block-diagonalizing labels of the states. The quantum numbers are (half-)integers,
    // thus doubling them yields an exact integer.
    std::vector<std::vector<int64_t>> labels_of_states;
    labels_of_states.reserve(blockdiagonalizing_labels.size());
    for (const auto &label : blockdiagonalizing_labels) {
        const auto &name = basis_t::sorter_type_to_quantum_number_name.at(label);
        auto &values = labels_of_states.emplace_back(basis->get_number_of_states());
        for (size_t state_index = 0; state_index < values.size(); ++state_index) {
            values[state_index] = std::llround(2 * basis->get_quantum_number(name, state_index));
        }
    }

    // Couplings between states with different labels vanish analytically. However, they can be
    // stored with rounding errors, e.g., in a symmetrized pair basis, whose couplings between
    // different symmetry sectors do not cancel exactly on all platforms (e.g., because of fused
    // multiply-adds). We drop these couplings so that they neither take up memory nor break the
    // block structure. As a consistency check that no physical coupling is dropped, we throw if a
    // coupling is larger than the rounding errors. Explicitly stored zeros are dropped as well.
    matrix.makeCompressed();
    const real_t tolerance =
        matrix.nonZeros() == 0 ? 0 : numerical_precision * matrix.coeffs().cwiseAbs().maxCoeff();
    matrix.prune([&](Eigen::Index row, Eigen::Index col, const scalar_t &value) {
        if (value == scalar_t{0}) {
            return false;
        }
        bool is_within_block =
            std::all_of(labels_of_states.begin(), labels_of_states.end(),
                        [&](const auto &values) { return values[row] == values[col]; });
        if (is_within_block) {
            return true;
        }
        if (std::abs(value) > tolerance) {
            throw std::runtime_error(fmt::format(
                "The Hamiltonian couples the states {} and {} with the coupling {:.3e} (tolerance "
                "{:.3e}), although they belong to different blocks. This should not happen.",
                row, col, std::abs(value), tolerance));
        }
        return false;
    });
}

template <typename Derived>
std::shared_ptr<const typename System<Derived>::basis_t> System<Derived>::get_basis() const {
    ensure_hamiltonian_constructed();
    return basis;
}

template <typename Derived>
std::shared_ptr<const typename System<Derived>::basis_t> System<Derived>::get_eigenbasis() const {
    ensure_hamiltonian_constructed();
    if (!this->is_diagonal()) {
        throw std::runtime_error("The Hamiltonian has not been diagonalized yet.");
    }
    return basis;
}

template <typename Derived>
Eigen::VectorX<typename System<Derived>::real_t> System<Derived>::get_eigenenergies() const {
    ensure_hamiltonian_constructed();
    if (!this->is_diagonal()) {
        throw std::runtime_error("The Hamiltonian has not been diagonalized yet.");
    }
    return matrix.diagonal().real();
}

template <typename Derived>
const Eigen::SparseMatrix<typename System<Derived>::scalar_t, Eigen::RowMajor> &
System<Derived>::get_matrix() const {
    ensure_hamiltonian_constructed();
    return matrix;
}

template <typename Derived>
Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>
System<Derived>::get_sorter(const std::vector<SorterType> &labels) const {
    // The first sorter type in labels is supposed to be the primary sorting criterion.
    ensure_hamiltonian_constructed();

    // Drop repeated sorter types. This also ensures that energy occurs at most once.
    std::vector<SorterType> unique_labels;
    unique_labels.reserve(labels.size());
    for (const auto &label : labels) {
        if (std::find(unique_labels.begin(), unique_labels.end(), label) == unique_labels.end()) {
            unique_labels.push_back(label);
        }
    }

    auto it = std::find(unique_labels.begin(), unique_labels.end(), SorterType::ENERGY);
    std::vector<SorterType> before_energy(unique_labels.begin(), it);
    bool contains_energy = (it != unique_labels.end());
    std::vector<SorterType> after_energy(contains_energy ? it + 1 : unique_labels.end(),
                                         unique_labels.end());

    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> sorter(matrix.rows());
    sorter.setIdentity();

    if (!after_energy.empty()) {
        basis->get_sorter_without_checks(after_energy, sorter);
    }

    if (contains_energy) {
        std::vector<real_t> energies_of_states;
        energies_of_states.reserve(matrix.rows());
        for (int i = 0; i < matrix.rows(); ++i) {
            energies_of_states.push_back(std::real(matrix.coeff(i, i)));
        }

        std::stable_sort(
            sorter.indices().data(), sorter.indices().data() + sorter.indices().size(),
            [&](int i, int j) { return energies_of_states[i] < energies_of_states[j]; });
    }

    if (!before_energy.empty()) {
        basis->get_sorter_without_checks(before_energy, sorter);
    }

    return sorter;
}

template <typename Derived>
std::vector<IndicesOfBlock>
System<Derived>::get_indices_of_blocks(const std::vector<SorterType> &labels) const {
    ensure_hamiltonian_constructed();

    std::set<SorterType> unique_labels(labels.begin(), labels.end());
    basis->perform_blocks_checks(unique_labels);

    IndicesOfBlocksCreator blocks_creator({0, static_cast<size_t>(matrix.rows())});

    if (!unique_labels.empty()) {
        basis->get_indices_of_blocks_without_checks(unique_labels, blocks_creator);
    }

    return blocks_creator.create();
}

template <typename Derived>
void System<Derived>::sort(const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> &sorter) {
    set_task_status("Sorting Hamiltonian...");

    if (!hamiltonian_requires_construction) {
        matrix = matrix.twistedBy(sorter.inverse());
    }
    basis = basis->get_sorted_copy(sorter);
}

template <typename Derived>
System<Derived> &System<Derived>::diagonalize(const DiagonalizerInterface<scalar_t> &diagonalizer,
                                              std::optional<real_t> min_eigenenergy,
                                              std::optional<real_t> max_eigenenergy, double rtol,
                                              bool sort_by_energy) {
    set_task_status("Preparing Hamiltonian...");

    ensure_hamiltonian_constructed();

    if (this->is_diagonal()) {
        if (sort_by_energy && !basis->is_sorted_by_energy()) {
            sort(get_sorter({SorterType::ENERGY}));
        }
        return *this;
    }

    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> eigenvectors;
    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> eigenenergies;

    // Sort the Hamiltonian according to the block structure
    if (!blockdiagonalizing_labels.empty()) {
        sort(get_sorter(blockdiagonalizing_labels));
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
                set_task_status("Diagonalizing Hamiltonian blocks...");
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
        set_task_status("Collecting eigenvectors...");
        for (int i = 0; i < matrix.outerSize(); ++i) {
            non_zeros_per_inner_index.push_back(matrix.outerIndexPtr()[i + 1] -
                                                matrix.outerIndexPtr()[i]);
        }
        num_rows += matrix.rows();
        num_cols += matrix.cols();
    }

    assert(static_cast<size_t>(num_rows) == basis->get_number_of_states());
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
            set_task_status("Combining eigenvectors...");
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
            set_task_status("Combining eigenenergies...");
            for (int i = 0; i < matrix.size(); ++i) {
                eigenenergies.insert(i + offset, i + offset) = matrix(i);
            }
            offset += matrix.size();
        }
        eigenenergies.makeCompressed();

        // Fix phase ambiguity
        std::vector<scalar_t> map_col_to_max(num_cols, 0);
        for (int row = 0; row < eigenvectors.outerSize(); ++row) {
            set_task_status("Normalizing eigenvector phases...");
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

        set_task_status("Applying eigenvector phases...");
        eigenvectors = eigenvectors * phase_matrix;
    }

    // Assign the eigenenergies to the eigenstates
    std::vector<real_t> energy_of_eigenstates;
    energy_of_eigenstates.reserve(static_cast<size_t>(num_cols));
    for (const auto &eigenenergies_block : eigenenergies_blocks) {
        energy_of_eigenstates.insert(energy_of_eigenstates.end(), eigenenergies_block.data(),
                                     eigenenergies_block.data() + eigenenergies_block.size());
    }

    // Assign the quantum numbers to the eigenstates, only the quantum numbers by which the
    // Hamiltonian is block-diagonal are well-defined and shared by all eigenstates of a block
    typename basis_t::quantum_numbers_t quantum_numbers_of_eigenstates;
    for (const auto &[label, name] : basis_t::sorter_type_to_quantum_number_name) {
        auto &quantum_numbers = quantum_numbers_of_eigenstates[name];
        quantum_numbers.reserve(static_cast<size_t>(num_cols));
        const bool is_blockdiagonalizing_label =
            std::find(blockdiagonalizing_labels.begin(), blockdiagonalizing_labels.end(), label) !=
            blockdiagonalizing_labels.end();
        for (size_t idx = 0; idx < blocks.size(); ++idx) {
            quantum_numbers.insert(
                quantum_numbers.end(), static_cast<size_t>(eigenvectors_blocks[idx].cols()),
                is_blockdiagonalizing_label ? basis->get_quantum_number(name, blocks[idx].start)
                                            : std::numeric_limits<real_t>::max());
        }
    }

    // Store the diagonalized hamiltonian and the eigenbasis
    set_task_status("Constructing eigenbasis...");
    matrix = eigenenergies;
    basis = basis->get_modified_copy(basis->get_coefficients() * eigenvectors,
                                     std::move(energy_of_eigenstates),
                                     std::move(quantum_numbers_of_eigenstates));

    hamiltonian_is_diagonal = true;
    if (sort_by_energy) {
        sort(get_sorter({SorterType::ENERGY}));
    }

    return *this;
}

template <typename Derived>
bool System<Derived>::is_diagonal() const {
    ensure_hamiltonian_constructed();

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

        // The states of the basis are eigenstates, so assign them the diagonal entries of the
        // Hamiltonian as energies
        Eigen::VectorX<real_t> diagonal = matrix.diagonal().real();
        basis = basis->get_modified_copy(
            basis->get_coefficients(),
            std::vector<real_t>(diagonal.data(), diagonal.data() + diagonal.size()),
            basis->get_quantum_numbers_of_states());

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
