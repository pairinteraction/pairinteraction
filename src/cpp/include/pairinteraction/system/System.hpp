// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/SparseCore>
#include <memory>
#include <optional>
#include <set>
#include <vector>

namespace pairinteraction {
enum class StateLabel : unsigned char;

template <typename Scalar>
class DiagonalizerInterface;

template <typename Derived>
class System {
public:
    using scalar_t = typename traits::CrtpTraits<Derived>::scalar_t;
    using real_t = typename traits::CrtpTraits<Derived>::real_t;
    using ketvec_t = typename traits::CrtpTraits<Derived>::ketvec_t;
    using basis_t = typename traits::CrtpTraits<Derived>::basis_t;

    std::shared_ptr<const basis_t> get_basis() const;
    std::shared_ptr<const basis_t> get_eigenbasis() const;
    Eigen::VectorX<real_t> get_eigenenergies() const;

    const Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &get_matrix() const;

    System<Derived> &diagonalize(const DiagonalizerInterface<scalar_t> &diagonalizer,
                                 std::optional<real_t> min_eigenenergy = {},
                                 std::optional<real_t> max_eigenenergy = {}, double rtol = 1e-6,
                                 bool sort_by_energy = true);
    bool is_diagonal() const;

protected:
    System(std::shared_ptr<const basis_t> basis);

    mutable std::shared_ptr<const basis_t> basis;
    mutable Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> matrix;
    mutable bool hamiltonian_requires_construction{true};
    mutable bool hamiltonian_is_diagonal{false};
    // The labels of the states that are conserved by the Hamiltonian, so that it is block-diagonal
    // with respect to them
    mutable std::set<StateLabel> blockdiagonalizing_labels;

private:
    void construct_hamiltonian() const;
    void prune_couplings_between_blocks() const;
    void sort(const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> &sorter);
};
} // namespace pairinteraction
