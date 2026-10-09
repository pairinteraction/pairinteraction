// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/interfaces/SorterBuilderInterface.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/SparseCore>
#include <memory>
#include <optional>
#include <set>
#include <vector>

namespace pairinteraction {
enum class SorterType : unsigned char;

template <typename Scalar>
class DiagonalizerInterface;

template <typename Derived>
class System : public SorterBuilderInterface {
public:
    using scalar_t = typename traits::CrtpTraits<Derived>::scalar_t;
    using real_t = typename traits::CrtpTraits<Derived>::real_t;
    using ketvec_t = typename traits::CrtpTraits<Derived>::ketvec_t;
    using basis_t = typename traits::CrtpTraits<Derived>::basis_t;

    System(std::shared_ptr<const basis_t> basis);

    std::shared_ptr<const basis_t> get_basis() const;
    std::shared_ptr<const basis_t> get_eigenbasis() const;
    Eigen::VectorX<real_t> get_eigenenergies() const;

    const Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &get_matrix() const;

    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>
    get_sorter(const std::vector<SorterType> &labels) const override;
    std::vector<IndicesOfBlock>
    get_indices_of_blocks(const std::vector<SorterType> &labels) const override;

    System<Derived> &
    transform(const Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &transformation);
    System<Derived> &
    transform(const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> &sorter);

    System<Derived> &diagonalize(const DiagonalizerInterface<scalar_t> &diagonalizer,
                                 std::optional<real_t> min_eigenenergy = {},
                                 std::optional<real_t> max_eigenenergy = {}, double rtol = 1e-6,
                                 bool sort_by_energy = true);
    bool is_diagonal() const;

protected:
    mutable std::shared_ptr<const basis_t> basis;
    mutable Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> matrix;
    mutable bool hamiltonian_requires_construction{true};
    mutable bool hamiltonian_is_diagonal{false};
    mutable std::vector<SorterType> blockdiagonalizing_labels;

    virtual void construct_hamiltonian() const = 0;
};
} // namespace pairinteraction
