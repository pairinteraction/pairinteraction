// SPDX-FileCopyrightText: 2026 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/basis/Basis.hpp"

#include <Eigen/SparseCore>
#include <memory>

namespace pairinteraction::utils {
template <typename BasisType>
Eigen::SparseMatrix<typename BasisType::scalar_t, Eigen::RowMajor>
get_energies_in_canonical_basis(const std::shared_ptr<const BasisType> &basis) {
    using scalar_t = typename BasisType::scalar_t;

    const auto num_kets = static_cast<Eigen::Index>(basis->get_number_of_kets());

    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> matrix(num_kets, num_kets);
    matrix.reserve(Eigen::VectorXi::Constant(num_kets, 1));
    for (Eigen::Index idx = 0; idx < num_kets; ++idx) {
        matrix.insert(idx, idx) = basis->get_ket(static_cast<size_t>(idx))->get_energy();
    }
    matrix.makeCompressed();

    return matrix;
}
} // namespace pairinteraction::utils
