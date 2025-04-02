// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/operator/OperatorAtom.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/enums/OperatorType.hpp"

#include <limits>

namespace pairinteraction {
template <typename Scalar>
OperatorAtom<Scalar>::OperatorAtom(std::shared_ptr<const basis_t> basis)
    : Operator<OperatorAtom<Scalar>>(std::move(basis)) {}

template <typename Scalar>
OperatorAtom<Scalar>::OperatorAtom(std::shared_ptr<const basis_t> basis, OperatorType type, int q)
    : Operator<OperatorAtom<Scalar>>(std::move(basis)) {
    if (type == OperatorType::ENERGY) {
        this->initialize_as_energy_operator();
    } else {
        this->initialize_from_matrix(
            this->basis->get_database().get_matrix_elements(this->basis, this->basis, type, q));
    }
}

template <typename Scalar>
OperatorAtom<Scalar>::OperatorAtom(std::shared_ptr<const basis_t> basis,
                                   Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &&matrix)
    : Operator<OperatorAtom<Scalar>>(std::move(basis)) {
    this->matrix = std::move(matrix);
}

// Explicit instantiations
template class OperatorAtom<double>;
template class OperatorAtom<std::complex<double>>;
} // namespace pairinteraction
