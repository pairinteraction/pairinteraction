// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisAtom.hpp"

#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/ket/KetAtom.hpp"

#include <cassert>

namespace pairinteraction {
template <typename Scalar>
BasisAtom<Scalar>::BasisAtom(Private /*unused*/, ketvec_t &&kets, std::string &&canonical_basis_id,
                             Database &database)
    : Basis<BasisAtom<Scalar>>(std::move(kets)), canonical_basis_id(std::move(canonical_basis_id)),
      database(database) {
    for (size_t i = 0; i < this->kets.size(); ++i) {
        ket_id_to_ket_index[this->kets[i]->get_id_in_database()] = i;
    }
}

template <typename Scalar>
Database &BasisAtom<Scalar>::get_database() const {
    return database;
}

template <typename Scalar>
const std::string &BasisAtom<Scalar>::get_species() const {
    return this->kets[0]->get_species();
}

template <typename Scalar>
int BasisAtom<Scalar>::get_ket_index_from_id(size_t ket_id) const {
    if (!ket_id_to_ket_index.contains(ket_id)) {
        return -1;
    }
    return ket_id_to_ket_index.at(ket_id);
}

template <typename Scalar>
const std::string &BasisAtom<Scalar>::get_canonical_basis_id() const {
    return canonical_basis_id;
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
BasisAtom<Scalar>::get_matrix_elements(std::shared_ptr<const Type> other, OperatorType type,
                                       int q) const {
    auto matrix_elements = this->get_database().get_matrix_elements_in_canonical_basis(
        this->shared_from_this(), other, type, q);
    matrix_elements =
        other->get_coefficients().adjoint() * matrix_elements * this->get_coefficients();

    assert(static_cast<size_t>(matrix_elements.rows()) == other->get_number_of_states());
    assert(static_cast<size_t>(matrix_elements.cols()) == this->get_number_of_states());

    return matrix_elements;
}

// Explicit instantiations
template class BasisAtom<double>;
template class BasisAtom<std::complex<double>>;
} // namespace pairinteraction
