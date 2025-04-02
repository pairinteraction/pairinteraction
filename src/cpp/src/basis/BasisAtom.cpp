// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisAtom.hpp"

#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/ket/KetAtom.hpp"

namespace pairinteraction {
template <typename Scalar>
BasisAtom<Scalar>::BasisAtom(Private /*unused*/, ketvec_t &&kets, std::string &&id_of_kets,
                             Database &database)
    : Basis<BasisAtom<Scalar>>(std::move(kets)), id_of_kets(std::move(id_of_kets)),
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
    if (ket_id_to_ket_index.count(ket_id) == 0) {
        return -1;
    }
    return ket_id_to_ket_index.at(ket_id);
}

template <typename Scalar>
const std::string &BasisAtom<Scalar>::get_id_of_kets() const {
    return id_of_kets;
}

template <typename Scalar>
Eigen::VectorX<Scalar> BasisAtom<Scalar>::get_matrix_elements(std::shared_ptr<const ket_t> ket,
                                                              OperatorType type, int q) const {
    auto final = this->get_canonical_state_from_ket(ket);
    auto matrix_elements =
        this->get_database().get_matrix_elements(this->shared_from_this(), final, type, q);

    assert(static_cast<size_t>(matrix_elements.rows()) == 1);
    assert(static_cast<size_t>(matrix_elements.cols()) == this->get_number_of_states());

    return matrix_elements.row(0);
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
BasisAtom<Scalar>::get_matrix_elements(std::shared_ptr<const Type> other, OperatorType type,
                                       int q) const {
    auto matrix_elements =
        this->get_database().get_matrix_elements(this->shared_from_this(), other, type, q);

    assert(static_cast<size_t>(matrix_elements.rows()) == other->get_number_of_states());
    assert(static_cast<size_t>(matrix_elements.cols()) == this->get_number_of_states());

    return matrix_elements;
}

// Explicit instantiations
template class BasisAtom<double>;
template class BasisAtom<std::complex<double>>;
} // namespace pairinteraction
