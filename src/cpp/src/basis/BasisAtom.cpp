// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisAtom.hpp"

#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/ket/KetAtom.hpp"

#include <Eigen/Sparse>
#include <cassert>
#include <vector>

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
    if (!ket_id_to_ket_index.contains(ket_id)) {
        return -1;
    }
    return ket_id_to_ket_index.at(ket_id);
}

template <typename Scalar>
const std::string &BasisAtom<Scalar>::get_id_of_kets() const {
    return id_of_kets;
}

template <typename Scalar>
Eigen::VectorX<Scalar> BasisAtom<Scalar>::get_amplitudes(std::shared_ptr<const ket_t> ket) const {
    if (get_ket_index_from_id(ket->get_id_in_database()) < 0) {
        return Eigen::VectorX<Scalar>::Zero(
            static_cast<Eigen::Index>(this->get_number_of_states()));
    }
    return this->Basis<BasisAtom<Scalar>>::get_amplitudes(ket);
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
BasisAtom<Scalar>::get_amplitudes(std::shared_ptr<const Type> other) const {
    // Fast path: identical ket ordering means K is the identity matrix
    if (this->get_id_of_kets() == other->get_id_of_kets()) {
        return other->get_coefficients().adjoint() * this->get_coefficients();
    }

    // General path: build ket-matching matrix K of shape (n_kets_other, n_kets_this)
    // K[j, i] = 1 if ket j of other matches ket i of this (by database ID)
    std::vector<Eigen::Triplet<Scalar>> triplets;
    const auto &other_kets = other->get_kets();
    for (size_t j = 0; j < other_kets.size(); ++j) {
        int i = get_ket_index_from_id(other_kets[j]->get_id_in_database());
        if (i >= 0) {
            triplets.emplace_back(static_cast<Eigen::Index>(j), i, Scalar(1));
        }
    }
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> K(
        static_cast<Eigen::Index>(other_kets.size()),
        static_cast<Eigen::Index>(this->get_number_of_kets()));
    K.setFromTriplets(triplets.begin(), triplets.end());
    K.makeCompressed();

    return other->get_coefficients().adjoint() * K * this->get_coefficients();
}

template <typename Scalar>
Eigen::VectorX<typename BasisAtom<Scalar>::real_t>
BasisAtom<Scalar>::get_overlaps(std::shared_ptr<const ket_t> ket) const {
    return get_amplitudes(ket).cwiseAbs2();
}

template <typename Scalar>
Eigen::SparseMatrix<typename BasisAtom<Scalar>::real_t, Eigen::RowMajor>
BasisAtom<Scalar>::get_overlaps(std::shared_ptr<const Type> other) const {
    return get_amplitudes(other).cwiseAbs2();
}

template <typename Scalar>
Eigen::VectorX<Scalar> BasisAtom<Scalar>::get_matrix_elements(std::shared_ptr<const ket_t> ket,
                                                              OperatorType type, int q) const {
    auto final = ket->template to_trivial_state<Scalar>();
    auto matrix_elements = this->get_database().get_matrix_elements_in_canonical_basis(
        this->shared_from_this(), final, type, q);
    matrix_elements =
        final->get_coefficients().adjoint() * matrix_elements * this->get_coefficients();

    assert(static_cast<size_t>(matrix_elements.rows()) == 1);
    assert(static_cast<size_t>(matrix_elements.cols()) == this->get_number_of_states());

    return matrix_elements.row(0);
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
