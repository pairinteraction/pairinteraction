// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisPair.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/tensor.hpp"

#include <cassert>
#include <memory>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
BasisPair<Scalar>::BasisPair(Private /*unused*/, ketvec_t &&kets,
                             map_range_t &&map_range_of_state_index2,
                             map_indices_t &&state_indices_to_ket_index,
                             std::shared_ptr<const BasisAtom<Scalar>> basis1,
                             std::shared_ptr<const BasisAtom<Scalar>> basis2)
    : Basis<BasisPair<Scalar>>(std::move(kets)),
      map_range_of_state_index2(std::move(map_range_of_state_index2)),
      state_indices_to_ket_index(std::move(state_indices_to_ket_index)), basis1(std::move(basis1)),
      basis2(std::move(basis2)) {}

template <typename Scalar>
const typename BasisPair<Scalar>::range_t &
BasisPair<Scalar>::get_index_range(size_t state_index1) const {
    return map_range_of_state_index2.at(state_index1);
}

template <typename Scalar>
std::shared_ptr<const BasisAtom<Scalar>> BasisPair<Scalar>::get_basis1() const {
    return basis1;
}

template <typename Scalar>
std::shared_ptr<const BasisAtom<Scalar>> BasisPair<Scalar>::get_basis2() const {
    return basis2;
}

template <typename Scalar>
int BasisPair<Scalar>::get_ket_index_from_tuple(size_t state_index1, size_t state_index2) const {
    if (!state_indices_to_ket_index.contains({state_index1, state_index2})) {
        return -1;
    }
    return state_indices_to_ket_index.at({state_index1, state_index2});
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
BasisPair<Scalar>::get_matrix_elements(std::shared_ptr<const Type> final_state, OperatorType type1,
                                       OperatorType type2, int q1, int q2) const {
    auto initial1 = this->get_basis1();
    auto initial2 = this->get_basis2();
    auto final1 = final_state->get_basis1();
    auto final2 = final_state->get_basis2();

    auto matrix_elements1 = initial1->get_database().get_matrix_elements_in_canonical_basis(
        initial1, final1, type1, q1);
    matrix_elements1 =
        final1->get_coefficients().adjoint() * matrix_elements1 * initial1->get_coefficients();
    auto matrix_elements2 = initial2->get_database().get_matrix_elements_in_canonical_basis(
        initial2, final2, type2, q2);
    matrix_elements2 =
        final2->get_coefficients().adjoint() * matrix_elements2 * initial2->get_coefficients();
    auto matrix_elements = utils::calculate_tensor_product_in_canonical_basis(
        this->shared_from_this(), final_state, matrix_elements1, matrix_elements2);
    assert(static_cast<size_t>(matrix_elements.rows()) == final_state->get_number_of_kets());
    assert(static_cast<size_t>(matrix_elements.cols()) == this->get_number_of_kets());

    return final_state->get_coefficients().adjoint() * matrix_elements * this->get_coefficients();
}

// Explicit instantiations
template class BasisPair<double>;
template class BasisPair<std::complex<double>>;
} // namespace pairinteraction
