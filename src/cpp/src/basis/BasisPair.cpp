// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisPair.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisPairCreator.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/tensor.hpp"

#include <cassert>
#include <memory>
#include <oneapi/tbb.h>
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
Eigen::VectorX<Scalar>
BasisPair<Scalar>::get_amplitudes(std::shared_ptr<const KetAtom> ket1,
                                  std::shared_ptr<const KetAtom> ket2) const {
    return get_amplitudes(ket1->template to_trivial_state<Scalar>(),
                          ket2->template to_trivial_state<Scalar>())
        .transpose();
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
BasisPair<Scalar>::get_amplitudes(std::shared_ptr<const Type> other) const {
    auto amplitudes_matrix1 = basis1->get_amplitudes(other->get_basis1());
    auto amplitudes_matrix2 = basis2->get_amplitudes(other->get_basis2());

    auto result = utils::calculate_tensor_product_in_canonical_basis(
        this->shared_from_this(), other, amplitudes_matrix1, amplitudes_matrix2);
    assert(static_cast<size_t>(result.rows()) == other->get_number_of_kets());
    assert(static_cast<size_t>(result.cols()) == this->get_number_of_kets());

    return other->get_coefficients().adjoint() * result * this->get_coefficients();
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
BasisPair<Scalar>::get_amplitudes(std::shared_ptr<const BasisAtom<Scalar>> other1,
                                  std::shared_ptr<const BasisAtom<Scalar>> other2) const {
    auto system1 = SystemAtom<Scalar>(other1);
    auto system2 = SystemAtom<Scalar>(other2);
    auto other = BasisPairCreator<Scalar>().add(system1).add(system2).create();
    assert(other->get_number_of_states() ==
           other1->get_number_of_states() * other2->get_number_of_states());
    assert(other->get_number_of_kets() ==
           other1->get_number_of_states() * other2->get_number_of_states());
    return get_amplitudes(other);
}

template <typename Scalar>
Eigen::VectorX<typename BasisPair<Scalar>::real_t>
BasisPair<Scalar>::get_overlaps(std::shared_ptr<const KetAtom> ket1,
                                std::shared_ptr<const KetAtom> ket2) const {
    return get_amplitudes(ket1, ket2).cwiseAbs2();
}

template <typename Scalar>
Eigen::SparseMatrix<typename BasisPair<Scalar>::real_t, Eigen::RowMajor>
BasisPair<Scalar>::get_overlaps(std::shared_ptr<const Type> other) const {
    return get_amplitudes(other).cwiseAbs2();
}

template <typename Scalar>
Eigen::SparseMatrix<typename BasisPair<Scalar>::real_t, Eigen::RowMajor>
BasisPair<Scalar>::get_overlaps(std::shared_ptr<const BasisAtom<Scalar>> other1,
                                std::shared_ptr<const BasisAtom<Scalar>> other2) const {
    return get_amplitudes(other1, other2).cwiseAbs2();
}

template <typename Scalar>
Eigen::VectorX<Scalar> BasisPair<Scalar>::get_matrix_elements(std::shared_ptr<const ket_t> /*ket*/,
                                                              OperatorType /*type*/,
                                                              int /*q*/) const {
    throw std::invalid_argument("It is required to specify one operator for each atom.");
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
BasisPair<Scalar>::get_matrix_elements(std::shared_ptr<const Type> /*final*/, OperatorType /*type*/,
                                       int /*q*/) const {
    throw std::invalid_argument("It is required to specify one operator for each atom.");
}

template <typename Scalar>
Eigen::VectorX<Scalar>
BasisPair<Scalar>::get_matrix_elements(std::shared_ptr<const ket_t> ket, OperatorType type1,
                                       OperatorType type2, int q1, int q2) const {
    // Construct a pair basis containing only the pair ket
    auto final = this->get_canonical_state_from_ket(ket);
    assert(final->get_number_of_states() == 1);

    return this->get_matrix_elements(final, type1, type2, q1, q2).row(0);
}

template <typename Scalar>
Eigen::VectorX<Scalar>
BasisPair<Scalar>::get_matrix_elements(std::shared_ptr<const KetAtom> ket1,
                                       std::shared_ptr<const KetAtom> ket2, OperatorType type1,
                                       OperatorType type2, int q1, int q2) const {
    // Construct a pair basis with the two single-atom kets
    auto final1 = this->get_basis1()->get_canonical_state_from_ket(ket1);
    auto final2 = this->get_basis2()->get_canonical_state_from_ket(ket2);
    auto system1 = SystemAtom<Scalar>(final1);
    auto system2 = SystemAtom<Scalar>(final2);
    auto final = BasisPairCreator<Scalar>().add(system1).add(system2).create();
    assert(final->get_number_of_states() == 1);
    assert(final->get_number_of_kets() == 1);

    return this->get_matrix_elements(final, type1, type2, q1, q2).row(0);
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
BasisPair<Scalar>::get_matrix_elements(std::shared_ptr<const Type> final, OperatorType type1,
                                       OperatorType type2, int q1, int q2) const {
    auto initial1 = this->get_basis1();
    auto initial2 = this->get_basis2();
    auto final1 = final->get_basis1();
    auto final2 = final->get_basis2();

    auto matrix_elements1 = initial1->get_database().get_matrix_elements_in_canonical_basis(
        initial1, final1, type1, q1);
    matrix_elements1 =
        final1->get_coefficients().adjoint() * matrix_elements1 * initial1->get_coefficients();
    auto matrix_elements2 = initial2->get_database().get_matrix_elements_in_canonical_basis(
        initial2, final2, type2, q2);
    matrix_elements2 =
        final2->get_coefficients().adjoint() * matrix_elements2 * initial2->get_coefficients();
    auto matrix_elements = utils::calculate_tensor_product_in_canonical_basis(
        this->shared_from_this(), final, matrix_elements1, matrix_elements2);
    assert(static_cast<size_t>(matrix_elements.rows()) == final->get_number_of_kets());
    assert(static_cast<size_t>(matrix_elements.cols()) == this->get_number_of_kets());

    return final->get_coefficients().adjoint() * matrix_elements * this->get_coefficients();
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
BasisPair<Scalar>::get_matrix_elements(std::shared_ptr<const BasisAtom<Scalar>> final1,
                                       std::shared_ptr<const BasisAtom<Scalar>> final2,
                                       OperatorType type1, OperatorType type2, int q1,
                                       int q2) const {
    // Construct a pair basis with the two single-atom bases
    auto system1 = SystemAtom<Scalar>(final1);
    auto system2 = SystemAtom<Scalar>(final2);
    auto final = BasisPairCreator<Scalar>().add(system1).add(system2).create();
    assert(final->get_number_of_states() ==
           final1->get_number_of_states() * final2->get_number_of_states());
    assert(final->get_number_of_kets() ==
           final1->get_number_of_states() * final2->get_number_of_states());

    return this->get_matrix_elements(final, type1, type2, q1, q2);
}

// Explicit instantiations
template class BasisPair<double>;
template class BasisPair<std::complex<double>>;
} // namespace pairinteraction
