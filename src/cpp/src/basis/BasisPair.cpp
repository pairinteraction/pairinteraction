// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
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
    if (state_indices_to_ket_index.count({state_index1, state_index2}) == 0) {
        return -1;
    }
    return state_indices_to_ket_index.at({state_index1, state_index2});
}

template <typename Scalar>
Eigen::VectorX<Scalar>
BasisPair<Scalar>::get_amplitudes(std::shared_ptr<const KetAtom> ket1,
                                  std::shared_ptr<const KetAtom> ket2) const {
    return get_amplitudes(basis1->get_canonical_state_from_ket(ket1),
                          basis2->get_canonical_state_from_ket(ket2))
        .transpose();
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
BasisPair<Scalar>::get_amplitudes(std::shared_ptr<const BasisAtom<Scalar>> other1,
                                  std::shared_ptr<const BasisAtom<Scalar>> other2) const {
    if (other1->get_id_of_kets() != basis1->get_id_of_kets() ||
        other2->get_id_of_kets() != basis2->get_id_of_kets()) {
        throw std::invalid_argument("The other objects must be expressed using the same kets.");
    }

    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();

    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> coefficients1 =
        basis1->get_coefficients().adjoint() * other1->get_coefficients();
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> coefficients2 =
        basis2->get_coefficients().adjoint() * other2->get_coefficients();

    oneapi::tbb::concurrent_vector<Eigen::Triplet<Scalar>> triplets;

    // Loop over the rows of the first coefficient matrix
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<Eigen::Index>(0, coefficients1.outerSize()),
        [&](const auto &range) {
            for (Eigen::Index row1 = range.begin(); row1 != range.end(); ++row1) {

                const auto &range_row2 = this->get_index_range(row1);

                // Loop over the rows of the second coefficient matrix that are energetically
                // allowed
                for (auto row2 = static_cast<Eigen::Index>(range_row2.min());
                     row2 < static_cast<Eigen::Index>(range_row2.max()); ++row2) {

                    Eigen::Index row = get_ket_index_from_tuple(row1, row2);
                    if (row < 0) {
                        continue;
                    }

                    // Loop over the non-zero column elements of the first coefficient matrix
                    for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it1(
                             coefficients1, row1);
                         it1; ++it1) {

                        Eigen::Index col1 = it1.col();
                        Scalar value1 = it1.value();

                        // Loop over the non-zero column elements of the second coefficient matrix
                        for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator
                                 it2(coefficients2, row2);
                             it2; ++it2) {

                            Eigen::Index col2 = it2.col();
                            Scalar value2 = it2.value();
                            Eigen::Index col = col1 * coefficients2.cols() + col2;

                            // Store the entry
                            Scalar value = value1 * value2;
                            if (std::abs(value) > numerical_precision) {
                                triplets.emplace_back(row, col, value);
                            }
                        }
                    }
                }
            }
        });

    // Construct the combined matrix from the triplets
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> matrix(this->get_number_of_kets(),
                                                        other1->get_number_of_states() *
                                                            other2->get_number_of_states());
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix.makeCompressed();

    return matrix.adjoint() * this->get_coefficients();
}

template <typename Scalar>
Eigen::VectorX<typename BasisPair<Scalar>::real_t>
BasisPair<Scalar>::get_overlaps(std::shared_ptr<const KetAtom> ket1,
                                std::shared_ptr<const KetAtom> ket2) const {
    return get_amplitudes(ket1, ket2).cwiseAbs2();
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
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
BasisPair<Scalar>::get_matrix_elements(std::shared_ptr<const Type> final, OperatorType type1,
                                       OperatorType type2, int q1, int q2) const {
    auto initial1 = this->get_basis1();
    auto initial2 = this->get_basis2();
    auto final1 = final->get_basis1();
    auto final2 = final->get_basis2();

    auto matrix_elements1 =
        initial1->get_database().get_matrix_elements(initial1, final1, type1, q1);
    auto matrix_elements2 =
        initial2->get_database().get_matrix_elements(initial2, final2, type2, q2);
    auto matrix_elements = utils::calculate_tensor_product(this->shared_from_this(), final,
                                                           matrix_elements1, matrix_elements2);
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

// Explicit instantiations
template class BasisPair<double>;
template class BasisPair<std::complex<double>>;
} // namespace pairinteraction
