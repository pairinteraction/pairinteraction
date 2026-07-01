// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisPair.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/tensor.hpp"

#include <algorithm>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <vector>

namespace pairinteraction {
namespace {
template <typename Scalar>
bool have_equivalent_atomic_basis(const std::shared_ptr<const BasisAtom<Scalar>> &basis_a,
                                  const std::shared_ptr<const BasisAtom<Scalar>> &basis_b) {
    if (basis_a == basis_b) {
        return true;
    }
    if (&basis_a->get_database() != &basis_b->get_database()) {
        return false;
    }
    if (basis_a->get_number_of_kets() != basis_b->get_number_of_kets() ||
        basis_a->get_number_of_states() != basis_b->get_number_of_states()) {
        return false;
    }
    if (basis_a->get_canonical_basis_id() != basis_b->get_canonical_basis_id()) {
        return false;
    }

    const auto &basis_a_kets = basis_a->get_kets();
    const auto &basis_b_kets = basis_b->get_kets();
    for (size_t i = 0; i < basis_a_kets.size(); ++i) {
        if (basis_a_kets[i]->get_id_in_database() != basis_b_kets[i]->get_id_in_database()) {
            return false;
        }
    }

    return basis_a->get_coefficients().isApprox(basis_b->get_coefficients());
}
} // namespace

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
std::shared_ptr<const typename BasisPair<Scalar>::Type>
BasisPair<Scalar>::merge(std::shared_ptr<const Type> other) const {
    if (!have_equivalent_atomic_basis(basis1, other->basis1) ||
        !have_equivalent_atomic_basis(basis2, other->basis2)) {
        throw std::invalid_argument(
            "Cannot merge pair bases with different ordered atomic states.");
    }

    ketvec_t merged_kets;
    merged_kets.reserve(this->kets.size() + other->kets.size());
    map_indices_t merged_state_indices;
    merged_state_indices.reserve(state_indices_to_ket_index.size() +
                                 other->state_indices_to_ket_index.size());

    auto append_unique_kets = [&](const Type &basis) {
        std::vector<std::vector<size_t>> state_indices_by_ket(basis.kets.size());
        for (const auto &[state_indices, ket_index] : basis.state_indices_to_ket_index) {
            state_indices_by_ket.at(ket_index) = state_indices;
        }
        for (size_t ket_index = 0; ket_index < basis.kets.size(); ++ket_index) {
            const auto &state_indices = state_indices_by_ket[ket_index];
            if (state_indices.size() != 2) {
                throw std::runtime_error("Invalid atomic state indices in pair basis.");
            }
            if (merged_state_indices.contains(state_indices)) {
                continue;
            }
            merged_state_indices.try_emplace(state_indices, merged_kets.size());
            merged_kets.push_back(basis.kets[ket_index]);
        }
    };
    append_unique_kets(*this);
    append_unique_kets(*other);

    const size_t number_of_states1 = basis1->get_number_of_states();
    const size_t number_of_states2 = basis2->get_number_of_states();
    std::vector<size_t> minimum_indices2(number_of_states1, number_of_states2);
    std::vector<size_t> maximum_indices2(number_of_states1, 0);
    std::vector<bool> has_index(number_of_states1, false);
    for (const auto &[state_indices, ket_index] : merged_state_indices) {
        static_cast<void>(ket_index);
        const size_t state_index1 = state_indices[0];
        const size_t state_index2 = state_indices[1];
        minimum_indices2[state_index1] = std::min(minimum_indices2[state_index1], state_index2);
        maximum_indices2[state_index1] = std::max(maximum_indices2[state_index1], state_index2);
        has_index[state_index1] = true;
    }

    map_range_t merged_ranges;
    merged_ranges.reserve(number_of_states1);
    for (size_t state_index1 = 0; state_index1 < number_of_states1; ++state_index1) {
        if (has_index[state_index1]) {
            merged_ranges.try_emplace(
                state_index1,
                range_t(minimum_indices2[state_index1], maximum_indices2[state_index1] + 1));
        } else {
            merged_ranges.try_emplace(state_index1, range_t(0, 0));
        }
    }

    return std::make_shared<const Type>(Private(), std::move(merged_kets), std::move(merged_ranges),
                                        std::move(merged_state_indices), basis1, basis2);
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
