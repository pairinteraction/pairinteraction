#include "pairinteraction/system/SystemCombined.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisCombined.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetCombined.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"
#include "pairinteraction/operator/OperatorCombined.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/SparseCore>
#include <algorithm>
#include <complex>
#include <limits>
#include <memory>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
SystemCombined<Scalar>::SystemCombined(std::shared_ptr<const basis_t> basis)
    : System<SystemCombined<Scalar>>(basis) {}

template <typename Scalar>
SystemCombined<Scalar> &SystemCombined<Scalar>::set_interatomic_distance(real_t distance) {
    this->interatomic_distance = distance;
    this->hamiltonian_requires_construction = true;
    return *this;
}

template <typename Scalar>
void SystemCombined<Scalar>::construct_hamiltonian() const {
    auto basis = this->hamiltonian->get_basis();
    auto basis1 = basis->get_basis1();
    auto basis2 = basis->get_basis2();

    // Construct the unperturbed Hamiltonian
    this->hamiltonian = std::make_unique<OperatorCombined<Scalar>>(basis, OperatorType::ENERGY);
    this->hamiltonian_is_diagonal = true;
    bool sort_by_quantum_number_f = basis->has_quantum_number_f();
    bool sort_by_quantum_number_m = basis->has_quantum_number_m();
    bool sort_by_parity = basis->has_parity();

    // Dipole-dipole interaction along the z-axis
    if (interatomic_distance != 0) {

        auto d1_plus = OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_DIPOLE, 1);
        auto d1_minus = OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_DIPOLE, -1);
        auto d1_zero = OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_DIPOLE, 0);
        auto d2_plus = OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_DIPOLE, 1);
        auto d2_minus = OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_DIPOLE, -1);
        auto d2_zero = OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_DIPOLE, 0);

        auto matrix_zero_zero =
            calculate_tensor_product(basis, d1_zero.get_matrix(), d2_zero.get_matrix());
        auto matrix_plus_minus =
            calculate_tensor_product(basis, d1_plus.get_matrix(), d2_minus.get_matrix());
        auto matrix_minus_plus =
            calculate_tensor_product(basis, d1_minus.get_matrix(), d2_plus.get_matrix());

        this->hamiltonian->get_matrix() +=
            (-2 * matrix_zero_zero - matrix_plus_minus - matrix_minus_plus) /
            std::pow(interatomic_distance, 3);
        this->hamiltonian_is_diagonal = false;
        sort_by_quantum_number_f = false;
        sort_by_parity = false;
    }

    // Store which labels can be used to block-diagonalize the Hamiltonian
    this->blockdiagonalizing_labels.clear();
    if (sort_by_quantum_number_f) {
        this->blockdiagonalizing_labels.push_back(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
    }
    if (sort_by_quantum_number_m) {
        this->blockdiagonalizing_labels.push_back(TransformationType::SORT_BY_QUANTUM_NUMBER_M);
    }
    if (sort_by_parity) {
        this->blockdiagonalizing_labels.push_back(TransformationType::SORT_BY_PARITY);
    }
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor> SystemCombined<Scalar>::calculate_tensor_product(
    const std::shared_ptr<const basis_t> &basis,
    const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix1,
    const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix2) {
    real_t precision = 10 * std::numeric_limits<real_t>::epsilon();

    std::vector<Eigen::Triplet<Scalar>> triplets;

    // Loop over the rows of the first matrix
    for (Eigen::Index outer_idx1 = 0; outer_idx1 < matrix1.outerSize(); ++outer_idx1) {

        // Loop over the non-zero column elements of the first matrix
        for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it1(matrix1,
                                                                                      outer_idx1);
             it1; ++it1) {

            Eigen::Index row1 = it1.row();
            Eigen::Index col1 = it1.col();
            Scalar value1 = it1.value();

            const auto &range_outer_idx2 = basis->get_index_range(outer_idx1);
            const auto &range_inner_idx2 = basis->get_index_range(it1.index());

            // Loop over the rows of the second matrix that are energetically allowed
            for (Eigen::Index outer_idx2 = static_cast<Eigen::Index>(range_outer_idx2.min());
                 outer_idx2 < static_cast<Eigen::Index>(range_outer_idx2.max()); ++outer_idx2) {

                if (!basis->are_valid_indices(outer_idx1, outer_idx2)) {
                    continue;
                }

                // Calculate the minimum and maximum values of the index pointer of the second
                // matrix
                Eigen::Index begin_idxptr2 = matrix2.outerIndexPtr()[outer_idx2];
                Eigen::Index end_idxptr2 = matrix2.outerIndexPtr()[outer_idx2 + 1];

                // The minimum value is chosen such that we start with an energetically allowed
                // column
                begin_idxptr2 +=
                    std::distance(matrix2.innerIndexPtr() + begin_idxptr2,
                                  std::lower_bound(matrix2.innerIndexPtr() + begin_idxptr2,
                                                   matrix2.innerIndexPtr() + end_idxptr2,
                                                   range_inner_idx2.min()));

                // Loop over the non-zero column elements of the second matrix that are
                // energetically allowed (we break the loop if the index pointer corresponds to
                // a column that is not energetically allowed)
                for (Eigen::Index idxptr2 = begin_idxptr2; idxptr2 < end_idxptr2; ++idxptr2) {

                    Eigen::Index col2 = matrix2.innerIndexPtr()[idxptr2];
                    if (!basis->are_valid_indices(col1, col2)) {
                        continue;
                    }
                    if (col2 >= static_cast<Eigen::Index>(range_inner_idx2.max())) {
                        break;
                    }
                    Eigen::Index row2 = outer_idx2;
                    Scalar value2 = matrix2.valuePtr()[idxptr2];

                    // Calculate the row and column index of the entry in the combined matrix
                    Eigen::Index row = basis->get_combined_index(row1, row2);
                    Eigen::Index col = basis->get_combined_index(col1, col2);

                    // Calculate the value of the entry in the combined matrix
                    Scalar value = value1 * value2;

                    // Store the entry
                    if (std::abs(value) > precision) {
                        triplets.emplace_back(row, col, value);
                    }
                }
            }
        }
    }

    // Construct the combined matrix from the triplets
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> matrix(basis->get_number_of_states(),
                                                        basis->get_number_of_states());
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix.makeCompressed();

    return matrix;
}

// Explicit instantiations
template class SystemCombined<float>;
template class SystemCombined<double>;
template class SystemCombined<std::complex<float>>;
template class SystemCombined<std::complex<double>>;
} // namespace pairinteraction
