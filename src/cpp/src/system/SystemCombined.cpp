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
#include <oneapi/tbb.h>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
SystemCombined<Scalar>::SystemCombined(std::shared_ptr<const basis_t> basis)
    : System<SystemCombined<Scalar>>(std::move(basis)) {}

template <typename Scalar>
SystemCombined<Scalar> &SystemCombined<Scalar>::set_distance(real_t distance) {
    this->distance = distance;
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
    if (distance != 0) {

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
            (-2 * matrix_zero_zero - matrix_plus_minus - matrix_minus_plus) / std::pow(distance, 3);
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
    real_t numerical_precision = 10 * std::numeric_limits<real_t>::epsilon();

    size_t number_of_states2 = basis->get_basis2()->get_number_of_states();

    oneapi::tbb::concurrent_vector<Eigen::Triplet<Scalar>> triplets;

    // Loop over the rows of the first matrix in parallel (outer index == row)
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<Eigen::Index>(0, matrix1.outerSize()), [&](const auto &range) {
            for (Eigen::Index row1 = range.begin(); row1 != range.end(); ++row1) {

                const auto &range_row2 = basis->get_index_range(row1);

                // Loop over the rows of the second matrix that are energetically allowed
                for (auto row2 = static_cast<Eigen::Index>(range_row2.min());
                     row2 < static_cast<Eigen::Index>(range_row2.max()); ++row2) {

                    size_t row_ket_id = row1 * number_of_states2 + row2;
                    if (!basis->has_ket_index(row_ket_id)) {
                        continue;
                    }
                    Eigen::Index row = basis->get_ket_index(row_ket_id);

                    // Loop over the non-zero column elements of the first matrix
                    for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it1(
                             matrix1, row1);
                         it1; ++it1) {

                        Eigen::Index col1 = it1.col();
                        Scalar value1 = it1.value();

                        // Calculate the minimum and maximum values of the index pointer of the
                        // second matrix
                        Eigen::Index begin_idxptr2 = matrix2.outerIndexPtr()[row2];
                        Eigen::Index end_idxptr2 = matrix2.outerIndexPtr()[row2 + 1];

                        // The minimum value is chosen such that we start with an energetically
                        // allowed column
                        const auto &range_col2 = basis->get_index_range(it1.index());
                        begin_idxptr2 +=
                            std::distance(matrix2.innerIndexPtr() + begin_idxptr2,
                                          std::lower_bound(matrix2.innerIndexPtr() + begin_idxptr2,
                                                           matrix2.innerIndexPtr() + end_idxptr2,
                                                           range_col2.min()));

                        // Loop over the non-zero column elements of the second matrix that are
                        // energetically allowed (we break the loop if the index pointer corresponds
                        // to a column that is not energetically allowed)
                        for (Eigen::Index idxptr2 = begin_idxptr2; idxptr2 < end_idxptr2;
                             ++idxptr2) {

                            Eigen::Index col2 = matrix2.innerIndexPtr()[idxptr2];
                            size_t col_ket_id = col1 * number_of_states2 + col2;
                            if (!basis->has_ket_index(col_ket_id)) {
                                continue;
                            }
                            if (col2 >= static_cast<Eigen::Index>(range_col2.max())) {
                                break;
                            }
                            Scalar value2 = matrix2.valuePtr()[idxptr2];
                            Eigen::Index col = basis->get_ket_index(col_ket_id);

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
