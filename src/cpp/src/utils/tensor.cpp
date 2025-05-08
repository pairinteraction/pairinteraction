// SPDX-FileCopyrightText: 2025 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/utils/tensor.hpp"

#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/SparseCore>
#include <algorithm>
#include <limits>
#include <memory>
#include <oneapi/tbb.h>

namespace pairinteraction::utils {

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
calculate_tensor_product(const std::shared_ptr<const BasisPair<Scalar>> &basis_initial,
                         const std::shared_ptr<const BasisPair<Scalar>> &basis_final,
                         const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix1,
                         const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix2) {
    oneapi::tbb::concurrent_vector<Eigen::Triplet<Scalar>> triplets;

    // Loop over the rows of the first matrix in parallel (outer index == row)
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<Eigen::Index>(0, matrix1.outerSize()), [&](const auto &range) {
            for (Eigen::Index row1 = range.begin(); row1 != range.end(); ++row1) {

                const auto &range_row2 = basis_final->get_index_range(row1);

                // Loop over the rows of the second matrix that are energetically allowed
                for (auto row2 = static_cast<Eigen::Index>(range_row2.min());
                     row2 < static_cast<Eigen::Index>(range_row2.max()); ++row2) {

                    Eigen::Index row = basis_final->get_ket_index_from_tuple(row1, row2);
                    if (row < 0) {
                        continue;
                    }

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
                        const auto &range_col2 = basis_initial->get_index_range(it1.index());
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
                            if (col2 >= static_cast<Eigen::Index>(range_col2.max())) {
                                break;
                            }

                            Eigen::Index col = basis_initial->get_ket_index_from_tuple(col1, col2);
                            if (col < 0) {
                                continue;
                            }

                            Scalar value2 = matrix2.valuePtr()[idxptr2];

                            // Store the entry
                            Scalar value = value1 * value2;
                            triplets.emplace_back(row, col, value);
                        }
                    }
                }
            }
        });

    // Construct the combined matrix from the triplets
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> matrix(basis_final->get_number_of_kets(),
                                                        basis_initial->get_number_of_kets());
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix.makeCompressed();

    return matrix;
}

// Explicit instantiations
template Eigen::SparseMatrix<double, Eigen::RowMajor>
calculate_tensor_product(const std::shared_ptr<const BasisPair<double>> &,
                         const std::shared_ptr<const BasisPair<double>> &,
                         const Eigen::SparseMatrix<double, Eigen::RowMajor> &,
                         const Eigen::SparseMatrix<double, Eigen::RowMajor> &);
template Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>
calculate_tensor_product(const std::shared_ptr<const BasisPair<std::complex<double>>> &,
                         const std::shared_ptr<const BasisPair<std::complex<double>>> &,
                         const Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> &,
                         const Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> &);

} // namespace pairinteraction::utils
