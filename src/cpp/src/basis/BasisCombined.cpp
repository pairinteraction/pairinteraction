#include "pairinteraction/basis/BasisCombined.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetCombined.hpp"
#include "pairinteraction/utils/Range.hpp"

#include <memory>
#include <oneapi/tbb.h>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
BasisCombined<Scalar>::BasisCombined(Private /*unused*/, ketvec_t &&kets, std::string &&id_of_kets,
                                     map_range_t &&map_range_of_state_index2,
                                     std::shared_ptr<const BasisAtom<Scalar>> basis1,
                                     std::shared_ptr<const BasisAtom<Scalar>> basis2)
    : Basis<BasisCombined<Scalar>>(std::move(kets), std::move(id_of_kets)),
      map_range_of_state_index2(std::move(map_range_of_state_index2)), basis1(std::move(basis1)),
      basis2(std::move(basis2)) {}

template <typename Scalar>
const typename BasisCombined<Scalar>::range_t &
BasisCombined<Scalar>::get_index_range(size_t state_index1) const {
    return map_range_of_state_index2.at(state_index1);
}

template <typename Scalar>
std::shared_ptr<const BasisAtom<Scalar>> BasisCombined<Scalar>::get_basis1() const {
    return basis1;
}

template <typename Scalar>
std::shared_ptr<const BasisAtom<Scalar>> BasisCombined<Scalar>::get_basis2() const {
    return basis2;
}

template <typename Scalar>
int BasisCombined<Scalar>::get_ket_index_from_tuple(size_t state_index1,
                                                    size_t state_index2) const {
    size_t ket_id = state_index1 * basis2->get_number_of_states() + state_index2;
    // TODO this method should not rely on get_ket_index_from_id
    return this->get_ket_index_from_id(ket_id);
}

template <typename Scalar>
Eigen::VectorX<Scalar>
BasisCombined<Scalar>::get_amplitudes(std::shared_ptr<const KetAtom<real_t>> ket1,
                                      std::shared_ptr<const KetAtom<real_t>> ket2) const {
    return get_amplitudes(basis1->get_canonical_state_from_ket(ket1),
                          basis2->get_canonical_state_from_ket(ket2))
        .transpose();
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
BasisCombined<Scalar>::get_amplitudes(std::shared_ptr<const BasisAtom<Scalar>> other1,
                                      std::shared_ptr<const BasisAtom<Scalar>> other2) const {
    if (other1->get_id_of_kets() != basis1->get_id_of_kets() ||
        other2->get_id_of_kets() != basis2->get_id_of_kets()) {
        throw std::invalid_argument("The other objects must be expressed using the same kets.");
    }

    real_t numerical_precision = 10 * std::numeric_limits<real_t>::epsilon();

    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> coefficients1 =
        basis1->get_coefficients().adjoint() * other1->get_coefficients();
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> coefficients2 =
        basis2->get_coefficients().adjoint() * other2->get_coefficients();

    std::vector<Eigen::Triplet<Scalar>> triplets;

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
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> matrix(this->get_number_of_states(),
                                                        other1->get_number_of_states() *
                                                            other2->get_number_of_states());
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix.makeCompressed();

    return matrix.adjoint() * this->get_coefficients();
}

template <typename Scalar>
Eigen::VectorX<typename BasisCombined<Scalar>::real_t>
BasisCombined<Scalar>::get_overlaps(std::shared_ptr<const KetAtom<real_t>> ket1,
                                    std::shared_ptr<const KetAtom<real_t>> ket2) const {
    return get_amplitudes(ket1, ket2).cwiseAbs2();
}

template <typename Scalar>
Eigen::SparseMatrix<typename BasisCombined<Scalar>::real_t, Eigen::RowMajor>
BasisCombined<Scalar>::get_overlaps(std::shared_ptr<const BasisAtom<Scalar>> other1,
                                    std::shared_ptr<const BasisAtom<Scalar>> other2) const {
    return get_amplitudes(other1, other2).cwiseAbs2();
}

// Explicit instantiations
template class BasisCombined<float>;
template class BasisCombined<double>;
template class BasisCombined<std::complex<float>>;
template class BasisCombined<std::complex<double>>;
} // namespace pairinteraction
