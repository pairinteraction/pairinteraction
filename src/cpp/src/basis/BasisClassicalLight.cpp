#include "pairinteraction/basis/BasisClassicalLight.hpp"

namespace pairinteraction {
template <typename Scalar>
BasisClassicalLight<Scalar>::BasisClassicalLight(Private /*unused*/, ketvec_t &&kets)
    : Basis<BasisClassicalLight<Scalar>>(std::move(kets)) {}

template <typename Scalar>
Eigen::VectorX<Scalar>
BasisClassicalLight<Scalar>::get_matrix_elements(std::shared_ptr<const ket_t> /*ket*/,
                                                 OperatorType /*type*/, int /*q*/) const {
    throw std::runtime_error("Not implemented.");
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
BasisClassicalLight<Scalar>::get_matrix_elements(std::shared_ptr<const Type> /*other*/,
                                                 OperatorType /*type*/, int /*q*/) const {
    throw std::runtime_error("Not implemented.");
}

// Explicit instantiations
template class BasisClassicalLight<double>;
template class BasisClassicalLight<std::complex<double>>;
} // namespace pairinteraction
