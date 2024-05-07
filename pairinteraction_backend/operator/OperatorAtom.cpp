#include "operator/OperatorAtom.hpp"
#include "basis/BasisAtom.hpp"
#include "database/Database.hpp"

template <typename Scalar>
OperatorAtom<Scalar>::OperatorAtom(std::shared_ptr<const basis_t> basis, OperatorType type, int q)
    : Operator<OperatorAtom<Scalar>>(basis), type(type), q(q) {
    *this = basis->get_database().get_operator(basis, type, q);
}

template <typename Scalar>
OperatorAtom<Scalar>::OperatorAtom(std::shared_ptr<const basis_t> basis, OperatorType type, int q,
                                   Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &&matrix)
    : Operator<OperatorAtom<Scalar>>(basis), type(type), q(q) {
    this->matrix = std::move(matrix);
}

// Explicit instantiations
template class OperatorAtom<float>;
template class OperatorAtom<double>;
template class OperatorAtom<std::complex<float>>;
template class OperatorAtom<std::complex<double>>;
