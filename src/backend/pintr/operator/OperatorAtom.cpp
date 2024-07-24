#include "pintr/operator/OperatorAtom.hpp"

#include "pintr/basis/BasisAtom.hpp"
#include "pintr/database/Database.hpp"
#include "pintr/enums/OperatorType.hpp"

#include <limits>

template <typename Scalar>
OperatorAtom<Scalar>::OperatorAtom(std::shared_ptr<const basis_t> basis)
    : Operator<OperatorAtom<Scalar>>(basis), type(OperatorType::ZERO),
      q(std::numeric_limits<int>::max()) {
    this->matrix = Eigen::SparseMatrix<Scalar, Eigen::RowMajor>(basis->get_number_of_states(),
                                                                basis->get_number_of_states());
}

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
