#include "operator/OperatorAtom.hpp"
#include "basis/BasisAtom.hpp"

template <typename Scalar>
OperatorAtom<Scalar>::OperatorAtom(const basis_t &basis, OperatorType type, int q)
    : Operator<OperatorAtom<Scalar>>(basis), type(type), q(q) {}

// Explicit instantiations
template class OperatorAtom<float>;
template class OperatorAtom<double>;
template class OperatorAtom<std::complex<float>>;
template class OperatorAtom<std::complex<double>>;
