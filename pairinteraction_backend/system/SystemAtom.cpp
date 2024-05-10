#include "system/SystemAtom.hpp"
#include "operator/OperatorAtom.hpp"

template <typename Scalar>
SystemAtom<Scalar>::SystemAtom(std::shared_ptr<const basis_t> basis)
    : System<SystemAtom<Scalar>>(basis) {}

// Explicit instantiations
template class SystemAtom<float>;
template class SystemAtom<double>;
template class SystemAtom<std::complex<float>>;
template class SystemAtom<std::complex<double>>;
