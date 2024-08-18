#include "pintr/basis/BasisClassicalLight.hpp"

namespace pintr {
template <typename Scalar>
BasisClassicalLight<Scalar>::BasisClassicalLight(Private /*unused*/, ketvec_t &&kets)
    : Basis<BasisClassicalLight<Scalar>>(std::move(kets)) {}

// Explicit instantiations
template class BasisClassicalLight<float>;
template class BasisClassicalLight<double>;
template class BasisClassicalLight<std::complex<float>>;
template class BasisClassicalLight<std::complex<double>>;
} // namespace pintr
