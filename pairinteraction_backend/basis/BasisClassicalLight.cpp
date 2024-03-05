#include "basis/BasisClassicalLight.hpp"

template <typename Scalar>
BasisClassicalLight<Scalar>::BasisClassicalLight(ketvec_t &&kets)
    : Basis<BasisClassicalLight<Scalar>>(std::move(kets)) {}

// Explicit instantiations
template class BasisClassicalLight<float>;
template class BasisClassicalLight<double>;
template class BasisClassicalLight<std::complex<float>>;
template class BasisClassicalLight<std::complex<double>>;
