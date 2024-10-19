#include "pairinteraction/basis/BasisClassicalLight.hpp"

namespace pairinteraction {
template <typename Scalar>
BasisClassicalLight<Scalar>::BasisClassicalLight(Private /*unused*/, ketvec_t &&kets,
                                                 std::string &&id_of_kets)
    : Basis<BasisClassicalLight<Scalar>>(std::move(kets), std::move(id_of_kets)) {}

// Explicit instantiations
template class BasisClassicalLight<float>;
template class BasisClassicalLight<double>;
template class BasisClassicalLight<std::complex<float>>;
template class BasisClassicalLight<std::complex<double>>;
} // namespace pairinteraction
