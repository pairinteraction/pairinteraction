#include "pairinteraction/basis/BasisCombined.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetCombined.hpp"

#include <memory>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
BasisCombined<Scalar>::BasisCombined(ketvec_t &&kets)
    : Basis<BasisCombined<Scalar>>(std::move(kets)) {}

// Explicit instantiations
template class BasisCombined<float>;
template class BasisCombined<double>;
template class BasisCombined<std::complex<float>>;
template class BasisCombined<std::complex<double>>;
} // namespace pairinteraction
