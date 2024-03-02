#include "basis/BasisAtom.hpp"

template <typename Scalar>
BasisAtom<Scalar>::BasisAtom(ketvec_t &&kets) : Basis<BasisAtom<Scalar>>(std::move(kets)) {}

// Explicit instantiations
template class BasisAtom<float>;
template class BasisAtom<double>;
template class BasisAtom<std::complex<float>>;
template class BasisAtom<std::complex<double>>;
