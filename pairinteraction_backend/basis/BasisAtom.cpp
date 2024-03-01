#include "basis/BasisAtom.hpp"

template <typename Scalar>
BasisAtom<Scalar>::BasisAtom(KetPtrVec &&kets)
    : kets(std::move(kets)), Basis<BasisAtom<Scalar>>() {}

// Explicit instantiations
template class BasisAtom<float>;
template class BasisAtom<double>;
template class BasisAtom<std::complex<float>>;
template class BasisAtom<std::complex<double>>;
