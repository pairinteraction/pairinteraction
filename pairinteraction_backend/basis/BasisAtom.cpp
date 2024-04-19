#include "basis/BasisAtom.hpp"

template <typename Scalar>
BasisAtom<Scalar>::BasisAtom(ketvec_t &&kets, Database &database)
    : Basis<BasisAtom<Scalar>>(std::move(kets)), database(database) {}

// Explicit instantiations
template class BasisAtom<float>;
template class BasisAtom<double>;
template class BasisAtom<std::complex<float>>;
template class BasisAtom<std::complex<double>>;
