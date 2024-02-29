#include "basis/BasisAtom.hpp"

template <typename T>
BasisAtom<T>::BasisAtom(std::vector<std::shared_ptr<const Ket<real_t<T>>>> &&kets)
    : Basis<T>(std::move(kets)) {}

// Explicit instantiations
template class BasisAtom<float>;
template class BasisAtom<double>;
template class BasisAtom<std::complex<float>>;
template class BasisAtom<std::complex<double>>;
