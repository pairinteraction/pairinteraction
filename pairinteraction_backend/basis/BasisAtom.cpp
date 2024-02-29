#include "basis/BasisAtom.hpp"

template <typename T>
BasisAtom<T>::BasisAtom(std::vector<std::shared_ptr<const Ket<real_t<T>>>> &&kets)
    : Basis<T>(std::move(kets)) {}

template <typename T>
BasisAtom<T>::IteratorAtom::IteratorAtom(const BasisAtom<T> &basis, size_t index)
    : Basis<T>::Iterator(basis, index) {}

template <typename T>
const KetAtom<real_t<T>> &BasisAtom<T>::IteratorAtom::operator*() const {
    return static_cast<const KetAtom<real_t<T>> &>(Basis<T>::Iterator::operator*());
}

template <typename T>
typename BasisAtom<T>::IteratorAtom BasisAtom<T>::begin() const {
    return IteratorAtom(*this, 0);
}

template <typename T>
typename BasisAtom<T>::IteratorAtom BasisAtom<T>::end() const {
    return IteratorAtom(*this, this->kets.size());
}

// Explicit instantiations
template class BasisAtom<float>;
template class BasisAtom<double>;
template class BasisAtom<std::complex<float>>;
template class BasisAtom<std::complex<double>>;
