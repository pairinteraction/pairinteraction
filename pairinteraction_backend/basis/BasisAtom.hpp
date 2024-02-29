#pragma once

#include "basis/Basis.hpp"
#include "basis/BasisAtomCreator.hpp"
#include "ket/KetAtom.hpp"

#include <complex>
#include <limits>

template <typename T>
class BasisAtomCreator;

/**
 * @class BasisAtom
 *
 * @brief Class for creating a basis of atomic kets.
 *
 * @tparam T Complex number type.
 */
template <typename T>
class BasisAtom : public Basis<T> {
public:
    // const KetAtom<real_t<T>> &Iterator::operator*() const override; // TODO can we make such a
    // thing work?
private:
    friend class BasisAtomCreator<T>;
    BasisAtom(std::vector<std::shared_ptr<const Ket<real_t<T>>>> &&kets);
};

extern template class BasisAtom<float>;
extern template class BasisAtom<double>;
extern template class BasisAtom<std::complex<float>>;
extern template class BasisAtom<std::complex<double>>;
