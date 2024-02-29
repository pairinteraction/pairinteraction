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
    class IteratorAtom : public Basis<T>::Iterator {
    public:
        IteratorAtom(const BasisAtom<T> &basis, size_t index);
        const KetAtom<real_t<T>> &operator*() const override;
    };
    IteratorAtom begin() const; // TODO this hides "Iterator begin() const", is it fine?
    IteratorAtom end() const;

private:
    friend class BasisAtomCreator<T>;
    BasisAtom(std::vector<std::shared_ptr<const Ket<real_t<T>>>> &&kets);
};

extern template class BasisAtom<float>;
extern template class BasisAtom<double>;
extern template class BasisAtom<std::complex<float>>;
extern template class BasisAtom<std::complex<double>>;
