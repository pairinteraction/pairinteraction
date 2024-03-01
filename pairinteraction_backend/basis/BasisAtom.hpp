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
 * @tparam Scalar Complex number type.
 */
template <typename Scalar>
class BasisAtom : public Basis<BasisAtom<Scalar>> {
public:
    using Real = real_t<Scalar>;
    using MyScalar = Scalar;
    using MyKet = KetAtom<real_t<Scalar>>;
    using KetPtrVec = std::vector<std::shared_ptr<const KetAtom<Real>>>;

    const KetPtrVec &get_kets() const { return kets; } // TODO to source or remove

    const KetAtom<Real> &get_ket(size_t index) const {
        return *kets[index];
    } // TODO to source or remove

private:
    friend class BasisAtomCreator<Scalar>;
    BasisAtom(KetPtrVec &&kets);
    KetPtrVec kets; // TODO move to Basis.hpp
};

extern template class BasisAtom<float>;
extern template class BasisAtom<double>;
extern template class BasisAtom<std::complex<float>>;
extern template class BasisAtom<std::complex<double>>;
