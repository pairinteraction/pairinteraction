#pragma once

#include "basis/Basis.hpp"
#include "ket/KetAtom.hpp"
#include "utils/Traits.hpp"

#include <complex>

// Specialize BasisTraits for BasisAtom
template <typename T>
class BasisAtomCreator;

template <typename Scalar>
class BasisAtom;

template <typename Scalar>
struct internal::BasisTraits<BasisAtom<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename internal::NumTraits<Scalar>::real_t;
    using ket_t = KetAtom<real_t>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
};

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
    using Type = BasisAtom<Scalar>;
    using ketvec_t = typename internal::BasisTraits<Type>::ketvec_t;

private:
    friend class BasisAtomCreator<Scalar>;
    BasisAtom(ketvec_t &&kets);
};

extern template class BasisAtom<float>;
extern template class BasisAtom<double>;
extern template class BasisAtom<std::complex<float>>;
extern template class BasisAtom<std::complex<double>>;
