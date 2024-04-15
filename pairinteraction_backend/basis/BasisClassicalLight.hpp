#pragma once

#include "basis/Basis.hpp"
#include "ket/KetClassicalLight.hpp"
#include "utils/Traits.hpp"

#include <complex>
#include <limits>

template <typename T>
class BasisClassicalLightCreator;

template <typename Scalar>
class BasisClassicalLight;

template <typename Scalar>
struct Traits::BasisTraits<BasisClassicalLight<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename Traits::NumTraits<Scalar>::real_t;
    using ket_t = KetClassicalLight<real_t>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
};

/**
 * @class BasisAtom
 *
 * @brief Class for creating a basis of atomic kets.
 *
 * @tparam T Complex number type.
 */
template <typename Scalar>
class BasisClassicalLight : public Basis<BasisClassicalLight<Scalar>> {
public:
    using Type = BasisClassicalLight<Scalar>;
    using ketvec_t = typename Traits::BasisTraits<Type>::ketvec_t;

private:
    friend class BasisClassicalLightCreator<Scalar>;
    BasisClassicalLight(ketvec_t &&kets);
};

extern template class BasisClassicalLight<float>;
extern template class BasisClassicalLight<double>;
extern template class BasisClassicalLight<std::complex<float>>;
extern template class BasisClassicalLight<std::complex<double>>;
