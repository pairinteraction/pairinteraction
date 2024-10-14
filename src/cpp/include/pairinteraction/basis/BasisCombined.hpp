#pragma once

#include "pairinteraction/basis/Basis.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>

namespace pairinteraction {
class Database;

template <typename Real>
class KetCombined;

template <typename Scalar>
class BasisCombined;

template <typename Scalar>
class BasisAtom;

template <typename Real>
class KetAtom;

template <typename Scalar>
struct traits::CrtpTraits<BasisCombined<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetCombined<real_t>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
};

template <typename Scalar>
class BasisCombined : public Basis<BasisCombined<Scalar>> {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using Type = BasisCombined<Scalar>;
    using real_t = typename traits::CrtpTraits<Type>::real_t;
    using ketvec_t = typename traits::CrtpTraits<Type>::ketvec_t;

public:
    BasisCombined(ketvec_t &&kets);
};

extern template class BasisCombined<float>;
extern template class BasisCombined<double>;
extern template class BasisCombined<std::complex<float>>;
extern template class BasisCombined<std::complex<double>>;
} // namespace pairinteraction
