#pragma once

#include "pairinteraction/operator/Operator.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <memory>
#include <vector>

namespace pairinteraction {
enum class OperatorType;

template <typename Scalar>
class BasisCombined;

template <typename Real>
class KetCombined;

template <typename T>
class OperatorCombined;

template <typename Scalar>
struct traits::CrtpTraits<OperatorCombined<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetCombined<Scalar>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
    using basis_t = BasisCombined<scalar_t>;
};

template <typename Scalar>
class OperatorCombined : public Operator<OperatorCombined<Scalar>> {
public:
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using Type = OperatorCombined<Scalar>;
    using basis_t = typename traits::CrtpTraits<Type>::basis_t;

    OperatorCombined(std::shared_ptr<const basis_t> basis);
    OperatorCombined(std::shared_ptr<const basis_t> basis, OperatorType type);
};

extern template class OperatorCombined<float>;
extern template class OperatorCombined<double>;
extern template class OperatorCombined<std::complex<float>>;
extern template class OperatorCombined<std::complex<double>>;
} // namespace pairinteraction
