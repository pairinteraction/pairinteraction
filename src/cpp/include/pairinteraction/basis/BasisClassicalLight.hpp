#pragma once

#include "pairinteraction/basis/Basis.hpp"
#include "pairinteraction/ket/KetClassicalLight.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>

namespace pairinteraction {
template <typename Scalar>
class BasisClassicalLightCreator;

template <typename Scalar>
class BasisClassicalLight;

template <typename Scalar>
struct traits::CrtpTraits<BasisClassicalLight<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
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
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    friend class BasisClassicalLightCreator<Scalar>;
    struct Private {};

public:
    using Type = BasisClassicalLight<Scalar>;
    using ket_t = typename traits::CrtpTraits<Type>::ket_t;
    using ketvec_t = typename traits::CrtpTraits<Type>::ketvec_t;

    BasisClassicalLight(Private /*unused*/, ketvec_t &&kets);

    Eigen::VectorX<Scalar> get_matrix_elements(std::shared_ptr<const ket_t> /*ket*/,
                                               OperatorType /*type*/, int /*q*/) const override;
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
    get_matrix_elements(std::shared_ptr<const Type> /*other*/, OperatorType /*type*/,
                        int /*q*/) const override;
};

extern template class BasisClassicalLight<double>;
extern template class BasisClassicalLight<std::complex<double>>;
} // namespace pairinteraction
