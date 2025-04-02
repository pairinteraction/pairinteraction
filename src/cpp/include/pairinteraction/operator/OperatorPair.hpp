// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/operator/Operator.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <memory>
#include <vector>

namespace pairinteraction {
enum class OperatorType;

template <typename Scalar>
class BasisPair;

template <typename Scalar>
class KetPair;

template <typename T>
class OperatorPair;

template <typename Scalar>
struct traits::CrtpTraits<OperatorPair<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetPair<Scalar>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
    using basis_t = BasisPair<scalar_t>;
};

template <typename Scalar>
class OperatorPair : public Operator<OperatorPair<Scalar>> {
public:
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using Type = OperatorPair<Scalar>;
    using basis_t = typename traits::CrtpTraits<Type>::basis_t;

    OperatorPair(std::shared_ptr<const basis_t> basis);
    OperatorPair(std::shared_ptr<const basis_t> basis, OperatorType type);
};

extern template class OperatorPair<double>;
extern template class OperatorPair<std::complex<double>>;
} // namespace pairinteraction
