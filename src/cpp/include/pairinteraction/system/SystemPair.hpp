// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/system/System.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <array>
#include <limits>
#include <memory>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
class OperatorPair;

template <typename Scalar>
class BasisPair;

template <typename Scalar>
class KetPair;

template <typename T>
class SystemPair;

template <typename Scalar>
class SystemAtom;

template <typename Scalar>
class BasisAtom;

template <typename Scalar>
class GreenTensor;

template <typename Scalar>
struct traits::CrtpTraits<SystemPair<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetPair<Scalar>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
    using basis_t = BasisPair<scalar_t>;
    using operator_t = OperatorPair<scalar_t>;
};

template <typename Scalar>
class SystemPair : public System<SystemPair<Scalar>> {
public:
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using Type = SystemPair<Scalar>;
    using real_t = typename traits::CrtpTraits<Type>::real_t;
    using basis_t = typename traits::CrtpTraits<Type>::basis_t;

    SystemPair(std::shared_ptr<const basis_t> basis);

    Type &set_interaction_order(int value);
    Type &set_distance_vector(const std::array<real_t, 3> &vector);
    Type &set_green_tensor(std::shared_ptr<const GreenTensor<Scalar>> &green_tensor);

private:
    int interaction_order{3};
    std::array<real_t, 3> distance_vector{0, 0, std::numeric_limits<real_t>::infinity()};
    std::shared_ptr<const GreenTensor<Scalar>> user_defined_green_tensor;

    void construct_hamiltonian() const override;
};

extern template class SystemPair<double>;
extern template class SystemPair<std::complex<double>>;
} // namespace pairinteraction
