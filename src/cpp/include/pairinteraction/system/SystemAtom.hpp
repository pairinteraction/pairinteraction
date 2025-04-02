// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/system/System.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <array>
#include <memory>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
class OperatorAtom;

template <typename Scalar>
class BasisAtom;

class KetAtom;

template <typename T>
class SystemAtom;

template <typename Scalar>
struct traits::CrtpTraits<SystemAtom<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetAtom;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
    using basis_t = BasisAtom<scalar_t>;
    using operator_t = OperatorAtom<scalar_t>;
};

template <typename Scalar>
class SystemAtom : public System<SystemAtom<Scalar>> {
public:
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using Type = SystemAtom<Scalar>;
    using real_t = typename traits::CrtpTraits<Type>::real_t;
    using basis_t = typename traits::CrtpTraits<Type>::basis_t;

    SystemAtom(std::shared_ptr<const basis_t> basis);

    Type &set_electric_field(const std::array<real_t, 3> &field);
    Type &set_magnetic_field(const std::array<real_t, 3> &field);
    Type &set_diamagnetism_enabled(bool enable);

private:
    std::array<Scalar, 3> electric_field_spherical{};
    std::array<Scalar, 3> magnetic_field_spherical{};
    bool diamagnetism_enabled{false};

    void construct_hamiltonian() const override;
};

extern template class SystemAtom<double>;
extern template class SystemAtom<std::complex<double>>;
} // namespace pairinteraction
