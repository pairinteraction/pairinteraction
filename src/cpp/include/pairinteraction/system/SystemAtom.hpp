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
    Type &set_ion_distance_vector(const std::array<real_t, 3> &vector);
    Type &set_ion_charge(real_t charge);
    Type &set_ion_interaction_order(int value);

private:
    std::array<real_t, 3> electric_field{0, 0, 0};
    std::array<real_t, 3> magnetic_field{0, 0, 0};
    bool diamagnetism_enabled{false};
    std::array<real_t, 3> ion_distance_vector{0, 0, std::numeric_limits<real_t>::infinity()};
    real_t ion_charge{1};
    int ion_interaction_order{3};

    void construct_hamiltonian() const override;
};

extern template class SystemAtom<double>;
extern template class SystemAtom<std::complex<double>>;
} // namespace pairinteraction
