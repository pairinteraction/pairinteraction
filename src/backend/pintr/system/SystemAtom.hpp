#pragma once

#include "pintr/system/System.hpp"
#include "pintr/utils/traits.hpp"

#include <array>
#include <memory>
#include <vector>

template <typename Scalar>
class OperatorAtom;

template <typename Scalar>
class BasisAtom;

template <typename Real>
class KetAtom;

// Specialize CrtpTraits for SystemAtom
template <typename T>
class SystemAtom;

template <typename Scalar>
struct traits::CrtpTraits<SystemAtom<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetAtom<real_t>;
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

    void set_electric_field(const std::array<real_t, 3> &field);
    void set_magnetic_field(const std::array<real_t, 3> &field);
    void enable_diamagnetism(bool enable);

private:
    std::array<Scalar, 3> electric_field_spherical{0, 0, 0};
    std::array<Scalar, 3> magnetic_field_spherical{0, 0, 0};
    bool diamagnetism_enabled{false};

    void construct_hamiltonian() const override;
};

extern template class SystemAtom<float>;
extern template class SystemAtom<double>;
extern template class SystemAtom<std::complex<float>>;
extern template class SystemAtom<std::complex<double>>;
